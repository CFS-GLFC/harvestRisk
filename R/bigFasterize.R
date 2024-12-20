##############################################################################################################################################
##############################################################################################################################################
##
## Date Created: Oct 25, 2021
## Auteur: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
## Name of script : "bigFasterize.R"
## Description : R program to handle and rasterize huge geodatabase layers (e.g. provincial FRI) by
##               disassembling, rasterizing and subsequently reassembling based on the field of interest
##
##############################################################################################################################################
##############################################################################################################################################

bigFasterize <- function(src_data = NULL,
                         layer = NULL,
                         trast = NULL,
                         nx = 10, ny = 10,
                         rfield = NULL,
                         outfile = NULL,
                         # fx = NULL,
                         # jquery = NULL,
                         ncores = parallel::detectCores() - 1,
                         verbose.arg = FALSE,
                         convert2integer = F,
                         ow = TRUE,
                         dtype = NULL) {

  if(is.null(outfile)) stop('Must assign a valid output file name.')

  suppressMessages(Require::Require(c('stringr','terra','sf','dplyr','gdalUtilities','filesstrings','DBI','RSQLite','parallel','parallelly')))

  ## Validate inputs (e.g. src_data, layer, rfield)
  message("Validating inputs...")
  if(!file.exists(src_data)) stop('src_data cannot be found')
  layer <- match.arg(layer, st_read(src_data, query = 'select tbl_name from sqlite_master where type = "table"', quiet = TRUE)$tbl_name)
  if(!is.null(rfield)) rfield <- match.arg(rfield,
                                           (st_read(src_data, quiet = TRUE,
                                                    query = paste0('PRAGMA table_info(', layer, ')')) %>% 
                                              filter(pk != 1) %>% 
                                              pull(name)))

  ## Create temp directory
  udir <- file.path(tempfile('rtiles'))
  if(dir.exists(udir)) file.remove(udir)
  dir.create(udir)

  #################################
  ## Establish EPSG of input feature layer
  iepsg <- (st_read(src_data, query = paste0('select * from ', layer, ' limit 1'), quiet = TRUE) %>% st_crs)$wkt

  ## Establish desired EPSG
  depsg <- crs(rast(trast))

  ## Establish desired spatial resolution
  dres <- res(rast(trast))

  if(iepsg == depsg) {
    message(str_c('EPSG of input layer and trast identical. No reprojection necessary.'), '\n')
  } else {
    message('Reprojection necessary.')
  }

  ## Establish DB connection
  dbcon <- dbConnect(RSQLite::SQLite(), src_data)
  dbExecute(dbcon, 'select load_extension("mod_spatialite")')
  
  ## Fetch data type of target column
  dtype <- dbGetQuery(dbcon, paste0('PRAGMA table_info (', layer, ')')) %>%
    filter(name == rfield) %>% pull(type)
  
  if(!is.null(rfield)) {
    
    if(is.null(layer)) stop('argument "layer" must be provided when "rfield" is non-null')
    
    ## Create reference table if rfield is not numeric or integer
    if(str_detect(dtype, 'TEXT')) {
      
      message('Field to rasterize is a STRING. Creating lookup table for proxy variable (INTEGER)...')

      ## create lookup table for corresponding unique integer values
      rtab <- data.frame(label = dbGetQuery(dbcon,
                                            paste0('select distinct ', rfield, ' from ', layer)) %>%
                           pull() %>% sort) %>%
        mutate(ID = c(0, 1:(length(label)-1))) %>%
        relocate(ID, .before=label)

      rtab <- bind_rows(rtab, data.frame(label = NA, ID = max(rtab$ID) + 1))
      
    } else {
      
      message('Field to rasterize is of type ', dtype, '. Determining optimal datatype for raster processing...')
      rtype <- dtOptimal(dbGetQuery(dbcon, paste0('select min(', rfield, ') as minval, max(', rfield, ') as maxval from ', layer)))
      
    } 
    
  } else {
    
    DBI::dbExecute(dbcon, paste('alter table', layer, 'add column col1 integer default 1'))
    rfield <- 'col1'
    dtype <- 'INTEGER'
    
  }
  
  dbDisconnect(dbcon)

  ## Divide (temporarily reprojected) target raster area into nx * ny smaller parts
  message('Subdiving target raster area into ', nx*ny, ' equal tiles...')

  source('R/splitRast.r')
  xfile <- paste0(udir, '/extents.rds')
  saveRDS(lapply(splitRast(r=rast(trast), nx=nx, ny=ny, path=udir,
                           rType='INT1U', fExt = ".tif", out='boxes'), as.vector),
          file=xfile)

  ###################################################################################
  ## Multi-core deploy rasterization of sub-extents (tiles)
  message(paste('Deploying clip + rasterize operation', '(n =', nx*ny, 'tiles) across', ncores, 'threads...'))
  
  cl <- makeClusterPSOCK(ncores, default_packages = c('sf','stringr','gdalUtilities','terra','dplyr','DBI','RSQLite'))
  vlist <- c('udir','layer','rfield','verbose.arg','outfile','iepsg','depsg', 'ow','dres','src_data','xfile')

  if(!is.null(rfield)) {
    if(str_detect(dtype, 'TEXT')) vlist <- c(vlist, 'rtab', 'dtype') else vlist <- c(vlist, 'dtype', 'rtype')
  } else {
    vlist <- c(vlist)
  }

  clusterExport(cl, varlist = vlist, envir = environment())

  ptime <- system.time({

    parLapply(cl, seq.int(from = 1, to = nx*ny), function(i) {

      ## Create temporary folder
      tdir <- file.path(udir, str_c('rtile_', i))
      if(file.exists(tdir)) unlink(tdir, recursive=T)
      dir.create(tdir)

      ## retrieve destination clip bounding box
      xbox <- readRDS(xfile)[[i]]

      ## transform clip area to input coordinates
      tex <- st_as_sfc(st_bbox(xbox))
      st_crs(tex) <- depsg
      tex <- st_bbox(st_transform(tex, iepsg))

      ## Clip mega vector to sub-extent where conditions met
      message('clipping large vector to sub-tile...')

      ## clip large vector feature layer to subtile using GDAL
      gdalUtilities::ogr2ogr(src_datasource_name = src_data,
                             layer = layer,
                             dst_datasource_name = paste0(tdir, '/rclip_', i, '.gpkg'),
                             spat = tex,
                             spat_srs = iepsg,
                             select = rfield,
                             t_srs = depsg,
                             overwrite = ow)

      #######################################
      ## Rasterize based on field of interest (optional)
      message('rasterizing...')

      if(file.exists(str_c(udir, '\\fRAST_', i, '.tif'))) {
        unlink(str_c(udir, '\\fRAST_', i, '.tif'))
      }

      dbcon <- dbConnect(RSQLite::SQLite(), paste0(tdir, '/rclip_', i, '.gpkg'))
      dbExecute(dbcon, 'select load_extension("mod_spatialite")')

      if(dbGetQuery(dbcon, paste('select count(', rfield, ') from ', layer)) > 0) {

          r <- rasterize(x = vect(paste0(tdir, '/rclip_', i, '.gpkg')),
                         y = rast(ext(xbox), resolution = dres, crs = depsg),
                         field = rfield)

          if(str_detect(dtype, 'TEXT')) {

            cortab <- levels(r)[[1]]
            names(cortab) <- c('ID','label')
            cortab <- mutate(cortab, OID = rtab$ID[match(cortab$label, pull(rtab, label))])

            classify(x = as.numeric(r),
                     rcl = select(cortab, ID, OID),
                     filename = paste0(tdir, '\\fRAST_', i, '.tif'),
                     overwrite = T)

          } else {
            
            terra::writeRaster(r, paste0(tdir, '\\fRAST_', i, '.tif'))
            
          }
          
      } else {

        init(rast(ext(xbox), resolution = dres, crs = depsg), NA,
             filename = paste0(tdir, '\\fRAST_', i, '.tif'),
             datatype = rtype,
             overwrite = ow)

      }

      dbDisconnect(dbcon)

      ## Ensure proper crop to cutline (avoid gdal errors)
      crop(rast(paste0(tdir, '\\fRAST_', i, '.tif')),
           rast(ext(xbox), resolution = dres, crs = depsg),
           filename = paste0(udir, '\\fRAST_', i, '.tif'), overwrite = T)

      ## Delete temporary folder
      unlink(tdir, recursive = T)

    })

    stopCluster(cl)

    ###################################################################
    ## Merge subcomponent rasters into one giant (provincial) raster
    rfiles <- list.files(udir, pattern='fRAST_', full.names=TRUE)

    ## create copy of rbox to populate with values
    if(file.exists(outfile)) {
      if(!ow) {
        stop('outfile already exists and ow = F')
      } else {
        unlink(outfile, recursive=T)
      }
    }

    message("Preparing aggregated output raster...")

    ## Build virtual raster
    gdalUtilities::gdalbuildvrt(gdalfile = rfiles,
                                output.vrt = paste0(udir, '/', 'tempmix.vrt'))
    ## Write to raster mosaic
    gdalUtilities::gdalwarp(srcfile = paste0(udir, '/', 'tempmix.vrt'),
                            dstfile = paste0(udir, '/', 'temprast.tif'),
                            te = ext(trast)[c(1,3,2,4)],
                            tr = res(trast),
                            r = 'near')

    orast <- rast(paste0(udir, '/', 'temprast.tif'))

    if(!is.null(rfield)) {

      if(str_detect(dtype, 'TEXT')) {

        message('Assigning labels to raster values...')

        ## set category labels
        names(rtab)[2] <- rfield
        set.cats(setMinMax(orast), value = rtab)

      }
    }

    if(convert2integer) {

      ## convert
      orast <- classify(orast,
                        rcl = levels(orast)[[1]] %>%
                          mutate(!!rfield := as.integer(UQ(rlang::sym(rfield)))))

    }

    ## write raster with corrected bbox to file
    writeRaster(orast, outfile, overwrite = ow)

    ## remove temp file
    unlink(udir, recursive = TRUE)

  })

  message(str_c('Processing complete. Elapsed time = ', round(as.numeric(ptime[3]) / 60, digits=1), ' minutes.'))
  
  return(invisible(NULL))

}
