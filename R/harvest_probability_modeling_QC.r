##############################################################################################################################################
##############################################################################################################################################
##
## Date Created: August 09, 2022
## Auteur: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
## Name of script : "harvest_probability_modeling_QC.r"
## Description : R script to extract variables serving as a proxy of merchantable volume from Quebec's données écoforestières.
##
##############################################################################################################################################
##############################################################################################################################################

source('R/00_global.R')

# ## 1.1) Download sous-domaines bioclimatiques (QC) ----
# download.file('https://diffusion.mffp.gouv.qc.ca/Diffusion/DonneeGratuite/Foret/DONNEES_FOR_ECO_SUD/Classification_ecologique/CLASSI_ECO_QC_GDB.zip',
#               destfile = paste0(gp$rawDataPath, '/CLASSI_ECO_QC_GDB.zip'))
# 
# unzip(file.path(gp$rawDataPath, '/CLASSI_ECO_QC_GDB.zip'), exdir = gp$rawDataPath)
# file.remove(file.path(gp$rawDataPath, 'CLASSI_ECO_QC_GDB.zip'))

## 1.1) Create hard links to millgate cost rasters ----
x <- list.files('D:/GitHub/millGateCosts/output', pattern = 'volume_wtd_rnetDist2mills', full.names = T)
y <- list.files('D:/GitHub/millGateCosts/output', pattern = 'volume_wtd_rnetDist2mills')

dir.create(paste0(gp$cleanDataPath, '/millGateCosts'))
lapply(1:length(x), function(i) file.link(from = x[i], to = paste0(gp$cleanDataPath, '/millGateCosts/', y[i])))

### 1.1.2) Create mosaic rasters
gdalUtilities::gdalwarp(srcfile = ,
                        dstfile = )
  
## 1.2) Download current QC Ecoforestry data ----
download.file('https://diffusion.mffp.gouv.qc.ca/Diffusion/DonneeGratuite/Foret/DONNEES_FOR_ECO_SUD/Cartes_ecoforestieres_perturbations/CARTE_ECO_MAJ_PROV_GPKG.zip',
              destfile = paste0(gp$rawDataPath, '/CARTE_ECO_MAJ_PROV_GPKG.zip'))

unzip(file.path(gp$rawDataPath, 'CARTE_ECO_MAJ_PROV_GPKG.zip'), exdir = gp$rawDataPath)
file.remove(file.path(gp$rawDataPath, 'CARTE_ECO_MAJ_PROV_GPKG.zip'))

## 1.3) Download forest harvest interventions (QC MFFP) ----
download.file('https://diffusion.mffp.gouv.qc.ca/Diffusion/DonneeGratuite/Foret/INTERVENTIONS_FORESTIERES/Recolte_et_reboisement/INTERV_FORES_PROV_GPKG.zip',
              destfile = paste0(gp$rawDataPath, '/INTERV_FORES_PROV_GPKG.zip'))

unzip(file.path(gp$rawDataPath, 'INTERV_FORES_PROV_GPKG.zip'), exdir = gp$rawDataPath)
file.remove(file.path(gp$rawDataPath, 'INTERV_FORES_PROV_GPKG.zip'))

## 1.4) Download QC managed forest boundaries ----
download.file('https://diffusion.mffp.gouv.qc.ca/Diffusion/DonneeGratuite/Foret/LIM_TERRITOIRE_FOREST_PUBLIC/UA/UA_GPKG.zip',
              destfile = paste0(gp$rawDataPath, '/UA_GPKG.zip'))

unzip(file.path(paste0(gp$rawDataPath, '/UA_GPKG.zip')), files = 'UA_GPKG/STF_UA.gpkg', exdir = gp$rawDataPath)
file.rename(file.path(gp$rawDataPath, 'UA_GPKG/STF_UA.gpkg'), paste0(gp$rawDataPath, '/STF_UA.gpkg'))
file.remove(file.path(gp$rawDataPath, 'UA_GPKG'), recursive = T)

### 1.4.1) Rasterize managed forest boundaries (raster template)
ua <- st_read(file.path(gp$rawDataPath, 'STF_UA.gpkg'), quiet = T)
terra::rasterize(ua,
                 rast(ext(st_bbox(ua)), res = gp$dres, crs = gp$dcrs),
                 filename = paste0(gp$rawDataPath, '/rtemplate.tif'))

## 1.5) Download HRDEM
curl::curl_download('https://datacube-prod-data-public.s3.ca-central-1.amazonaws.com/store/elevation/mrdem/mrdem-30/mrdem-30-dsm.tif',
              destfile = paste0(gp$rawDataPath, '/mrdem-30-dsm.tif'),
              quiet = FALSE)

##################################################################
## 2.1) Extract an_origine from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data = file.path(gp$dataPath, 'CARTE_ECO_MAJ_PROV.gpkg'), layer = 'pee_maj', 
             trast = rast(file.path(gp$dataPath, 'STF_UA.tif')),
             rfield = 'an_origine', nx = 10, ny = 10,
             outfile = paste0(gp$outputPath, '/an_origine.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
             ncores = parallelly::availableCores() - 1,
             convert2integer = T, 
             ow = T, verbose.arg = F)

## 2.2) Extract origine from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data = file.path(gp$dataPath, 'CARTE_ECO_MAJ_PROV.gpkg'), layer = 'pee_maj', 
             trast = rast(file.path(gp$dataPath, 'STF_UA.tif')),
             rfield = 'origine', nx = 10, ny = 10,
             outfile = paste0(gp$outputPath, '/origine.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
             ncores = parallelly::availableCores() - 1,
             ow = T, verbose.arg = TRUE)

## 2.3) Extract age class from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data = file.path(gp$dataPath, 'CARTE_ECO_MAJ_PROV.gpkg'), layer = 'pee_maj', 
             trast = rast(file.path(gp$dataPath, 'STF_UA.tif')),
             rfield = 'cl_age', nx = 10, ny = 10,
             outfile = paste0(gp$outputPath, '/cl_age.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
             ncores = parallelly::availableCores() - 1,
             ow = T, verbose.arg = TRUE)

## 2.4) Extract density class from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data = file.path(gp$dataPath, 'CARTE_ECO_MAJ_PROV.gpkg'), layer = 'pee_maj', 
             trast = rast(file.path(gp$dataPath, 'STF_UA.tif')),
             rfield = 'cl_dens', nx = 10, ny = 10,
             outfile = paste0(gp$outputPath, '/cl_dens.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
             ncores = parallelly::availableCores() - 1,
             ow = T, verbose.arg = F)

## 2.5) Create and subsequently rasterize unique FRI polygon identifiers from CARTE_ECO_MAJ_PROV ----
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(gp$dataPath, 'CARTE_ECO_MAJ_PROV.gpkg'))
dbExecute(con, 'select load_extension("mod_spatialite")')
dbExecute(con, 'ALTER TABLE PEE_MAJ ADD COLUMN uniqueID INTEGER;')
dbExecute(con, 'UPDATE PEE_MAJ SET uniqueID = ROWID;')
dbDisconnect(con)

bigFasterize(src_data = file.path(gp$dataPath, 'CARTE_ECO_MAJ_PROV.gpkg'), layer = 'pee_maj', 
             trast = rast(file.path(gp$dataPath, 'STF_UA.tif')),
             rfield = 'uniqueID', nx = 10, ny = 10,
             outfile = paste0(gp$outputPath, '/uniqueID.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
             ncores = parallelly::availableCores() - 1,
             ow = T, verbose.arg = F)

## 2.6) Process 30x30m DEM ----

### 2.6.1) Crop & project DEM to study area projection & extent ----
system.time({
  gdalUtilities::gdalwarp(srcfile = file.path(gp$dataPath, 'mrdem-30-dsm.tif'),
                          dstfile = paste0(gp$dataPath, '/dem_crop.tif'),
                          te = st_bbox(ext(rast(file.path(gp$dataPath, 'STF_UA.tif')))),
                          te_srs = crs(rast(file.path(gp$dataPath, 'STF_UA.tif'))),
                          t_srs = crs(rast(file.path(gp$dataPath, 'STF_UA.tif'))),
                          tr = res(rast(file.path(gp$dataPath, 'mrdem-30-dsm.tif'))),
                          overwrite = T)
}) # 72 sec

### 2.6.2) Compute surface roughness ----
system.time({
  roughness <- terrain(x = rast(paste0(gp$dataPath, '/dem_crop.tif')),
                       v = "roughness",
                       neighbors = 8,
                       filename = paste0(gp$dataPath, '/roughness.tif'),
                       overwrite = T)
})

### 2.6.3) Compute slope ----
system.time({
  slope <- terrain(x = rast(paste0(gp$dataPath, '/dem_crop.tif')),
                   v = "slope",
                   neighbors = 8,
                   filename = paste0(gp$dataPath, '/slope.tif'),
                   overwrite = T)
}) #110 sec

### 2.6.4) Resample DEM, roughness & slope to desired extent & resolution ----
for(i in c('dem_crop.tif','slope.tif','roughness.tif')) {

 message('i = ', i)
  
  system.time({
    gdalUtilities::gdalwarp(srcfile = file.path(gp$dataPath, i),
                            dstfile = paste0(gp$dataPath, '/', i, '_', gp$dres, 'm.tif'),
                            te = st_bbox(ext(rast(file.path(gp$dataPath, 'STF_UA.tif')))),
                            tr = res(rast(file.path(gp$dataPath, 'STF_UA.tif'))),
                            r = 'bilinear',
                            srcnodata = 'None')
  })
  
}
  
gdalUtilities::gdalbuildvrt(gdalfile = c(dem = paste0(gp$dataPath, '/dem_crop.tif'),
                                         roughness = paste0(gp$dataPath, '/roughness_', gp$dres, 'm.tif'),
                                         slope = paste0(gp$dataPath, '/roughness_', gp$dres, 'm.tif')),
                            output.vrt = paste0(gp$dataPath, '/temp.vrt'),
                            separate = T)

gdalUtilities::gdalwarp(srcfile = file.path(gp$dataPath, 'temp.vrt'),
                        dstfile = paste0(gp$dataPath, '/dem_variables_', gp$dres, 'm.tif'))

file.remove(file.path(gp$dataPath, 'temp.vrt'))

## 2.7) Extract unique cutblock identifier from INTERV_FORES_PROV ----
rasterize(x = st_read(file.path(gp$dataPath, 'INTERV_FORES_PROV.gpkg'), 
                      query = 'select geom from INTERV_FORES_PROV where an_origine > 2015',
                      quiet = TRUE) %>%
            dplyr::mutate(cutblock_id = 1:nrow(.), .before = geom),
          y = rast(file.path(gp$dataPath, 'STF_UA.tif')),
          field = 'cutblock_id',
          filename = paste0(gp$outputPath, '/cutblock_id.INTERV_FORES_PROV.tif'))
