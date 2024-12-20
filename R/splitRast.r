#' r <- rast('data/input/QC/rbox.tif')
#' nx <- 10
#' ny <- 10
#' buffer <- c(0, 0)
#' 
#' 
#' 
#' #' Split and re-merge \code{RasterLayer}(s)
#' #'
#' #' \code{splitRaster} divides up a raster into an arbitrary number of pieces (tiles).
#' #' Split rasters can be recombined using \code{do.call(merge, y)} or \code{mergeRaster(y)},
#' #' where \code{y <- splitRaster(x)}.
#' #'
#' #' This function is parallel-aware, using the same mechanism as used in \pkg{raster}.
#' #' Specifically, if you start a cluster using \code{\link{beginCluster}},
#' #' then this function will automatically use that cluster.
#' #' It is always a good idea to stop the cluster when finished, using \code{\link{endCluster}}.
#' #'
#' #' @param r       The raster to be split.
#' #'
#' #' @param nx      The number of tiles to make along the x-axis.
#' #'
#' #' @param ny      The number of tiles to make along the y-axis.
#' #'
#' #' @param buffer  Numeric vector of length 2 giving the size of the buffer along the x and y axes.
#' #'                If these values less than or equal to \code{1} are used, this
#' #'                is interpreted as the number of pixels (cells) to use as a buffer.
#' #'                Values between \code{0} and \code{1} are interpreted as proportions
#' #'                of the number of pixels in each tile (rounded up to an integer value).
#' #'                Default is \code{c(0, 0)}, which means no buffer.
#' #'
#' #' @param path    Character specifying the directory to which the split tiles will be saved.
#' #'                If missing, the function will write to memory.
#' #' @param cl      A cluster object. Optional. This would generally be created using
#' #'                parallel::makeCluster or equivalent. This is an alternative way, instead
#' #'                of \code{beginCluster()}, to use parallelism for this function, allowing for
#' #'                more control over cluster use.
#' #' @param rType   Data type of the split rasters. Defaults to FLT4S.
#' #'
#' #' @param fExt    file extension (e.g., \code{".grd"} or \code{".tif"}) specifying the file format.
#' #'
#' #' @return \code{splitRaster} returns a list (length \code{nx*ny}) of cropped raster tiles.
#' #'
#' #' @seealso \code{\link{do.call}}, \code{\link[raster]{merge}}.
#' #'
#' #' @author Alex Chubaty and Yong Luo
#' #' @export
#' #' @importFrom magrittr %>%
#' #' @importFrom parallel clusterApplyLB
#' #' @importFrom raster crop crs<- extent getCluster returnCluster writeRaster
#' #' @importFrom raster xmax xmin xres ymax ymin yres
#' #' @importFrom Require checkPath
#' #' @rdname splitRast
#' #'
#' #' @example inst/examples/example_splitRast.R
#' #'
#' setGeneric(
#'   "splitRast",
#'   function(r, nx = 1, ny = 1, buffer = c(0, 0), path = NA, cl, rType = "FLT4S", fExt = ".grd") {
#'     standardGeneric("splitRast")
#'   })
#' 
#' #' @export
#' #' @rdname splitRast
#' setMethod(
#'   "splitRast",
#'   signature = signature(r = "SpatRaster"),
#'   definition = 
    
splitRast <- function(r, nx, ny, buffer=c(0,0), path, cl, rType, fExt,
                      out='tiles') {
  
    if (!is.numeric(nx) | !is.numeric(ny) | !is.numeric(buffer)) {
      stop("nx, ny, and buffer must be numeric")
    }
  
    if (!is.integer(nx)) nx <- as.integer(nx)
    if (!is.integer(ny)) ny <- as.integer(ny)
    if (is.integer(buffer)) buffer <- as.numeric(buffer)
    
    if (!is.na(path)) {
      Require::checkPath(path, create = TRUE)
    }
    
    if (missing(cl)) {
      cl <- tryCatch(getCluster(), error = function(e) NULL)
      on.exit(if (!is.null(cl)) returnCluster(), add = TRUE)
    }
    
    if (length(buffer) > 2) {
      warning("buffer contains more than 2 elements - only the first two will be used.")
      buffer <- buffer[1:2]
    } else if (length(buffer) == 1) {
      buffer <- c(buffer, buffer)
    }
    if (buffer[1] < 1) {
      buffer[1] <- ceiling((buffer[1] * (xmax(r) - xmin(r)) / nx) / xres(r)) # nolint
    }
    if (buffer[2] < 1) {
      buffer[2] <- ceiling((buffer[2] * (ymax(r) - ymin(r)) / ny) / yres(r)) # nolint
    }
    
    # ext <- extent(r)
    # extents <- vector("list", length = nx * ny)
    # n <- 1L
    # for (i in seq_len(nx) - 1L) {
    #   for (j in seq_len(ny) - 1L) {
    #     x0 <- ext@xmin + i * ((ext@xmax - ext@xmin) / nx) - buffer[1] * xres(r) # nolint
    #     x1 <- ext@xmin + (i + 1L) * ((ext@xmax - ext@xmin) / nx) + buffer[1] * xres(r) # nolint
    #     y0 <- ext@ymin + j * ((ext@ymax - ext@ymin) / ny) - buffer[2] * yres(r) # nolint
    #     y1 <- ext@ymin + (j + 1L) * ((ext@ymax - ext@ymin) / ny) + buffer[2] * yres(r) # nolint
    #     extents[[n]] <- extent(x0, x1, y0, y1)
    #     n <- n + 1L
    #   }
    # }
    
    rext <- ext(r)
    extents <- vector("list", length = nx * ny)
    n <- 1L
    for (i in seq_len(nx) - 1L) {
      for (j in seq_len(ny) - 1L) {
        x0 <- xmin(rext) + i * ((xmax(rext) - xmin(rext)) / nx) - buffer[1] * xres(r) # nolint
        x1 <- xmin(rext) + (i + 1L) * ((xmax(rext) - xmin(rext)) / nx) + buffer[1] * xres(r) # nolint
        y0 <- ymin(rext) + j * ((ymax(rext) - ymin(rext)) / ny) - buffer[2] * yres(r) # nolint
        y1 <- ymin(rext) + (j + 1L) * ((ymax(rext) - ymin(rext)) / ny) + buffer[2] * yres(r) # nolint
        extents[[n]] <- ext(x0, x1, y0, y1)
        n <- n + 1L
      }
    }
    
    if(out=='boxes') return(extents)
    
    tiles <- if (!is.null(cl)) {
      clusterApplyLB(cl = cl, x = seq_along(extents), fun = cropRast, 
                     e = extents, r = r, path = path, rType = rType, 
                     fExt = fExt)
    } else {
      lapply(X = seq_along(extents), FUN = cropRast, e = extents, r = r, path = path, rType = rType, fExt = fExt)
    }
    
    return(tiles)
    
  }

cropRast <- function(i, e, r, path, rType, fExt) {
  ri <- terra::crop(r, e[[i]], datatype = rType)
  crs(ri) <- crs(r)
  if (is.na(path)) {
    return(ri)
  } else {
    filename <- file.path(path, paste0(names(r), "_tile", i, fExt))
    terra::writeRaster(ri, filename, overwrite = TRUE, datatype = rType)
    return(rast(filename))
  }
}
