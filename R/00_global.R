source("R/bigFasterize.R")
source('R/splitRast.R')
source('R/dtOptimal.R')

Require::Require(c('terra','sf','DBI','curl'))

options(timeout = 5000)

gp <- list(rawDataPath = 'data/raw',
           cleanDataPath = 'data/clean',
           outputPath = 'output',
           ## desired project coordinate reference system
           dcrs = sf::st_crs(6623)$wkt,
           ## desired output raster resolution
           dres = 100)

if(!dir.exists(gp$rawDataPath)) dir.create(gp$rawDataPath, recursive = T)
if(!dir.exists(gp$cleanDataPath)) dir.create(gp$cleanDataPath, recursive = T)
if(!dir.exists(gp$outputPath)) dir.create(gp$outputPath, recursive = T)
