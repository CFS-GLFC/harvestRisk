## incomplete and unreproducible script !!
Require::Require(c('terra','sf','dplyr'))
source('D:/RHome/R/00_1_TEIB_functions.R')
setwd('E:/TEIB/ECCC-S&T/data')

## Select roads in existence as of 2015 in QC, buffer and write to file
st_read("D:/RHome/data/raw/QC/AQ_ReseauPlus.gpkg", 
        query = "select geom, An_Classi from Reseau_routier where An_Classi NOT IN('2016','2017','2018','2019','2020','2021','2022')") %>%
  st_transform(., crs(rast('D:/RHome/data/input/QC/FIEqc/rRGA.tif'))) %>%
  st_buffer(., 25) %>%
  st_write(., 'D:/tmp/QC_roads_2015.gpkg', delete_layer = T)

## rasterize roads
vForm2(idsn = 'D:/tmp/QC_roads_2015.gpkg', ilayer = 'QC_roads_2015',
       orast = 'QC_road_network_2015.tif',
       rbox = 'D:/RHome/data/input/QC/FIEqc/rRGA.tif.tif')

## calculate distance to nearest road (2015)
distance(ifel(rast('QC_road_network_2015.tif') == 1, 1, NA), 
         filename = 'dist2nearestRoad_2015.tif',
         overwrite = T)

## calculate distance to nearest recent (<= 5 yrs) cutblock
distance(ifel(rast('an_origine.INTERV_FORES_PROV.tif') > 2010 &
                rast('an_origine.INTERV_FORES_PROV.tif') < 2016, 1, NA),
         unit = 'km', filename = 'dist2nearestRecentCut_5yr.tif')

## done for all pre-existing layers...
rForm2(irast = 'E:/TEIB/ECCC-S&T/data/an_origine.INTERV_FORES_PROV.tif',
       orast = 'G:/ECCC_ST/an_origine.INTERV_FORES_PROV.tif',
       rbox = 'D:/RHome/data/input/QC/FIEqc/rRGA.tif',
       r = 'near',
       ow = TRUE)

## remove recent fires from input layers
frast <- st_read("D:/RHome/data/raw/QC/FEUX_PROV_10/FEUX_PROV_10.gdb", quiet = T,
             query="select exercice from FEUX_PROV") %>%
  mutate(exercice = as.integer(exercice)) %>%
  filter(exercice > 2015) %>%
  st_transform(crs(rast('D:/RHome/data/input/QC/FIEqc/rRGA.tif'))) %>%
  rasterize(., rast('D:/RHome/data/input/QC/FIEqc/rRGA.tif'))

mask(
    ## binary harvested or not
  c(ifel(rast('an_origine.INTERV_FORES_PROV.tif') > 2015, 1, 0),
    ## unique cutblock identifier
    rast('cutblock_id.INTERV_FORES_PROV.tif'),
    ## year of disturbance
    rast('an_origine.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
    ## age class
    rast('cl_age.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
    ## density class
    rast('cl_dens.pee_maj.CARTE_ECO_MAJ_PROV.tif'),
    ## distance to nearest road in 2015
    rast('dist2nearestRoad_2015.tif'),
    ## distance to nearest recent cutblock (5 yrs)
    rast('dist2nearestRecentCut_5yr.tif'),
    ## volume-weighted distance to nearest mill
    rast('QC_accessCostDist_wtdxvol_mosaic.tif', lyr = 'rdist2nearestMill')
  ),
  mask = frast, inverse = T,  filename = 'rsandwich.tif', overwrite = T)

rsand <- rast('rsandwich.tif')
names(rsand) <- c('harvested','cutblock_id','stand_origin_year','age_class','density_class',
                  'dist2nearestRoad_2015','dist2nearestRecentCut_5yr','vol_wtd_dist2nearestMill')

## Cycle through every RGA and crop all inputs to multi-band RGA-specific raster sandwich
tenvec <- st_read('D:/RHome/data/input/QC/FIEqc/TDR_MOD.gpkg', 'RGA', quiet = T)

for(rga in tenvec$UNITE) {
  
  cat('RGA', rga, '\n')
  
  ## Crop raster to RGA
  crop(rsand, filter(tenvec, UNITE == rga), 
       mask = T, filename = paste0('rstack_RGA_', rga, '.tif'), overwrite = T)
  
  ## Save data.frame as RDS
  terra::as.data.frame(rast(paste0('rstack_RGA_', rga, '.tif')), xy=T) %>%
    filter(., complete.cases(.)) %>%
    mutate(harvested = factor(harvested)) %>%
    base::saveRDS(., paste0('HarvProb_model_inputData_RGA_', rga, '.rds'))
  
}

x <- lapply(tenvec$UNITE, function(rga) {
  readRDS(paste0('HarvProb_model_inputData_RGA_', rga, '.rds')) %>% head()
})


#### TO DO: ####
## Add Slope
## Add/calculate clumping/autocorrelation metric