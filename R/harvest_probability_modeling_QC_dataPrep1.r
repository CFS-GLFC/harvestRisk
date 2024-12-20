##############################################################################################################################################
##############################################################################################################################################
##
## Date Created: August 09, 2022
## Auteur: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
## Name of script : "harvest_probability_modeling_QC_data_prep.r"
## Description : R script to extract candidate variables serving to explain observed
##               variation in harvesting patterns (QC ecoforestry databases).
##
##############################################################################################################################################
##############################################################################################################################################

setwd("D:/RHome")
source("D:/RHome/R/functions/BigFasterize.R")
library(terra)
library(sf)

## 1) Extract an_origine from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data='D:/RHome/data/raw/QC/DONNEES_ECOFORESTIERES/D5. 2015/CARTE_ECO_MAJ_PROV/CARTE_ECO_MAJ_PROV.gpkg',
  layer='pee_maj', 
  trast=rast('D:/RHome/data/input/QC/FIEqc/rRGA.tif'),
  rfield='an_origine', nx=10, ny=10,
  outfile='E:/TEIB/ECCC-S&T/data/an_origine.pee_maj.CARTE_ECO_MAJ_PROV.tif',
  ncores=parallel::detectCores() - 1,
  convert2integer=F, 
  ow=T, verbose.arg=F)

## 2) Extract origine from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data='D:/RHome/data/raw/QC/DONNEES_ECOFORESTIERES/CARTE_ECO_MAJ_PROV.gpkg',
             layer='pee_maj', trast='D:/RHome/data/input/QC/rbox.tif',
             rfield='origine', nx=10, ny=10,
             outfile='E:/TEIB/ECCC-S&T/data/origine.pee_maj.CARTE_ECO_MAJ_PROV.tif',
             ncores=parallel::detectCores() - 1,
             ow=T, verbose.arg=TRUE)

## 3) Extract age class from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data='D:/RHome/data/raw/QC/DONNEES_ECOFORESTIERES/D5. 2015/CARTE_ECO_MAJ_PROV/CARTE_ECO_MAJ_PROV.gpkg',
             layer='pee_maj', trast='D:/RHome/data/input/QC/rbox.tif',
             rfield='cl_age', nx=10, ny=10,
             outfile='E:/TEIB/ECCC-S&T/data/cl_age.pee_maj.CARTE_ECO_MAJ_PROV.tif',
             ncores=parallel::detectCores() - 1,
             ow=T, verbose.arg=TRUE)

## 4) Extract density class from CARTE_ECO_MAJ_PROV ----
bigFasterize(src_data='D:/RHome/data/raw/QC/DONNEES_ECOFORESTIERES/D5. 2015/CARTE_ECO_MAJ_PROV/CARTE_ECO_MAJ_PROV.gpkg',
             layer='pee_maj', trast='D:/RHome/data/input/QC/rbox.tif',
             rfield='cl_dens', nx=10, ny=10,
             outfile='E:/TEIB/ECCC-S&T/data/cl_dens.pee_maj.CARTE_ECO_MAJ_PROV.tif',
             ncores=parallel::detectCores() - 1,
             ow=T, verbose.arg = F)

## 5) Extract unique cutblock identifier from INTERV_FORES_PROV ----
Require::Require(c('terra','sf'))
rasterize(x = st_read('D:/RHome/data/raw/QC/DONNEES_ECOFORESTIERES/Perturbations/INTERV_FORES_PROV.gpkg', 
                      query='select geom from INTERV_FORES_PROV where an_origine > 2015') %>%
            dplyr::mutate(cutblock_id = 1:nrow(.), .before = geom),
          y = rast('D:/RHome/data/input/QC/rbox.tif'),
          field = 'cutblock_id',
          filename = 'E:/TEIB/ECCC-S&T/data/cutblock_id.INTERV_FORES_PROV.tif')
