###############################################################################################
##################
####### Project: Bernadouze | Remote sensing on vegetation
####### Title: predictive maps generated from JSDM presence-absence model
####### Date: 07/07/2026
############################################################################################

rm(list=ls(all=TRUE))

###set working directory ----

setwd("WORKING_DIRECTORY_LOCATION/Data")

###packages ----

library(dplyr)
library(tidyverse)
library(terra)

###data import ----

#JSDM parameters
load("JSDM data/mod_JSDM_4chains_PA.RData")
load("JSDM data/JSDM_param_outputs.RData")
load("JSDM data/JSDM_theta_outputs.RData")

param_sites <- sitecoefs.mean
param_species <- spcoefs.mean
probit_theta_latent <- probit_theta_latent.mean
theta_latent <- theta_latent.mean

#site borders
map_extent <- terra::rast("UAV data/Raw data/BAND1_R_BERNADOUZE_202306_CUTFINAL.tif") #a raster to obtain the site extent
map_extent <- terra::project(map_extent,"+proj=longlat")

borders <- terra::as.polygons(map_extent > -Inf) #borders of the site as shp
borders <- terra::project(borders,"+proj=longlat")
ext_cut <- ext(1.422, 1.4245, 42.801, 42.8033) #optionnal : an extent can also be defined manually

#xy coordinates of plots
xy <- terra::vect("Plot data/coords_perm.shp", crs=crs(borders))

#VIs maps (predictive variables)
GCC_map <- terra::rast("UAV data/UAV_GCC_2023_06.tif")
GCC_map <- terra::project(GCC_map,"+proj=longlat")
GCC_map <- terra::crop(GCC_map, ext_cut, snap="near")

BCC_map <- terra::rast("UAV data/UAV_BCC_2023_06.tif")
BCC_map <- terra::project(BCC_map,"+proj=longlat")
BCC_map <- terra::crop(BCC_map, ext_cut, snap="near")

NDRE_map <- terra::rast("UAV data/UAV_NDRE_2023_06.tif")
NDRE_map <- terra::project(NDRE_map,"+proj=longlat")
NDRE_map <- terra::crop(NDRE_map, ext_cut, snap="near")

ARI_map <- terra::rast("UAV data/UAV_ARI_2023_06.tif")
ARI_map <- terra::project(ARI_map,"+proj=longlat")
ARI_map <- terra::crop(ARI_map, ext_cut, snap="near")

MEAN_BAND1_R <- terra::rast("UAV data/Raw data/BAND1_R_BERNADOUZE_202306_CUTFINAL.tif") #red band of multispectral image 
MEAN_BAND1_R <- terra::project(MEAN_BAND1_R,"+proj=longlat")
MEAN_BAND1_R <- terra::crop(MEAN_BAND1_R, ext_cut, snap="near")

scaled_clim_var <- c(GCC_map,BCC_map,NDRE_map,MEAN_BAND1_R,ARI_map,)
names(scaled_clim_var) <- c("GCC","BCC","NDRE","BAND1_R","ARI")

###interpolation of spatial variables ----

#interpolate spatial parameters (alpha, W1 and W2) at site scale using the ArcGIS PRO function "Spline with tension"
#other methods can be used for interpolation

###processing the results of interpolation ----

#process RST alpha raster
alpha_rst <- terra::rast("UAV data/UAV_RST_alpha.tif") #interpolated raster of alpha site parameter 

#raster restricted to site borders
alpha_rst <- terra::project(alpha_rst, terra::crs(borders))
alpha_rst  <- terra::mask(alpha_rst, borders)

#center interpolated site effect 
alpha_rst_centered <- terra::app(alpha_rst, fun=scale, scale=FALSE)
names(alpha_rst_centered) <- names(alpha_rst)

alpha_rst_centered_xy <- terra::extract(alpha_rst_centered, xy)[,"UAV_RST_alpha"]
plot(alpha_rst_centered_xy, param_sites[1:9,]$alphas,
     xlab="alpha interpolated by RST",
     ylab="alpha estimated by JSDM",
     main="Random site effect")
abline(a=0, b=1, col='red')
terra::writeRaster(alpha_rst_centered, "RST_alpha_centered.tif",
                   gdal=c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE)

#process RST W1 raster
W1_rst <- terra::rast("UAV data/UAV_RST_W1.tif") #interpolated raster of W1 site parameter 

#raster restricted to site borders
W1_rst <- terra::project(W1_rst, terra::crs(borders))
W1_rst  <- terra::mask(W1_rst, borders)

#center interpolated W1 site effect 
W1_rst_centered <- terra::app(W1_rst, scale, scale=FALSE)
names(W1_rst_centered) <- names(W1_rst)

W1_rst_centered_xy <- terra::extract(W1_rst_centered, xy)[,"UAV_RST_W1"]
plot(W1_rst_centered_xy, param_sites[1:9,]$W1,
     xlab="W1 interpolated by RST",
     ylab="W1 estimated by JSDM",
     main="First latent axis")
abline(a=0, b=1, col='red')
terra::writeRaster(W1_rst_centered, "RST_W1_centered.tif",
                   gdal=c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE)

#process RST W2 raster
W2_rst <- terra::rast("UAV data/UAV_RST_W2.tif") #interpolated raster of W2 site parameter 

#raster restricted to site borders
W2_rst <- terra::project(W2_rst, terra::crs(borders))
W2_rst  <- terra::mask(W2_rst, borders)

#center interpolated W2 site effect 
W2_rst_centered <- terra::app(W2_rst, scale, scale=FALSE)
names(W2_rst_centered) <- names(W2_rst)

W2_rst_centered_xy <- terra::extract(W2_rst_centered, xy)[,"UAV_RST_W2"]
plot(W2_rst_centered_xy, param_sites[1:9,]$W2,
     xlab="W2 interpolated by RST",
     ylab="W2 estimated by JSDM",
     main="Second latent axis")
abline(a=0, b=1, col='red')
terra::writeRaster(W2_rst_centered, "RST_W2_centered.tif",
                   gdal=c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE)

###predictive map of probabilities of occurrence (theta) -----

nsp <- nrow(param_species) #number of species

#same dimensions and crs for all rasters
rst_alpha_resamp <- terra::resample(rst_alpha, scaled_clim_var)
rst_W1_resamp <- terra::resample(rst_W1, scaled_clim_var)
rst_W2_resamp <- terra::resample(rst_W2, scaled_clim_var)

crs(scaled_clim_var) <- crs(rst_alpha_resamp)

#function to compute probit_theta at site scale 
predfun <- function(scaled_clim_var, params_species, rst_alpha, rst_W1, rst_W2, species.range){
  lambda_1 <- as.matrix(params_species[,"lambda_1"])
  lambda_2 <- as.matrix(params_species[,"lambda_2"])
  beta <- as.matrix(params_species[,1:6])
  
  print(species.range[1])
  print(species.range[2])
  #Xbeta_1
  np <- terra::nlyr(scaled_clim_var)
  Xbeta_1 <- terra::rast(ncols=dim(rst_alpha)[2], nrows=dim(rst_alpha)[1],
                         ext=terra::ext(rst_alpha), crs=terra::crs(rst_alpha),
                         resolution=terra::res(rst_alpha))
  terra::values(Xbeta_1) <- rep(beta[1,1][[1]], terra::ncell(Xbeta_1))
  for (p in 1:np) {
    Xbeta_1 <- Xbeta_1 + scaled_clim_var[[p]]*beta[1,p+1] 
  }
  #Wlambda_1
  Wlambda_1 <- rst_W1*lambda_1[1] + rst_W2*lambda_2[1]
  #probit_theta_1
  probit_theta_1 <- Xbeta_1 + Wlambda_1 + rst_alpha
  probit_theta <- probit_theta_1
  remove(list=c("probit_theta_1","Wlambda_1"))
  #Other species
  for (j in (species.range[1]+1):species.range[2]) {
    print(paste("Traitement de l'espèce ", j))
    if (j > nrow(beta)) {
      print(paste("Alerte : j =", j, " dépasse les limites de beta (nrow(beta) = ", nrow(beta), ")"))} 
    #Xbeta_j
    Xbeta_j <- Xbeta_1
    terra::values(Xbeta_j) <- rep(beta[j,1][[1]], terra::ncell(Xbeta_j))
    for (p in 1:np) {
      Xbeta_j <- Xbeta_j + scaled_clim_var[[p]]*beta[j,p+1] 
    }
    #Wlambda_j
    Wlambda_j <- rst_W1*lambda_1[j] + rst_W2*lambda_2[j]  
    #probit_theta_j
    probit_theta_j <- Xbeta_j + Wlambda_j + rst_alpha
    probit_theta <- c(probit_theta, probit_theta_j)
    remove(list=c("probit_theta_j", "Xbeta_j", "Wlambda_j"))
  }
  names(probit_theta) <- make.names(params_species$CI[species.range[1]:species.range[2]])
  return(probit_theta)
}

#compute theta
npart <- 1
first.species <- seq(1, nsp, by=floor(nsp/npart)+1)

for (n in 1:(npart)){
  probit_theta <- predfun(scaled_clim_var, param_species,
                          rst_alpha_resamp, rst_W1_resamp, rst_W2_resamp,
                          species.range=c(first.species[n],
                                          min(nsp,first.species[n]+floor(nsp/npart))))
  terra::tmpFiles(remove=TRUE)
  terra::writeRaster(probit_theta,
                     filename="UAV data/UAV_RST_probit_theta.tif",
                     filetype="GTiff",
                     gdal=c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE)
  #compute and save SpatRaster of probabilities of presence theta
  terra::app(probit_theta, pnorm, cores=2,
             filename="UAV data/UAV_RST_theta.tif", overwrite=TRUE, wopt=list(gdal=c("COMPRESS=LZW", "PREDICTOR=2"), filetype="GTiff"))
  remove(probit_theta)
}

###predictive map of species richness -----

theta <- terra::rast("UAV data/UAV_2023_RST_theta.tif")

species_richness <- terra::app(theta, sum)
names(species_richness) <- "species_richness"

terra::writeRaster(species_richness,
                   filename="UAV data/UAV_richness_species.tif", 
                   filetype="GTiff", gdal=c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE)
