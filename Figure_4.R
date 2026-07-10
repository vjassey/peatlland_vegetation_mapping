###############################################################################################
##################
####### Project: Bernadouze | Remote sensing on vegetation
####### Title: Figure 4
####### Date: 07/07/2026
############################################################################################

rm(list=ls(all=TRUE))

###set working directory ----

setwd("WORKING_DIRECTORY_LOCATION/Data")

###packages ----

library(dplyr)
library(tidyverse)
library(ggplot2)
library(terra)

###data import ----

#species occurrences raster (theta)
theta <- terra::rast("UAV data/UAV_2023_RST_theta.tif")

#species occurrence raster low threshold for 89% confidence interval calculation
theta_low <- terra::rast("UAV data/UAV_2023_RST_low_0.89_theta.tif")

#species occurrence raster high threshold for 89% confidence interval calculation
theta_high <- terra::rast("UAV data/UAV_2023_RST_high_0.89_theta.tif")

###Figure 4 (a) and (c) : plant species richness maps ----

species_richness <- terra::app(theta, sum)
names(species_richness) <- "species_richness"

raster_df <- as.data.frame(species_richness, xy = TRUE, na.rm = TRUE)

png(file = "pred_richness_map.png", #directory you want to save the file in
    width = 12000, #width of the plot in inches
    height = 10000, #height of the plot in inches
    res = 1000) #resolution of the plot

pred_map <- ggplot(raster_df, aes(x = x, y = y, fill = species_richness)) +
  geom_tile() +  #use geom_tile() to plot raster data
  scale_fill_viridis_c() +  #color scale for the values
  theme_classic() +
  coord_equal() +  #maintain aspect ratio
  theme(legend.position = "right") +
  labs(title = "Plant species richness map", 
       x = "Longitude", y = "Latitude", fill = "Plant species richness")
print(pred_map)

dev.off()

###Figure 4 (b) and (d) : plant species richness confidence interval maps with a 89% threshold ----

species_richness_low <- terra::app(theta_low, sum)
species_richness_high <- terra::app(theta_high, sum)

CI_rich <- species_richness_high - species_richness_low
names(CI_rich) <- "species_richness"

raster_df <- as.data.frame(CI_rich, xy = TRUE, na.rm = TRUE)

#plot
jpeg(file = "richness_credibility_int_89_percent.jpeg",
     width = 1200,
     height = 800,
     quality = 100)

pred_map <- ggplot(raster_df, aes(x = x, y = y, fill = species_richness)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic() +
  coord_equal() +
  theme(legend.position = "right") +
  labs(title = "Plant species richness 89% credibility interval", 
       x = "Longitude", y = "Latitude", fill = "Credibility interval (89%) of plant species richness")
print(pred_map)

dev.off()
