###############################################################################################
##################
####### Project: Bernadouze | Remote sensing on vegetation
####### Title: Figure 5
####### Date: 07/07/2026
############################################################################################

rm(list=ls(all=TRUE))

###set working directory ----

setwd("WORKING_DIRECTORY_LOCATION/Data")

###packages ----

library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)
library(terra)
library(parallel)
library(doParallel)
library(sf)
library(stats)
library(corrplot)

###data import ----

#data to process PCOA on predicted species occurrences
UAV_theta <- terra::rast("UAV data/UAV_2023_RST_theta.tif")
theta_17 <- terra::rast("Sentinel-2 data/S2_201706_RST_theta.tif")

#data to process heatmaps of species versus PCOA axes correlations
load("UAV data/UAV_2023_30000_samples_pcoa_results.RData") #load pcoa_result matrix
all_df <- as.data.frame(terra::values(UAV_theta))

###Figure 5 (a) and (c) : plant communities maps (PCoA axes scores) ----

#(a) UAV ----

theta_df <- as.data.frame(terra::values(UAV_theta))

#remove NAs rows and sample some pixels of the raster
na_rows <- apply(theta_df, 1, function(x) any(is.na(x))) #identify NAs id rows
valid_coords_df <- theta_df[!na_rows, ] #extract rows without NAs

gc()
valid_coords_samp <- valid_coords_df[sample(nrow(valid_coords_df), 30000),] #sampling for pcoa calculation
samp_indices <- as.numeric(rownames(valid_coords_samp)) #extract sampled pixels coordinates 
samp_coords <- terra::xyFromCell(theta, samp_indices)

#Bray-Curtis ditance matrix applied on valid pixel sampled
gc()
cl <- makeCluster(18)
registerDoParallel(cl)
dist_bray <- vegdist(valid_coords_samp, method = "bray")
stopCluster(cl)

#PCOA calculation on Bray-Curtis matrix
set.seed(123) #for reproducibility
gc()
pcoa_result <- wcmdscale(dist_bray, k = 3) #PCOA with 3 axis

#preparation of the PCOA data and spatial data (global map extent)
raster_chunk <- theta[[1]] #extract 1rst layer to obtain the map extent
theta_coarse <- terra::aggregate(raster_chunk, fact = 20)  #resample map resolution by a factor 20 to make calculation easier

pcoa_coords <- as.data.frame(pcoa_result)
samp_indices <- as.numeric(rownames(pcoa_coords))
samp_coords <- terra::xyFromCell(theta, samp_indices)

pcoa_data <- data.frame(
  x = samp_coords[, 1],
  y = samp_coords[, 2],
  axis1 = pcoa_coords[, 1],
  axis2 = pcoa_coords[, 2],
  axis3 = pcoa_coords[, 3]
)

projcrs_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
pcoa_sf <- st_as_sf(x = pcoa_data,                         
                    coords = c("x", "y"),
                    crs = projcrs_wgs84)

bbox <- terra::ext(theta_coarse) #raster extent
xmin <- as.numeric(bbox[1]) #extract bbox coordinates
ymin <- as.numeric(bbox[3])
xmax <- as.numeric(bbox[2])
ymax <- as.numeric(bbox[4])

bbox_sf <- sf::st_as_sfc(
  sf::st_bbox(
    c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
    crs = 4326  #coordinate system WGS84
  )
)

#extract theta_coarse raster parameters
res <- res(theta_coarse)  #resolution (x, y)
ext <- ext(theta_coarse)  #extent
crs_theta <- crs(theta_coarse)  #CRS

#define tiles size (number of pixels per tile)
tile_size <- c(40, 40)

#create tiles grid aligned with pixels
x_coords <- seq(ext[1], ext[2], by = res[1] * tile_size[1])
y_coords <- seq(ext[3], ext[4], by = res[2] * tile_size[2])

#construct polygons corresponding to the tiles
tile_coords <- map(x_coords, function(x) {
  map(y_coords, function(y) {
    matrix(c(
      x, y,
      x + res[1] * tile_size[1], y,
      x + res[1] * tile_size[1], y + res[2] * tile_size[2],
      x, y + res[2] * tile_size[2],
      x, y
    ), ncol = 2, byrow = TRUE)
  })
}) |> unlist(recursive = FALSE)

tile_polygons <- map(tile_coords, ~ st_polygon(list(.x))) |>
  st_sfc() |>
  st_sf() |>
  st_set_crs(crs_theta)

#cut theta_coarse following the tiles pattern
#theta_coarse <- mask(theta_coarse, theta_coarse)  #remove NAs
chunk_rasters <- lapply(tile_polygons$st_sfc.map.tile_coords...st_polygon.list..x...., function(poly) {
  tile <- crop(theta_coarse, vect(poly)) |> mask(vect(poly))
  return(tile)
})

#filter and keep tiles non empty
tile_status <- lapply(chunk_rasters, function(raster) {
  stats <- global(raster, fun = function(x) sum(!is.na(x)))
  data.frame(non_na_count = stats[[1]], is_empty = stats[[1]] == 0)
}) |> do.call(rbind, args = _)

non_empty_tiles <- tile_polygons[tile_status$is_empty == FALSE, ]
non_empty_chunk_rasters <- chunk_rasters[tile_status$is_empty == FALSE]

#visualize the grid
theta_sf <- st_as_sf(as.polygons(theta_coarse)) |> st_set_crs(crs_theta)
ggplot() +
  geom_sf(data = theta_sf, fill = "transparent", color = "red", linewidth = 1) +
  geom_sf(data = non_empty_tiles, fill = "green", alpha = 0.3) +
  labs(title = "Tuiles non vides (vert) alignées sur theta_coarse (rouge)") +
  theme_minimal()

#create a function to apply the iwd interpolation to each tile
idw_for_axis_and_chunk <- function(axis_name, theta_coarse) {
  formula <- as.formula(paste(axis_name, "~ 1"))
  idw_model <- gstat(data = pcoa_sf, formula = formula, set = list(idp = 5), nmax = 10)
  prediction <- predict(idw_model, newdata = st_as_stars(theta_coarse))
  return(prediction)
}

#test the idw function on a single tile 
test_result <- idw_for_axis_and_chunk("axis1", non_empty_chunk_rasters[[1]])
plot(test_result)

#apply the idw_for_axis_and_chunk() function to all the tiles
cl <- makeCluster(18)
registerDoParallel(cl)
axes <- c("axis1", "axis2", "axis3")
results <- list() # List to stock all the results

gc()
for (i in 1:length(non_empty_chunk_rasters)) {
  chunk_raster <- non_empty_chunk_rasters[[i]]
  
  chunk_results <- list(
    axis1 = idw_for_axis_and_chunk("axis1", chunk_raster),
    axis2 = idw_for_axis_and_chunk("axis2", chunk_raster),
    axis3 = idw_for_axis_and_chunk("axis3", chunk_raster)
  )
  
  names(chunk_results) <- axes
  results[[i]] <- chunk_results
  gc()
}
names(results) <- paste0("tile_", seq_along(non_empty_chunk_rasters))

#function to assemble the tiles from an axis
assemble_axis_tiles <- function(axis_name, results_list, theta_coarse_ref) {
  #extract all tiles for the specified axis and convert to SpatRaster
  axis_tiles <- lapply(results_list, function(tile_result) {
    stars_object <- tile_result[[axis_name]] #stars object
    spatraster <- rast(stars_object) #conversion into SpatRaster object
    return(spatraster)
  })
  #merge the tiles one by one
  merged_raster <- axis_tiles[[1]]  # Start with the 1rst tile
  for (i in 2:length(axis_tiles)) {
    merged_raster <- merge(merged_raster, axis_tiles[[i]])
  }
  #refocus on the extent of theta_coarse_ref
  merged_raster <- resample(merged_raster, theta_coarse_ref, method = "bilinear")
  merged_raster <- crop(merged_raster, ext(theta_coarse_ref))
  return(merged_raster)
}

#merge all tiles rasters for each axis with the assemble_axis_tiles() function
global_rasters <- list()
for (axis in c("axis1", "axis2", "axis3")) {
  print(paste("Assemblage de l'axe :", axis))
  global_rasters[[axis]] <- assemble_axis_tiles(axis, results, theta_coarse)
}

#save global rasters for each axis
dir.create("DIRECTORY_NAME", showWarnings = FALSE)

for (axis in names(global_rasters)) {
  writeRaster(
    global_rasters[[axis]],
    file.path(output_dir, paste0("UAV_2023_", axis, "_global.tif")),
    filetype = "GTiff",
    overwrite = TRUE
  )
}

#normalize each axis and create a RGB raster
global_rasters_axis1 <- terra::rast("UAV_2023_axis1_global.tif")
global_rasters_axis2 <- terra::rast("UAV_2023_axis2_global.tif")
global_rasters_axis3 <- terra::rast("UAV_2023_axis3_global.tif")

normalize <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  (x - x_min) / (x_max - x_min) * 255
}

pcoa_axis1_norm <- normalize(terra::values(global_rasters_axis1[[1]]))
pcoa_axis2_norm <- normalize(terra::values(global_rasters_axis2[[1]]))
pcoa_axis3_norm <- normalize(terra::values(global_rasters_axis3[[1]]))

#create RGB raster
rgb_raster <- terra::rast(
  nlyr = 3,
  nrows = nrow(global_rasters_axis1),
  ncols = ncol(global_rasters_axis1),
  ext = terra::ext(global_rasters_axis1),
  crs = terra::crs(global_rasters_axis1)
)

terra::values(rgb_raster)[, 1] <- pcoa_axis1_norm
terra::values(rgb_raster)[, 2] <- pcoa_axis2_norm
terra::values(rgb_raster)[, 3] <- pcoa_axis3_norm

#plot RGB raster
terra::plotRGB(rgb_raster,
  r = 1, g = 2, b = 3,
  stretch = "hist",
  main = "Plant communities map (PCOA scores) from UAV-derived JSDM predictions")

#save raster RGB for later
terra::writeRaster(rgb_raster,
                   filename = paste0("UAV_2023_communities_pcoa_axis.tif"),
                   filetype = "GTiff",
                   gdal = c("COMPRESS=LZW", "PREDICTOR=2"),
                   overwrite = TRUE)

#(c) Sentinel-2 ----

normalize <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  (x - x_min) / (x_max - x_min) * 255
}

#set theta df 2017
theta_df_17 <- as.data.frame(terra::values(theta_17))
na_rows_17 <- apply(theta_df_17, 1, function(x) any(is.na(x)))
valid_coords_df_17 <- theta_df_17[!na_rows_17,]
valid_coords_df_17$npx <- rownames(valid_coords_df_17)
valid_coords_df_17$year <- 2017

#Select same pixels for all years

nstart <- 2018
nend <- 2025
all_df <- valid_coords_df_17

for (year in nstart:nend){
  theta <- terra::rast(paste0("S2_",year,"06_RST_theta_rast.tif"))
  theta_df <- as.data.frame(terra::values(theta))
  na_rows <- apply(theta_df, 1, function(x) any(is.na(x)))
  valid_coords_df <- theta_df[!na_rows, ]
  valid_coords_df$npx <- rownames(valid_coords_df)
  valid_coords_df$year <- year
  all_df <- bind_rows(all_df, valid_coords_df)
}

dist_bray <- vegdist(all_df[,-c(34,35)], method = "bray")

set.seed(123) #for reproctubility
pcoa_result <- wcmdscale(dist_bray, k = 3)

pcoa_axes <- as.data.frame(pcoa_result)
pcoa_axes$npx <- all_df$npx
pcoa_axes$year <- all_df$year

nstart <- 2017
nend <- 2025

for (year in nstart:nend){
  theta <- terra::rast(paste0("S2_",year,"06_RST_theta_rast.tif"))
  theta_df <- as.data.frame(terra::values(theta))
  current_year <- as.integer(year)
  #select rows for each year
  pcoa_axes_select <- filter(pcoa_axes, pcoa_axes$year==current_year)
  print(dim(pcoa_axes_select))
  #give back npx for names rows
  rownames(pcoa_axes_select) <- pcoa_axes_select$npx
  pcoa_axes_select <- pcoa_axes_select[,-c(4,5)]
  #process raster axes
  pcoa_axis1 <- rep(NA, nrow(theta_df))
  pcoa_axis2 <- rep(NA, nrow(theta_df))
  pcoa_axis3 <- rep(NA, nrow(theta_df))
  
  pcoa_axis1[!na_rows] <- pcoa_axes_select[, 1]
  pcoa_axis2[!na_rows] <- pcoa_axes_select[, 2]
  pcoa_axis3[!na_rows] <- pcoa_axes_select[, 3]
  
  pcoa_raster_axis1 <- theta
  pcoa_raster_axis2 <- theta
  pcoa_raster_axis3 <- theta
  
  terra::values(pcoa_raster_axis1) <- pcoa_axis1
  terra::values(pcoa_raster_axis2) <- pcoa_axis2
  terra::values(pcoa_raster_axis3) <- pcoa_axis3
  
  pcoa_axis1_norm <- normalize(pcoa_axis1)
  pcoa_axis2_norm <- normalize(pcoa_axis2)
  pcoa_axis3_norm <- normalize(pcoa_axis3)
  
  rgb_raster <- terra::rast(
    nlyr = 3,
    nrows = nrow(theta),
    ncols = ncol(theta),
    ext = terra::ext(theta),
    crs = terra::crs(theta))
  
  terra::values(rgb_raster) <- matrix(NA, ncol = 3, nrow = terra::ncell(theta))
  
  terra::values(rgb_raster)[, 1] <- pcoa_axis1_norm
  terra::values(rgb_raster)[, 2] <- pcoa_axis2_norm
  terra::values(rgb_raster)[, 3] <- pcoa_axis3_norm
  #save rgb raster
  terra::writeRaster(rgb_raster,
                     filename = paste0("S2_",year,"06_communities_pcoa_axes_all_years.tif"),
                     filetype = "GTiff",
                     gdal = c("COMPRESS=LZW", "PREDICTOR=2"),
                     overwrite = TRUE)
}

rgb_raster <- terra::rast("S2_202306_communities_pcoa_axes_all_years_to2025.tif")
terra::plotRGB(rgb_raster,
  r = 1, g = 2, b = 3,
  stretch = "hist",  #stretching color values based on a histogram
  main = "2023 plant communities map (PCOA scores) from Sentinel-2-derived JSDM predictions")

###Figure 5 (b) and (d) : heatmaps revealing the correlations between plant species and PCoA axes ----

pcoa_coords <- as.data.frame(scores(pcoa_result))
pcoa_coords$year <- all_df$year
pcoa_coords <- filter(pcoa_coords,year==2023)
pcoa_coords <- pcoa_coords[,-4]

#extract plant species data (all columns except the last ones)
species_data <- filter(all_df, year==2023)
species_data <- species_data[, -c(1,35,36)]
species_data <- drop_na(species_data)

#calculate the correlations between PCOA axes and predicted plant species occurrences
correlation_matrix <- cor(species_data,pcoa_coords)
correlation_matrix <- as.data.frame(correlation_matrix)
correlation_matrix <- drop_na(correlation_matrix)
correlation_matrix <- arrange(correlation_matrix)

#plot correlation matrix
rownames(correlation_matrix) <- correlation_matrix$X
correlation_matrix <- as.matrix(correlation_matrix[,2:4])
rownames(correlation_matrix) <- gsub("\\."," ", rownames(correlation_matrix))

heatmap(correlation_matrix, margins=c(5,8), col = COL2(diverging = "RdYlBu", n = 1000))
colorlegend(xlim=c(17,20), ylim=c(0,33), COL2(diverging = "RdYlBu", n = 1000), c(seq(-1,1,0.5)), align="l", vertical=TRUE, addlabels=TRUE)
