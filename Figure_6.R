###############################################################################################
##################
####### Project: Bernadouze | Remote sensing on vegetation
####### Title: Figure 6
####### Date: 07/07/2026
############################################################################################

rm(list=ls(all=TRUE))

###set working directory ----

setwd("WORKING_DIRECTORY_LOCATION/Data")

###packages ----

library(dplyr)
library(tidyverse)
library(vegan)
library(terra)
library(ggplot2)
library(ggrepel)
library(stats)
library(magrittr)

###data import ----

env_mean <- read.csv("Sentinel-2 data/S2_env_year.csv") #environmental parameters for each year
S2_2017 <-terra::rast("Sebtinel-2 data/S2_201706_RST_theta.tif") #Sentinel-2 species occurrences raster (theta) from June 2017

###Figure 6 : RDA between predicted plant species occurrences and year ----

S2_theta.df <- as.data.frame(S2_2017)
WTD_select <- filter(env_mean,year==2017)
S2_theta.df$year <- 2017

t_start <- 2018
t_final <- 2024
stack_S2_theta <- c()

for (t in t_start:t_final){
  rast2add <- terra::rast(paste0("Sentinel-2 data/S2_",t,"06_RST_theta.tif"))
  
  stack_S2_theta <- append(stack_S2_theta, rast2add)
  
  rast2add.df <- as.data.frame(rast2add)
  rast2add.df$year <- t
  S2_theta.df <- rbind(S2_theta.df, rast2add.df)
}

theta_env.df <- cbind(S2_theta.df, env_mean[,-c(1,2)])
theta_env.df$year <- as.factor(theta_env.df$year)

rda_theta <- rda(theta_env.df[,1:33] ~ year, theta_env.df)

S2_theta.df %<>%
  mutate(realID=as.character(rownames(theta_env.df)))

fortif_mod <- fortify(rda_theta)
fortif_mod$label <- gsub("\\."," ", fortif_mod$label)

fortif_centroids <- subset(fortif_mod, score=="centroids")
fortif_centroids$year <- as.factor(seq(2017,2024,by=1))

fortif_species <- subset(fortif_mod, score=="species")

#calculate mean rainfall and WTD sd by year

env_year <- aggregate(theta_env.df[,c(36,40)], list(theta_env.df$year), mean)
colnames(env_year) <- c("year","Rainf","WTD_sd")

centroid_scores <- fortif_centroids[c("year","RDA1","RDA2")]
centroid_scores <- left_join(centroid_scores,env_year,by="year")

env.w <- hclust(dist(scale(centroid_scores[,c(4,5)])), "ward.D2") #make clusters of years

gr <- cutree(env.w, k=4) #cut dendogram to yield 4 groups
grl <- levels(factor(gr))
centroid_scores$group <- gr

centroides_plot <- centroid_scores
species_plot <- fortif_species

#define a hydraulic profile for each year cluster

group_means <- aggregate(centroid_scores[, c("Rainf", "WTD_sd")],
                         by = list(Group = centroid_scores$group),
                         FUN = mean)

#classify clusters from driest to wettest
group_means <- group_means[order(group_means$Rainf), ]

group_means$profile <- factor(c("Driest", "Moderately dry", "Moderately wet", "Wettest"),
                              levels = c("Driest", "Moderately dry", "Moderately wet", "Wettest"))

group_profile_mapping <- data.frame(
  group = 1:length(grl),
  profile = group_means$profile[match(1:length(grl), group_means$Group)]
)

centroides_plot$profile <- group_profile_mapping$profile[centroid_scores$group]

#plot
perc <- round(100*(summary(rda_theta)$cont$importance[2, 1:2]), 2)

ggplot() +
  geom_point(data = centroides_plot,
             aes(x = RDA1, y = RDA2, color = profile), shape = "square", size = 4) +
  geom_text_repel(data = centroides_plot,
                  aes(x = RDA1, y = RDA2, label = year),
                  size = 4, color = "black") +
  geom_text_repel(data = species_plot,
                  aes(x = RDA1, y = RDA2, label = label),
                  size = 3.5, color = "darkgrey", max.overlaps = 30) +
  labs(title = "RDA - Year clusters and effect on predicted plant species communities",
       x = "RDA1", y = "RDA2") +
  scale_color_brewer(palette = "Set1", name = "Hydraulic profile") +
  xlab(paste0("RDA1 (", perc[1], "%)")) +
  ylab(paste0("RDA2 (", perc[2], "%)")) +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))