###############################################################################################
##################
####### Project: Bernadouze | Remote sensing on vegetation
####### Title: Figure 3
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
library(ggvegan)
library(ggrepel)
library(magrittr)

###data import ----

#plot data
perm <- read.csv("Plot data/UAV_ALL_DATA.csv") #all data from plots (location, ID, vegetation occurrence, VIs, etc)
perm_abu_veg <- perm[,c(1:7,64:119)] #select only plant occurrences from perm dataframe
perm_abu_veg.h <- decostand(perm_abu_veg, "hellinger")
perm_evt_VIs <- perm[,1:54] #select only VIs data from perm dataframe
perm_evt_VIs <- perm_evt_VIs %>% mutate_at(scale, .vars = vars(-colnames(perm_evt_VIs[c(1:7, 9, 11, 13, 15, 17)])))

#UAV validation data
R2_df_all <- read.csv("UAV_R2_df_cross_validation.csv") #R² results per species for all methods (cross validation and null models)
R2_LOPO <- read.csv("UAV_LOPO_results_random_splits_90_iterations.csv") #R² results per species for leave-one-out method
R2_LOPO$method <- "leave-one-out"

R2_LOPO_simple <- R2_LOPO[,-c(3,4)] #only keep R², species and method data
colnames(R2_LOPO_simple) <- colnames(R2_df_all)
R2_df_all <- rbind(R2_df_all, R2_LOPO_simple) #combine all methods in a dataframe

R2_df_all$name_species <- str_replace_all(R2_df_all$name_species, pattern="\\.", replacement=" ")
R2_df_all$name_species <- as.factor(R2_df_all$name_species)
R2_df_all$method <- as.factor(R2_df_all$method)

#random forest models results
iteration_summary <- read.csv("RF_mod_all_iterations.csv") #results of all iterations of random forest models per species

###Figure 3 (a) : RDA between vegetation communities and selected VIs ----

perm %<>%
  mutate(realID=as.character(11:37))

rda <- rda(perm_abu_veg.h ~ GCC_MEAN + BCC_MEAN + NDRE_MEAN + MEAN_BAND1_R + ARI_MEAN, perm_evt_VIs[,c(8,10,12,14,16,18:54)])

fortif_data <- fortify(rda)
fortif_data%<>%
  left_join(rename(perm[,c("Date","realID")],label=realID))
fortif_data$label <- sub("[.]", " ", fortif_data$label)

perc <- round(100*(summary(rda)$cont$importance[2, 1:2]), 2)

autoplot(rda, layer="biplot", arrow.col = "black", data=fortif_data, label=TRUE, title="First axes of the RDA model between the vegetation abundance matrix and the selected VIs")+
  expand_limits(x=c(-1.5,1.5)) +
  geom_point(data=subset(fortif_data, fortif_data$score == "sites"), aes(RDA1, RDA2, colour=Date)) +
  geom_text_repel(data=subset(fortif_data, fortif_data$score == "species"), aes(RDA1,RDA2,label=label), colour="darkgrey") +
  xlab(paste0("RAD1 (",perc[1], "%)")) +
  ylab(paste0("RAD2 (",perc[2], "%)")) +
  theme_classic() +
  scale_colour_manual(values = cols, limits=c("02/06/2023", "16/06/2023", "04/07/2023"))

###Figure 3 (b) : global R² coefficient for each validation method ----

dodge2 <- position_dodge(width = 0.5)
my_comparisons <- list(c("cross-validation", "leave-one-out"), c("cross-validation", "null model RA3"), c("cross-validation", "null model sim5"), c("leave-one-out", "null model RA3"), c("leave-one-out", "null model sim5"), c("null model RA3", "null model sim5"))
my_symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("***", "***", "**", "*", "ns"))

ggplot(data = R2_df_all, aes(x = as.factor(method), y = r2_species, fill = method))+
  geom_boxplot(width =0.4, position = dodge2)+
  stat_compare_means(comparisons = my_comparisons, symnum.args = my_symnum.args, method="t.test")+ #add pairwise comparisons p-value
  scale_fill_manual(values = c("#D01556FF", "#FFA400FF", "#50A860FF", "lightblue"))+
  labs(title="Global UAV R² for each validation method",
       x = "Method",
       y = "R²")+
  theme_classic()

###Figure 3 (c) : global performance metrics of random forest models ----

species_data <- iteration_summary %>%
  select(species, accuracy, tss, auc) %>%
  pivot_longer(
    cols = c(accuracy, tss, auc), 
    names_to = "Metric", 
    values_to = "Value"
  ) %>%
  mutate(Panel_Group = "By Species") #Column identifier for faceting

all_species_combined <- species_data %>%
  mutate(species = "All Species Combined",
    Panel_Group = "Overall Benchmark") #Column identifier for faceting

plot_data <- bind_rows(species_data, all_species_combined)
plot_data$Panel_Group <- factor(plot_data$Panel_Group, levels = c("By Species", "Overall Benchmark"))
plot_data$Metric <- factor(plot_data$Metric, levels = c("accuracy", "tss", "auc"))

plot_combined_species_data <- filter(plot_data, plot_data$Panel_Group == "Overall Benchmark")

dodge2 <- position_dodge(width = 0.5)
ggplot(plot_combined_species_data, aes(x = species, y = Value, fill=Metric)) + 
  geom_boxplot(width = 0.4, position = dodge2) +
  facet_wrap(.~ Metric, scales = "free") +
  ggtitle("Performance metrics from UAV Random forest model for all plant species") +
  xlab("species name") +
  ylab("Percentage / Score") +
  theme_classic() +
  scale_fill_manual(values = c("#D01556FF", "#FFA400FF", "#50A860FF")) +
  theme(legend.position="none")