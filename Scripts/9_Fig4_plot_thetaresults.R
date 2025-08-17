### This script plots effects of spatial downsampling (channel reductions) on
### group-results in the theta frequency bands 
### (correlation analyses between 256 channel map and others)
### Christina Stier, 2025

## R version 4.2.2 (2022-10-31)
## RStudio 2023.3.0.386 for macOS

rm(list = ls())

install.packages("corrplot")
install.packages("igraph")
install.packages("qgraph")
install.packages("car")
install.packages("compute.es")
install.packages("effects")
install.packages("compute.es")
install.packages("ggplot2")
install.packages("multcomp")
install.packages("pastecs")
#install.packages("WRS", repos="http://R-Forge.R-project.org") #nicht installiert!
install.packages("psych")
install.packages("Hmisc")
install.packages("Rcmdr")
install.packages("splines")
install.packages("gridExtra")
install.packages("grid")	
install.packages("ggpubr")
install.packages("cowplot")
install.packages("lme4")
install.packages("lme4",
                 repos=c("http://lme4.r-forge.r-project.org/repos",
                         getOption("repos")[["CRAN"]]))
install.packages("optimx")
install.packages("plyr")
install.packages("doBy")
install.packages("boot")
install.packages("lmPerm")
install.packages('R.matlab')
install.packages('abind')

library(qgraph)
library(corrplot)
library(igraph)
require(igraph)
library(compute.es)
library(effects)
library(ggplot2)
library(multcomp)
library(pastecs)
library(psych)
library(Hmisc)
#library(Rcmdr)
library(car)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(lme4)
library(optimx)
library(plyr)
library(doBy)
library(boot)
library(lmPerm)
library(R.matlab)
library(plyr)
library(abind)
library(reshape2)
library(tidyr)

setwd("~/Channels/R/files")

# get demographics for the main (HD-EEG) sample
dataset = as.data.frame(read.csv("~/sciebo/Channels/R/files/results_theta_vertices.csv", sep = ';'))
dataset$Metric = as.factor(dataset$Metric)
dataset$Head.model = as.factor(dataset$Head.model)
dataset$Channels = as.factor(dataset$Channels)
str(dataset)

#################################### split data
dataset$Head.model <- factor(dataset$Head.model, levels = c("individual", "canonical"))
conn = dataset[dataset$Metric == 'Connectivity (theta)',]
pow = dataset[dataset$Metric == 'Power (theta)',]


ggplot(conn, aes(x = Head.model, y = Channels, fill = r)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  geom_text(aes(label = sprintf("%.2f (p %s)", r, p_sign)), size = 7) +
  scale_fill_gradient2(
    low      = "grey",
    mid      = "white",
    high     = "orange",
    midpoint = 0.28,
    limits=c(-0.12,0.97),
    #midpoint = mean(pow$r),
    name     = "spatial correlation (r)"
  ) +
  theme_minimal() +
  labs(
    # title = "Spatial Correlation (Spearman’s r)\nacross Channel Densities",
    x     = "head model",
    y     = "channel density (theta connectivity)"
  ) +
  theme(text = element_text(size = 20),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y   = element_text(size = 20),
        #  plot.title    = element_text(hjust = 0.5)
)
ggsave(file="results_densities_theta_conn.png",  dpi = 300, limitsize = TRUE, width = 10, height = 9)	



# Plot power
ggplot(pow, aes(x = Head.model, y = Channels, fill = r)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  geom_text(aes(label = sprintf("%.2f (p %s)", r, p_sign)), size = 7) +
  scale_fill_gradient2(
    low      = "grey",
    mid      = "white",
    high     = "orange",
    midpoint = 0.28,
    limits=c(-0.12,0.97),
    #midpoint = mean(pow$r),
    name     = "spatial correlation (r)"
  ) +
  theme_minimal() +
  labs(
   # title = "Spatial Correlation (Spearman’s r)\nacross Channel Densities",
    x     = "head model",
    y     = "channel density (theta power)"
  ) +
  theme(text = element_text(size = 20),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y   = element_text(size = 20),
  #  plot.title    = element_text(hjust = 0.5)
  )

ggsave(file="results_densities_theta_power.png",  dpi = 300, limitsize = TRUE, width = 10, height = 9)	
