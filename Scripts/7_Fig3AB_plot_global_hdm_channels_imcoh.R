
### This script plots global results across channel densities for connectivity
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
# install.packages("WRS", repos="http://R-Forge.R-project.org") #nicht installiert!
install.packages("psych")
install.packages("Hmisc")
install.packages("Rcmdr")
install.packages("splines")
install.packages("gridExtra")
install.packages("grid")	
install.packages("ggpubr")
install.packages("cowplot")
# install.packages("lme4")
# install.packages("lme4",
# repos=c("http://lme4.r-forge.r-project.org/repos",
# getOption("repos")[["CRAN"]]))
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
# library(WRS)
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


setwd("~/Documents/Projects/Channels/R/files")

data = readMat('global_analysis_allchannels_coh_img.mat')
a = array(unlist(data$output.metric), dim = c(96,5))
all = as.data.frame(a)

names(all)[1] = 'head model'
names(all)[2] = 'channeldensity'
names(all)[3] = 'frequency'
names(all)[4] = 'cohend'
names(all)[5] = 'p'

# change/reorder levels of variables
levels(all$channeldensity)[6] = '256'
levels(all$channeldensity)[7] = '19 (classic 10-20)'
levels(all$channeldensity)[8] = '25 (IFCN 10-20)'

# convert from factor to numberic values
all$cohend = as.numeric(levels(all$cohend))[all$cohend]
all$p = as.numeric(levels(all$p))[all$p]

# add variable indicating signficiance
all$p_sign = rep(' ', length(all$p))
all$p_sign[all$p < 0.05] = '*'

# change order of freqs
all$frequency_n = factor(all$frequency, levels=c('delta','theta','alpha','beta1', 'beta2', 'gamma'))
all$channeldensity_n = factor(all$channeldensity, levels=c("256","192","128","64","48","32", "25 (IFCN 10-20)", "19 (classic 10-20)"))

# split dataset according to headmodel
all_indiv = all[all$`head model` == 'individual', ]
all_canon = all[all$`head model` == 'canonical', ]

# loop over frequency and plot cohen's d for each resolution

# use wrap for each frequency
data = all_indiv

g1 = ggplot(data, aes(x=channeldensity_n, y=cohend, label=p_sign)) +
  geom_segment( aes(x=channeldensity_n, xend=channeldensity, y=0, yend=cohend), color="grey") +
  geom_point( color="orange", size=2) +
  scale_y_continuous(limits = c(0, 0.95)) +
  geom_text(size=6, vjust = 0, nudge_y = 0.01) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor=element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 55, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.text= element_text(size=14),
    axis.title = element_text(size=14)
  ) +
  xlab("individual headmodel") +
  ylab("cohen d (connectivity)") + 
  facet_grid(.~frequency_n, switch="both") 
  
  
data = all_canon

g2 =  ggplot(data, aes(x=channeldensity_n, y=cohend, label=p_sign)) +
  geom_segment( aes(x=channeldensity_n, xend=channeldensity, y=0, yend=cohend), color="grey") +
  geom_point( color="aquamarine3", size=2) +
  scale_y_continuous(limits = c(0, 0.95)) +
  geom_text(size=6, vjust = 0, nudge_y = 0.01) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor=element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 55, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.text= element_text(size=14),
    axis.title = element_text(size=14)
  ) +
  xlab("canoncial headmodel") +
  ylab("cohen d (connectivity)") + 
  facet_grid(.~frequency_n, switch="both") 


plot = ggarrange(g1, g2, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

ggsave(file="global_d_channels_imcoh_grey.png", plot=plot, dpi = 300, limitsize = TRUE, width = 12, height = 9)	


