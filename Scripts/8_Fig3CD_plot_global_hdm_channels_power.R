
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

data = readMat('global_analysis_allchannels_power.mat')
a = array(unlist(data$output.metric), dim = c(96,5))
all = as.data.frame(a)

names(all)[1] = 'head model'
names(all)[2] = 'channeldensity'
names(all)[3] = 'frequency'
names(all)[4] = 'cohend'
names(all)[5] = 'p'

# change/reorder levels of variables
levels(all$channeldensity)[6] = '256'
levels(all$channeldensity)[7] = '19 (10-20 classic)'
levels(all$channeldensity)[8] = '25 (10-20 ext.)'

# convert from factor to numberic values
all$cohend = as.numeric(levels(all$cohend))[all$cohend]
all$p = as.numeric(levels(all$p))[all$p]

# change order of freqs
all$frequency_n = factor(all$frequency, levels=c('delta','theta','alpha','beta1', 'beta2', 'gamma'))
all$channeldensity_n = factor(all$channeldensity, levels=c("256","192","128","64","48","32", "25 (10-20 ext.)", "19 (10-20 classic)"))

# split dataset according to headmodel
all_indiv = all[all$`head model` == 'individual', ]
all_canon = all[all$`head model` == 'canonical', ]

# prepare loop over frequency for correction of the p-values for 8 channel sets
freq = c('delta', 'theta', 'alpha', 'beta1', 'beta2', 'gamma') 
pcorr_list_indiv = c()
pcorr_list_canon = c()

for (f in 1:6){
  all_indiv_f = all_indiv[all_indiv$frequency_n == freq[f], ]
  p_corr_indiv = as.data.frame(p.adjust(all_indiv_f$p, method = "fdr", n = length(all_indiv_f$p)))
  pcorr_list_indiv = rbind(pcorr_list_indiv, p_corr_indiv)
  
  all_canon_f = all_canon[all_canon$frequency_n == freq[f], ]
  p_corr_canon = as.data.frame(p.adjust(all_canon_f$p, method = "fdr", n = length(all_canon_f$p)))
  pcorr_list_canon = rbind(pcorr_list_canon, p_corr_canon)
}

names(pcorr_list_indiv)[1] = 'p_fdr'
names(pcorr_list_canon)[1] = 'p_fdr'

# sort according to frequency and apply fdr correction for 8 comparisons
all_indiv2 = all_indiv[order(all_indiv$frequency_n),]
all_indiv2 = cbind(all_indiv2, pcorr_list_indiv)

all_canon2 = all_canon[order(all_canon$frequency_n),]
all_canon2 = cbind(all_canon2, pcorr_list_canon)

# add variable indicating signficiance
all_indiv2$p_sign = rep(' ', length(all_indiv2$p))
all_indiv2$p_sign[all_indiv2$p_fdr < 0.05] = '*'

all_canon2$p_sign = rep(' ', length(all_canon2$p))
all_canon2$p_sign[all_canon2$p_fdr < 0.05] = '*'

# use wrap for each frequency
data = all_indiv2

g1 = ggplot(data, aes(x=channeldensity_n, y=cohend, label=p_sign)) +
  geom_segment( aes(x=channeldensity_n, xend=channeldensity, y=0, yend=cohend), color="grey") +
  geom_point( color="red", size=2) +
  scale_y_continuous(limits = c(0, 1.25)) +
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
  ylab("cohen d (power)") + 
  facet_grid(.~frequency_n, switch="both") 


data = all_canon2

g2 =  ggplot(data, aes(x=channeldensity_n, y=cohend, label=p_sign)) +
  geom_segment( aes(x=channeldensity_n, xend=channeldensity, y=0, yend=cohend), color="grey") +
  geom_point( color="darkgreen", size=2) +
  scale_y_continuous(limits = c(0, 1.25)) +
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
  ylab("cohen d (power)") + 
  facet_grid(.~frequency_n, switch="both") 


plot = ggarrange(g1, g2, 
                 labels = c("A", "B"),
                 ncol = 1, nrow = 2)

ggsave(file="global_d_channels_power_grey_fdr.png", plot=plot, dpi = 300, limitsize = TRUE, width = 12, height = 9)	


