### This script does z-transformation of correlation coefficients for comparison
### of group-results across resolutions (vertex, desikan, yeo)
### Christina Stier, 2025

## R version 4.2.2 (2022-10-31)
## RStudio 2023.3.0.386 for macOS

rm(list = ls())

setwd("~/Projects/Channels/R/files")

# get correlation coefficients (reference EEG metric at 256 channels, 3 resolutions)
# imcoh first
# r_list = as.data.frame(read.csv("~/Documents/Projects/Channels/R/files/corr_comparison_imcoh.csv", sep = ';'))
 r_list = as.data.frame(read.csv("~/Documents/Projects/Channels/R/files/corr_comparison_power.csv", sep = ';'))
str(r_list)

# follow code from https://studysites.uk.sagepub.com/dsur/study/answers.htm Chapter 6 on correlations (comparing independent rs)
# change code to my needs
zdiff  <-function(r1, r2, n1, n2)
{zd<-(atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3))
p <-1 - pnorm(abs(zd)) # gives one-tailed p
p_tt = p*2 # compute two-tailed p with p*2 if needed
stats = c(zd, p_tt)
print(stats)
}

# compute z(r) for each coefficient first (Fisher's z transformation => is equal to r's inverse hyperbolic tangent (artanh))
#r_list$z_vertex = atanh(r_list$r_vertex)
#r_list$z_desikan = atanh(r_list$r_desikan)
#r_list$z_yeo = atanh(r_list$r_yeo)

# now compute pairwise comparison: z_vertex vs. z_desikan; z_vertex vs. z_yeo
n_vertex = 2004
n_des = 68
n_yeo = 14

res_desikan = list()
res_yeo = list()

for ( c in 1:14){
  r_v = r_list$r_vertex[c]
  r_d = r_list$r_desikan[c]
  r_y = r_list$r_yeo[c]
  
  out_des = zdiff(r_v, r_d, n_vertex, n_des)
  out_yeo = zdiff(r_v, r_y, n_vertex, n_yeo)
  res_desikan[[c]] = out_des
  res_yeo[[c]] = out_yeo
}

full_res_desikan = do.call(rbind, res_desikan)
full_res_desikan = as.data.frame(full_res_desikan)
names(full_res_desikan)[1] = 'zdiff_des'
names(full_res_desikan)[2] = 'p_des'

full_res_yeo = do.call(rbind, res_yeo)
full_res_yeo = as.data.frame(full_res_yeo)
names(full_res_yeo)[1] = 'zdiff_yeo'
names(full_res_yeo)[2] = 'p_yeo'

r_list = cbind(r_list, full_res_desikan, full_res_yeo)

r_list$p_des_fdr[1:7] = p.adjust(r_list$p_des[1:7], method = "fdr", n = length(r_list$p_des[1:7]))
r_list$p_des_fdr[7:14] = p.adjust(r_list$p_des[7:14], method = "fdr", n = length(r_list$p_des[7:14])) # replaced length(r_list$p_des) as we wanted to correct only within hdms
r_list$p_yeo_fdr[1:7] = p.adjust(r_list$p_yeo[1:7], method = "fdr", n = length(r_list$p_yeo[1:7])) # replaced length(r_list$p_yeo) as we wanted to correct only within hdms
r_list$p_yeo_fdr[7:14] = p.adjust(r_list$p_yeo[7:14], method = "fdr", n = length(r_list$p_yeo[7:14])) 

# write.table(r_list, "comparison_r_imcoh_onetailed_fdr7.csv", sep = ";", row.names = FALSE) # save as csv
# write.table(r_list, "comparison_r_imcoh_twotailed_fdr7.csv", sep = ";", row.names = FALSE) # save as csv
 write.table(r_list, "comparison_r_power_twotailed_fdr7.csv", sep = ";", row.names = FALSE) # save as csv
# write.table(r_list, "comparison_r_power_onetailed_fdr7.csv", sep = ";", row.names = FALSE) # save as csv

