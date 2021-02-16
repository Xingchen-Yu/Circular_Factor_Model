library(circular)
library(mvnfast)
source(file="./source/ymat_spit.R")
set.seed(2021)
out = ymat_sim(100,700,1,50,1,1/2)
save(file='circular_data.Rdata',out)

set.seed(2021)
ymat = ymat_sim_eu(100,700,2)
save(file='euclidean_2d_data.Rdata',ymat)

##spherical 1d data
waic_eu_vec = rep(0,3)
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_eu2d_data_1d_eu.RData")
waic_eu_vec[1] = waic_group[1] - waic_group[3]
plot(likeli_chain)
rm(list=setdiff(ls(),c('waic_eu_vec')))
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_eu2d_data_circ.RData")
waic_eu_vec [2] = waic_group[1] - waic_group[3]
plot(likeli_chain[burnin:iter])
rm(list=setdiff(ls(),c('waic_eu_vec')))
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_eu2d_data_2d_eu.RData")
waic_eu_vec [3] = waic_group[1] - waic_group[3]
plot(likeli_chain)
rm(list=setdiff(ls(),c('waic_eu_vec')))

waic_eu_vec

##spherical 1d data
waic_circ_vec = rep(0,3)
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_circ_data_1d_eu.RData")
waic_circ_vec[1] = waic_group[1] - waic_group[3]
plot(likeli_chain)
rm(list=setdiff(ls(),c('waic_circ_vec')))
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_circ_data_circ.RData")
waic_circ_vec [2] = waic_group[1] - waic_group[3]
plot(likeli_chain[burnin:iter])
rm(list=setdiff(ls(),c('waic_circ_vec')))
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_circ_data_2d_eu.RData")
waic_circ_vec [3] = waic_group[1] - waic_group[3]
plot(likeli_chain)
rm(list=setdiff(ls(),c('waic_circ_vec')))

waic_circ_vec
############################
###
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_circ_data_2d_eu.RData")
beta_mean = beta_rs/n_pos

dev.new(height=15,width=15,noRStudioGD = T)
par(mar = c(10, 10, 4, 4))
plot(beta_mean,xlab=bquote(~beta[i~','~1]),cex=2, ylab=bquote(~beta[i~','~2]),cex.axis=2,cex.lab=2,mgp=c(6,2,0),lwd=2,pch=16,col='gray28')

##
load("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/simulation_results/paper1_eu2d_data_2d_eu.RData")
beta_mean = beta_rs/n_pos

dev.new(height=15,width=15,noRStudioGD = T)
par(mar = c(10, 10, 4, 4))
plot(beta_mean,xlab=bquote(~beta[i~','~1]),cex=2, ylab=bquote(~beta[i~','~2]),cex.axis=2,cex.lab=2,mgp=c(6,2,0),lwd=2,pch=16,col='gray28')








