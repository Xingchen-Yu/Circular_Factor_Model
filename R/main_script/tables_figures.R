setwd("D:/Study_Backedup/UCSC/depository/Circular_Factor_Model")
library(stringr)
library(plotrix)

####helper function###
#####convert character to lower case except the first letter####
pol_lower = function(x){
  p_loc = str_locate(x, "\\(")[, 1]
  to_lower = sapply(1:length(x),function(i) substring(x[i],2,p_loc[i]-1))
  sapply(1:length(x), function(i) gsub(to_lower[i],tolower(to_lower[i]),x[i]))
}
######convert to paper format table
to_table = function(x){
  temp = cbind(paste0(x[,2]," (",x[,1],",",x[,3],")"))
  rownames(temp) = rownames(x)
  temp
}
###############tables and figures related to House 116##################################
load(file=paste0('H116_pol.Rdata'),verbose=T)
load(file='rank_squad_dim1',verbose = T)

squad = c(195,216,210,270)
pol = pol_lower(pol)
Rep_gang = c(grep('\\Amash',pol),
             grep('\\Massie',pol),
             grep('\\Gaetz',pol))


col_stream = rep('red',length(pol))
col_stream[grep('\\(D',pol)] = 'blue'
col_stream[grep('\\(I',pol)] = 'green'
#### load the posterior samples of ideal points in the euclidean model for House 116
load(file='H116_beta_master_eu.Rdata',verbose=T)
rank_master_eu = apply(beta_master,2,rank)
rank_all_eu = round(t(apply(rank_master_eu,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all_eu) = pol
rm(beta_master)
####load the posterior samples of ideal points in the spherical model for House 116
load(file='H116_beta_master_sph.Rdata',verbose=T)
rank_master = apply(beta_master,2,rank)
rank_all = round(t(apply(rank_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all) = pol
rm(beta_master)
########################
rownames(beta_master) = pol
library(tidyverse)
beta_master = as_tibble(t(beta_master))
write_csv(beta_master,file='./H116_posterior_beta_estimates.csv')
beta_master = read_csv(file='./H116_posterior_beta_estimates.csv')


#######Table 1###################################
table_1 = cbind(to_table(rank_all_eu[squad,]),to_table(rank_squad_dim1))
colnames(table_1) = c("Euclidean (1D)","Euclidean (2D)")
table_1
#######Table 2###################################
table_2 = to_table(rank_all[squad,])
colnames(table_2) = 'Rank Order (Circular)'
table_2 
#######Table 3###################################
table_3 = cbind(to_table(rank_all_eu[Rep_gang,]),
                to_table(rank_all[Rep_gang,]))
colnames(table_3) = c("Euclidean (1D)","Circular")
table_3 
#######Figure 3###################################
x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
plot(rank_all_eu[,2],rank_all[,2],xlim=c(-5,nrow(rank_all)),ylim=c(-5,nrow(rank_all)),col=col_stream,
     xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(rank_all_eu[,2],rank_all[,2],label = pol, cex = 1.5)

rm(list=setdiff(ls(), c("pol_lower","to_table")))
#####colorless Figure 3#####################################
x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
col_stream = rep('gray28',length(pol))
col_stream[grep('\\(D',pol)] = 'gray68'
col_stream[grep('\\(I',pol)] = 'gray68'

pch_stream = rep(17,length(rank_all_eu))
pch_stream[grep('\\(D',pol)] = 16
pch_stream[grep('\\(I',pol)] = 16 ## no independent in h116

plot(rank_all_eu[,2],rank_all[,2],xlim=c(-5,nrow(rank_all)),ylim=c(-5,nrow(rank_all)),col=col_stream,pch=pch_stream,
     xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',mgp=c(6,2,0),cex=1.5)

abline(1:500,1:500,col='gray',lwd=3)
legend(1,400, legend=c("Republican Party", "Democratic Party"),
       col=c("gray28", "gray68"), pch=c(17,16), cex=2,bty = "n")

identify(rank_all_eu[,2],rank_all[,2],label = pol, cex = 1.5)

rm(list=setdiff(ls(), c("pol_lower","to_table")))
#################Figures related to House 112###############################################
load(file=paste0('H112_pol.Rdata'),verbose=T)
pol = pol_lower(pol)
col_stream = rep('red',length(pol))
col_stream[grep('\\(D',pol)] = 'blue'
col_stream[grep('\\(I',pol)] = 'green'
####load the posterior samples of ideal points in the euclidean model for 112 House####
load(file='H112_beta_master_eu.Rdata',verbose=T)
rank_master_h112_eu = apply(beta_master,2,rank)
rank_all_h112_eu = round(t(apply(rank_master_h112_eu,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all_h112_eu) = pol
rm(beta_master)
####load the posterior samples of ideal points in the spherical model for 112 House####
load(file='H112_beta_master_sph.Rdata',verbose=T)
rank_master_h112 = apply(beta_master,2,rank)
rank_all_h112 = round(t(apply(rank_master_h112,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all_h112) = pol
rm(beta_master)

########Figure 4############################
howmany = 16

median_diff = abs(rank_all_h112[,2] - rank_all_h112_eu[,2])
top20 = which(rank(median_diff, ties.method='first')%in% c((length(pol)-howmany+1):length(pol)))
print(rev(tail(sort(median_diff),howmany)))
# median_diff[top20]
# pol[top20]
sph_top20 = rank_all_h112[top20,]
eu_top20 = rank_all_h112_eu[top20,]

l_to_s = order(-median_diff[top20])

sph_top20 = sph_top20[l_to_s,]
eu_top20 = eu_top20[l_to_s,]
###Remove Democrat legislator Kucinich since we are only interested in Republicans####
sph_top20 = sph_top20[-12,]
eu_top20 = eu_top20[-12,]

x11(height=15,width=15)
par(mar = c(14, 10, 2,2),mgp=c(6,2,0))
plot(seq(1,15), seq(1,15), type="n", ylim=c(200, 450), axes=FALSE, xlab="", ylab="Rank Order",cex.lab=2)
axis(1, at=seq(1,15), labels=rep("",15), las=2,cex.axis=1.5)
axis(2,las=2,cex.axis=1.5)
box()
plotCI(x=seq(1,15), y=eu_top20[,2], li=eu_top20[,1], ui=eu_top20[,3], add=TRUE, slty=2, pch=17, col=grey(0.6), scol=grey(0.6),cex=2)
plotCI(x=seq(1,15), y=sph_top20[,2], li=sph_top20[,1], ui=sph_top20[,3], add=TRUE, pch=16,cex=2)

legend(1.5, 385, c("Euclidean 1D","Circlar"), pch=c(17,16), lty=c(2,1), bty="n", col=c(grey(0.6), "black"),cex=1.5)

whobold <- c(1,2,3,5,6,7,8,11,12,13,14,15)
fontxaxis <- rep(1,15)
fontxaxis[whobold] <- 2
mtext(side=1, line=1, at=seq(1,15), text=rownames(sph_top20), font=fontxaxis,las=2,cex=1.5)

#######Figure 6###################################
x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
plot(rank_all_h112_eu[,2],rank_all_h112[,2],xlim=c(-5,length(pol)+20),ylim=c(-5,length(pol)+20),col=col_stream,
     xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(rank_all_h112_eu[,2],rank_all_h112[,2],label = pol, cex = 1.5)
#######colourless Figure 6#####################################

x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
col_stream = rep('gray28',length(pol))
col_stream[grep('\\(D',pol)] = 'gray68'
col_stream[grep('\\(I',pol)] = 'gray68'

pch_stream = rep(17,length(rank_all_h112_eu))
pch_stream[grep('\\(D',pol)] = 16
pch_stream[grep('\\(I',pol)] = 16 ## no independent in h116

plot(rank_all_h112_eu[,2],rank_all_h112[,2],xlim=c(-5,nrow(rank_all_h112)),ylim=c(-5,nrow(rank_all_h112)),col=col_stream,pch=pch_stream,
     xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',mgp=c(6,2,0),cex=1.5)

abline(1:500,1:500,col='gray',lwd=3)
legend(1,400, legend=c("Republican Party", "Democratic Party"),
       col=c("gray28", "gray68"), pch=c(17,16), cex=2,bty = "n")

identify(rank_all_h112_eu[,2],rank_all_h112[,2],label = pol, cex = 1.5)

rm(list=setdiff(ls(), c("pol_lower","to_table")))
########Tau_yes and Tau_no#########################################
library(circular)
library(plotrix)
hn = 112
if(hn==112){
  load("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/H112_workspace_grouped_waic_dic.Rdata")
}else{
  load("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/H116_workspace_grouped_waic_dic.Rdata")
}
yes_mean = apply(yes_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
no_mean = apply(no_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
# 
# 
# x11(width=10, height=10)
# bryant_yes = as.circular(yes_mean,units="radians",modulo="2pi",rotation="clock")
# bryant_no = as.circular(no_mean,units="radians",modulo="2pi",rotation="clock")
# plot.circular(bryant_yes[2,1],shrink=0.9,axes=F,units="radians",zero=pi/2,col='red',cex=1.5,pch=15)
# par(new=T)
# plot.circular(bryant_no[2,1],shrink=0.9,axes=F,units="radians",zero=pi/2,cex=1.5,col='blue',pch=17)
# 


plot_range = seq(1,nc,100)
nnn = length(plot_range)
if(hn==112){
plot_range[nnn] = nc
}else{
  nnn = nnn + 1
}

for(i in 1:(nnn-1)){
  aaa = plot_range[i]

  if(i==(nnn-1)){
    bbb = nc
  }else{
    bbb = plot_range[i+1]-1
  }

  plot_i = aaa:bbb
  
  pdf(file=paste0("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/h",hn,"/bill_",aaa,'_to_',bbb,'.pdf'),width=30, height=10)
  plot(plot_i, plot_i, type="n", ylim=c(-pi, pi), axes=T,xlab="Bill number", ylab="",cex.lab=2,cex.axis=1.5)
  
  plotCI(x=plot_i, y=yes_mean[2,plot_i], li=yes_mean[1,plot_i], ui=yes_mean[3,plot_i],slty=1, add=T,pch=17,cex.lab=2,cex.axis=1.5,lwd=2,bty='n',
         col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), scol=rgb(red = 0, green = 0, blue = 1, alpha = 0.5),cex=1.5,sfrac=0.001)
  plotCI(x=plot_i, y=no_mean[2,plot_i], li=no_mean[1,plot_i], ui=no_mean[3,plot_i], add=TRUE, pch=16,lwd=2,slty=2,
         col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5), scol=rgb(red = 1, green = 0, blue = 0, alpha = 0.5),cex=1.5,sfrac=0.001)
  dev.off()
}

#########################################
setwd("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model")
library(wnominate)
library(circular)
library(plotrix)
library(stringr)
pol_lower = function(x){
  p_loc = str_locate(x, "\\(")[, 1]
  to_lower = sapply(1:length(x),function(i) substring(x[i],2,p_loc[i]-1))
  sapply(1:length(x), function(i) gsub(to_lower[i],tolower(to_lower[i]),x[i]))
}

hn = 112
if(hn==112){
  load("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/H112_workspace_grouped_waic_dic.Rdata")
}else{
  load("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/H116_workspace_grouped_waic_dic.Rdata")
}
# yes_mean = apply(yes_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
# no_mean = apply(no_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))

yes_mean = apply(yes_master,1,mean.circular)
no_mean = apply(no_master,1, mean.circular)
# pol[135] = "GARCIA (D IL-4)"
pol = pol_lower(pol)

out = ymat_spit(hn=hn)
if(hn==116){
  print('Only session 1 (first 700 votes) of the 116 House data will be analyzed')
}
ymat = out[[1]]

# l_range = 0.5
# u_range = 2.5
# yes_0_no_pi = T
# if(yes_0_no_pi==T){
#   temp_mean1 = yes_mean
#   temp_mean2 = no_mean
# }else{
#   temp_mean1 = no_mean
#   temp_mean2 = yes_mean
# }

# yes_check = sapply(temp_mean1[2,],function(x) x<l_range & x>-l_range)
# no_check = sapply(temp_mean2[2,],function(x) x>u_range | x<= -u_range)
# check = which(yes_check+no_check == 2)
# out = apply(ymat[,check],2,function(x) table(x))
# if(is.list(out)==F){
#   check_wo_unani = check[1:ncol(out)]
# }else{
#   check_wo_unani = check[which(unlist(lapply(out,length))==2)]
# }
# out_new = apply(ymat[,check_wo_unani],2,function(x) table(x))
# nnn = length(check_wo_unani)
# name_list = vector('list',nnn)
# 
# for(i in 1:nnn){
#   Y_temp = ymat[,check_wo_unani[i]]
#   check_value = c(0,1)[as.numeric(which(table(Y_temp)==min(table(Y_temp))))]
#   name_list[[i]] = pol[which(Y_temp==check_value)]
# }
# names(name_list) = paste0('Bill ', check_wo_unani)
# 
# name_list
# bill_number = 531
# c(0,1)[which.min(table_list[[bill_number]])]
# names(which(ymat[,bill_number]==c(0,1)[which.min(table_list[[bill_number]])]))

table_list = apply(ymat,2,function(x) table(x,exclude = NA))

close = which(unlist(lapply(table_list,function(x) any(x/sum(x)<0.05) + length(x))==3)==T)
length(close)
plot_index = close[which(acos(cos(yes_mean[close]-no_mean[close]))>2.1)] ### 2.1
# plot_index = close[which(acos(cos(yes_mean[2,close]-no_mean[2,close]))>2.4)]
nnn = length(plot_index)
name_list = vector('list',nnn)



# bryant_yes = as.circular(yes_mean[2,],units="radians",modulo="2pi",rotation="clock")
# bryant_no = as.circular(no_mean[2,],units="radians",modulo="2pi",rotation="clock")
# for(i in 1:length(close)){
#   if(i %in% seq(1,length(close),16)){
#     x11(width = 15,height=15)
#     par(mfrow=c(4,4))
#   }
#   plot_k = close[i]
#   y_temp = table(ymat[,plot_k])
#   plot.circular(bryant_yes[plot_k],shrink=0.9,axes=F,units="radians",zero=pi/2,col='blue',cex=1.5,pch=15,main=paste0('No:',y_temp[1],", Yes: ",y_temp[2]))
#   par(new=T)
#   plot.circular(bryant_no[plot_k],shrink=0.9,axes=F,units="radians",zero=pi/2,cex=1.5,col='red',pch=17)
# 
# }

# names(which(Y_temp==c(0,1)[which.min(table_list[[bill_number]])]))
for(i in 1:nnn){
  Y_temp = ymat[,plot_index[i]]
  # check_value = c(0,1)[as.numeric(which(table(Y_temp)==min(table(Y_temp))))]
  
  name_list[[i]] = pol_lower(names(which(Y_temp==c(0,1)[which.min(table_list[[plot_index[i]]])])))
  
}
names(name_list) = paste0('Bill ', plot_index)
name_list

# dev.new(width=15, height=15,noRStudioGD = T)
# # par(mar=c(7,7,2,2))
# beta_mean = apply(beta_master,1,median)
# names(beta_mean) = pol


# ymat[squad,48]
# acos(cos(yes_mean[2,48]-no_mean[2,48]))


# bryant_yes = as.circular(yes_mean[2,],units="radians",modulo="2pi",rotation="clock")
# bryant_no = as.circular(no_mean[2,],units="radians",modulo="2pi",rotation="clock")
for(i in 1:nnn){
  plot_k = plot_index[i]
  pdf(file=paste0("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/h",hn,"/Bill_",plot_k,'.pdf'),width = 15,height=15)
  bryant = circular::as.circular(as.numeric(beta_mean),units="radians",modulo="2pi",rotation="clock")
  circular::plot.circular(bryant[grep("\\(D",pol)],shrink=1,axes=F,units="radians",zero=pi/2,col='gray38',cex=1.5)
  par(new=T)
  circular::plot.circular(bryant[-grep("\\(D",pol)],shrink=1,axes=F,units="radians",zero=pi/2,col='gray68',cex=1.5,pch=17)
  
  b_name = beta_mean[name_list[[i]]]
  for(jj in 1:length(b_name)){
    b_name_jj = b_name[jj]
    text(sin(b_name_jj),cos(b_name_jj),labels=names(b_name_jj),cex=1.5,offset = 1, font=1)
  }
  
  
  # par(new=T)
  # plot.circular(bryant_yes[plot_k],shrink=1,axes=F,units="radians",zero=pi/2,cex=4,cex.main=3)
  text(sin(yes_mean[2,plot_k]),cos(yes_mean[2,plot_k]),labels=bquote(symbol("\326")),cex=4,offset = 0, font=1)
  # par(new=T)
  # plot.circular(bryant_no[plot_k],shrink=1,axes=F,units="radians",zero=pi/2,cex=4,col='orange',pch=16)
  text(sin(no_mean[2,plot_k]),cos(no_mean[2,plot_k]),labels=bquote("X"),cex=4,offset = 0, font=1)
  
  
  dev.off()
}


kobe = read.csv('F:/Study_Backedup/UCSC/LatentFactorModel/PAP.csv',header=T)

kobe_voteview = read.csv(paste0('F:/Study_Backedup/UCSC/LatentFactorModel/H',hn,'_rollcalls.csv'),header=T)
# head(kobe_voteview[which(kobe_voteview$bill_number=='HR1540'),])
# str(kobe)
# str(kobe_voteview)

kobe_hn = kobe[intersect(which(kobe$filter_House==1),which(kobe$cong==hn)),]
# str(kobe_hn)

# kobe_voteview$vote_desc[plot_index]
# kobe_voteview[plot_index,c("vote_result","vote_desc")]

# kobe_voteview[plot_index,c("bill_number")] 


# which(kobe_hn$bill %in% kobe_voteview[plot_index[-5],c("bill_number")] )

# kobe_hn_subset = kobe_hn[which(kobe_hn$bill %in% kobe_voteview[plot_index[-5],c("bill_number")] ),]

# condition_1 = which(kobe_hn$sess_count %in% kobe_voteview[plot_index,c('clerk_rollnumber')])
# kobe_hn[condition_1,"sess"] == kobe_voteview[plot_index,c('session')]
if(hn==112){

index1 = which(kobe_voteview[plot_index,c('session')]==1)
index2 = which(kobe_voteview[plot_index,c('session')]==2)

kobe_hn_subset1 = kobe_hn[which(kobe_hn$sess==1),][which(kobe_hn[which(kobe_hn$sess==1),]$sess_count %in% kobe_voteview[plot_index[index1],c('clerk_rollnumber')]),]
kobe_hn_subset2 = kobe_hn[which(kobe_hn$sess==2),][which(kobe_hn[which(kobe_hn$sess==2),]$sess_count %in% kobe_voteview[plot_index[index2],c('clerk_rollnumber')]),]

kobe_hn_subset = rbind(kobe_hn_subset1,kobe_hn_subset2)

cbind(kobe_hn_subset$sess_count,kobe_voteview[plot_index,'clerk_rollnumber'])
cbind(kobe_hn_subset$bill,kobe_voteview[plot_index,'bill_number'])
head(kobe_hn_subset)
head(kobe_voteview)


output = cbind(kobe_voteview[plot_index,c("bill_number","vote_result","vote_desc")],kobe_hn_subset[,14:17])
}else{
  output = kobe_voteview[plot_index,c("bill_number","vote_result","vote_desc")]
}


write.csv(file=paste0("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/h",hn,"/bill_detail.csv"),output)

################################################
###101 and 111 house
hn = 111
load(paste0("F:/Study_Backedup/UCSC/LatentFactorModel/House/DIC_WAIC_GROUP/sph/H",hn,"_workspace_grouped_waic_dic.Rdata"))
rank_master = apply(beta_master,2,rank)
rank_all = round(t(apply(rank_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all) = pol

col_stream = rep('red',length(pol))
col_stream[grep('\\(D',pol)] = 'blue'
col_stream[grep('\\(I',pol)] = 'green'
rm(list=setdiff(ls(), c("col_stream","pol","rank_all",'hn')))


load(paste0("F:/Study_Backedup/UCSC/LatentFactorModel/House/DIC_WAIC_GROUP/eu/H",hn,"_1_nd_eu_workspace.Rdata"))
rank_master_eu = apply(beta_master,2,rank)
rank_all_eu = round(t(apply(rank_master_eu,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all_eu) = pol

rm(list=setdiff(ls(), c("col_stream","pol","rank_all",'rank_all_eu')))
x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
plot(rank_all_eu[,2],rank_all[,2],xlim=c(-5,nrow(rank_all)),ylim=c(-5,nrow(rank_all)),col=col_stream,
     xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(rank_all_eu[,2],rank_all[,2],label = pol, cex = 1.5)
################################################
hn = 116
library(wnominate)
if(hn==112){
  load("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/H112_workspace_grouped_waic_dic.Rdata")
}else{
  load("F:/Study_Backedup/UCSC/LatentFactorModel/House/h116_112_tau_yes_tau_no_additional_runs/H116_workspace_grouped_waic_dic.Rdata")
}
setwd("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model")
out = ymat_spit(hn)
ymat = out[[1]]
pol = pol_lower(pol)
sum_bill = apply(ymat,2,function(x) sum(x,na.rm = T))
dem = grep('\\(D',pol)
gop = grep('\\(R',pol)
l_dem = length(dem)
l_gop = length(gop)
# ind = grep('\\(I',pol)
sum_bill_dem = apply(ymat,2,function(x) sum(x[dem],na.rm = T))
sum_bill_gop = apply(ymat,2,function(x) sum(x[gop],na.rm = T))
sum_by_party = cbind(sum_bill_dem,sum_bill_gop)
colnames(sum_by_party) = c('Dem','GOP')
diff_party = abs(apply(sum_by_party,1,diff))
bill_split_party = which(diff_party == max(diff_party))
sum_by_party[bill_split_party,]
party_n = length(bill_split_party)
for(i in 1:party_n){
  if(i %in% seq(1,party_n,4)){
    x11(height=15,width=15)
    par(mfrow=c(2,2))
    par(mar = c(10, 10, 4, 4))
  }
  
  hist(kappa_master[bill_split_party[i],],cex.axis=2,cex.lab=2,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
  abline(v=25,lwd=2)
  
}

sum_bill[order(sum_bill)]
una_bill = which(sum_bill==424)
table(ymat[,32])
# table(ymat[,81])
# table(ymat[,149])
which(is.na(ymat[,32])==T)
which(ymat[,una_bill]==0)



una_n= length(una_bill)
# for(i in 1:una_n){
#   if(i %in% seq(1,una_n,9)){
#     x11(height=15,width=15)
#     par(mfrow=c(3,3))
#     par(mar = c(10, 10, 4, 4))
#   }
#     
#   hist(kappa_master[una_bill[i],],cex.axis=2,cex.lab=2,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
#   abline(v=25,lwd=2)
#   
# }

circ_bill = c(496,531)
range_kappa = range(kappa_master[una_bill[1],],kappa_master[circ_bill[1],],kappa_master[circ_bill[2],],bill_split_party[1])
x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
hist(kappa_master[una_bill[1],],xlim = c(0,range_kappa[2]),cex.axis=2,cex.lab=3,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
abline(v=25,lwd=3)

x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
hist(kappa_master[bill_split_party[1],],xlim = c(0,range_kappa[2]),cex.axis=2,cex.lab=3,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
abline(v=25,lwd=3)

x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
hist(kappa_master[circ_bill[1],],xlim = c(0,range_kappa[2]),cex.axis=2,cex.lab=3,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
abline(v=25,lwd=3)

x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
hist(kappa_master[circ_bill[2],],xlim = c(0,range_kappa[2]),cex.axis=2,cex.lab=3,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
abline(v=25,lwd=3)



circ_bill = c(1553)

x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
hist(kappa_master[circ_bill[1],],cex.axis=2,cex.lab=3,cex.main=2,main='',xlab = bquote(~kappa[j]),mgp=c(6,2,0),ylab='')
abline(v=25,lwd=3)


library(circular)
# beta_mean = apply(beta_master,1,median)
beta_mean = apply(beta_master,1,circular::mean.circular)
names(beta_mean) = pol

# yes_mean = apply(yes_master,1,function(x) quantile.circular(x,probs = c(0.025,0.5,0.975)))
# no_mean = apply(no_master,1,function(x) quantile.circular(x,probs = c(0.025,0.5,0.975)))

yes_mean2 = apply(yes_master,1,mean.circular)
no_mean2 = apply(no_master,1,mean.circular)

# bryant_yes = as.circular(yes_mean[2,],units="radians",modulo="2pi",rotation="clock")
# bryant_no = as.circular(no_mean[2,],units="radians",modulo="2pi",rotation="clock")
plot_k = 1553 ## unanimous bills for h116 in paper are 5(split by party) and 32(unanimous)
  x11(width = 15,height=15)
  bryant = circular::as.circular(as.numeric(beta_mean),units="radians",modulo="2pi",rotation="clock")
  circular::plot.circular(bryant[grep("\\(D",pol)],shrink=1,axes=F,units="radians",zero=pi/2,col='gray28',cex=1.5)
  par(new=T)
  circular::plot.circular(bryant[-grep("\\(D",pol)],shrink=1,axes=F,units="radians",zero=pi/2,col='gray68',cex=1.5,pch=17)
  if(plot_k==496){
    b_name = beta_mean[name_list$`Bill 496`]
  }else if(plot_k==531){
    b_name = beta_mean[name_list$`Bill 531`]
  }else if(plot_k==1553){
    b_name = beta_mean[name_list$`Bill 1553`]
  }
  if(plot_k %in% c(496,531,1553)){
    for(jj in 1:length(b_name)){
      b_name_jj = b_name[jj]
      text(sin(b_name_jj),cos(b_name_jj),labels=names(b_name_jj),cex=1.5,offset = 1, font=1)
    }
  }
  text(sin(yes_mean2[plot_k]),cos(yes_mean2[plot_k]),labels=bquote(symbol("\326")),cex=4,offset = 0, font=1)
  text(sin(no_mean2[plot_k]),cos(no_mean2[plot_k]),labels=bquote("X"),cex=4,offset = 0, font=1)
  
  # yes_median = median.circular(yes_master[plot_k,])
  # no_median = median.circular(no_master[plot_k,])
  # text(sin(yes_median),cos(yes_median),labels=bquote(symbol("\326")),cex=4,offset = 0, font=1)
  # text(sin( no_median ),cos( no_median ),labels=bquote("X"),cex=4,offset = 0, font=1)
# yes_median = median.circular(yes_master[plot_k,])
# 
  
  yes_mean_3 = apply(yes_master,1,function(x) quantile.circular(x,probs = c(0.025,0.5,0.975)))
  no_mean_3 = apply(no_master,1,function(x) quantile.circular(x,probs = c(0.025,0.5,0.975)))
# yes_mean[plot_k]
# median.circular(c(pi-0.1,pi-0.2,-pi+0.1,-pi+0.2))
# 
# median.circular(no_master[496,])
# mean.circular(no_master[496,])
# no_mean[2,496]
  quantile.circular(yes_master[plot_k,],probs = c(0.025,0.5,0.975))
  median.circular(yes_master[plot_k,])
  
  quantile.circular(no_master[plot_k,],probs = c(0.025,0.5,0.975))
  median.circular(no_master[plot_k,])
  
  no_mean[,plot_k]
  
  
  