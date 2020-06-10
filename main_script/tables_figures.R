
library(stringr)
library(plotrix)

####helper function###
#####convert character to lower case except the first letter####
pol_lower = function(x){
  p_loc = str_locate(x, "\\(")[, 1]
  to_lower = sapply(1:length(x),function(i) substring(x[i],2,p_loc[i]-1))
  sapply(1:length(x), function(i) gsub(to_lower[i],tolower(to_lower[i]),x[i]))
}
######convert to paper format table#####
to_table = function(x){
  temp = cbind(paste0(x[,2]," (",x[,1],",",x[,3],")"))
  rownames(temp) = rownames(x)
  temp
}
###############
load(file=paste0('H116_pol.Rdata'),verbose=T)
load(file='rank_squad_dim1',verbose = T)

squad = c(195,216,210,270)
pol = pol_lower(pol)
Rep_gang = c(grep('\\Amash',pol),
             grep('\\Massie',pol),
             grep('\\Gaetz',pol))
####load the posterior samples of ideal points in the euclidean model####
load(file=paste0('H',hn,"_beta_master_eu.Rdata"),verbose=T)
rank_master_eu = apply(beta_master,2,rank)
rank_all_eu = round(t(apply(rank_master_eu,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all_eu) = pol
rm(beta_master)
####load the posterior samples of ideal points in the spherical model####
load(file=paste0('H',hn,"_beta_master_sph.Rdata"),verbose=T)
rank_master = apply(beta_master,2,rank)
rank_all = round(t(apply(rank_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all) = pol
rm(beta_master)
########################


#######Table 1######
table1 = cbind(to_table(rank_all_eu[squad,]),to_table(rank_squad_dim1))
table1
#######Table 2#######
table_2 = to_table(rank_all[squad,])
colnames(table_2) = 'Rank Order (Circular)'
table_2 
#######Table 3#######
table_3 = cbind(to_table(rank_all_eu[Rep_gang,]),
                to_table(rank_all[Rep_gang,]))
colnames(table_3) = c("Euclidean (1D)","Circular")
table_3 
#######Figure 3#######
dev.new(height=15,width=15,noRStudioGD = T)
par(mar = c(10, 10, 4, 4))
plot(rank_all_eu[,2],rank_all[,2],xlim=c(-5,nrow(rank_all)),ylim=c(-5,nrow(rank_all)),col=col_stream,xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(rank_all_eu[,2],rank_all[,2],label = pol, cex = 1.5)

rm(list=setdiff(ls(), c("pol_lower","to_table")))





