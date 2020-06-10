
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
###############tables and figures related to House 116##################################
load(file=paste0('H116_pol.Rdata'),verbose=T)
load(file='rank_squad_dim1',verbose = T)

squad = c(195,216,210,270)
pol = pol_lower(pol)
Rep_gang = c(grep('\\Amash',pol),
             grep('\\Massie',pol),
             grep('\\Gaetz',pol))

col_stream = rep('red',nr)
col_stream[dem] = 'blue'
col_stream[ind] = 'green'
####load the posterior samples of ideal points in the euclidean model for House 116####
load(file='H116_beta_master_eu.Rdata',verbose=T)
rank_master_eu = apply(beta_master,2,rank)
rank_all_eu = round(t(apply(rank_master_eu,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all_eu) = pol
rm(beta_master)
####load the posterior samples of ideal points in the spherical model for House 116####
load(file='H116_beta_master_sph.Rdata',verbose=T)
rank_master = apply(beta_master,2,rank)
rank_all = round(t(apply(rank_master,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_all) = pol
rm(beta_master)
########################


#######Table 1###################################
table1 = cbind(to_table(rank_all_eu[squad,]),to_table(rank_squad_dim1))
table1
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

#################Figures related to House 112###############################################
load(file=paste0('H112_pol.Rdata'),verbose=T)
col_stream = rep('red',nr)
col_stream[dem] = 'blue'
col_stream[ind] = 'green'
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
top20 = which(rank(median_diff, ties.method='first')%in% c((nr-howmany+1):nr))
print(rev(tail(sort(median_diff),howmany)))
median_diff[top20]
pol[top20]

sph_top20 = rank_all_h112[top20,]
eu_top20 = rank_all_h112_eu[top20,]

l_to_s = order(-median_diff[top20])
sph_top20 = sph_top20[l_to_s,]
eu_top20 = eu_top20[l_to_s,]
###Remove Democrat legislator Kucinich since we are interested only in Republicans####
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

#######Figure 5###################################
x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
plot(rank_all_h112_eu[,2],rank_all_h112[,2],xlim=c(-5,nr+20),ylim=c(-5,nr+20),col=col_stream,
     xlab='Euclidean Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(rank_all_h112_eu[,2],rank_all_h112[,2],label = pol, cex = 1.5)





