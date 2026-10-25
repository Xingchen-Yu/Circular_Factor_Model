library(wnominate)
library(circular)
library(plotrix)
# setwd("F:\\Study_Backedup\\UCSC\\LatentFactorModel")
house_check = F
if(house_check==T){
  setwd("F:\\Study_Backedup\\UCSC\\LatentFactorModel\\House\\House_beta\\")
  h_s = 'H'
}else{
  setwd("F:\\Study_Backedup\\UCSC\\LatentFactorModel\\Senate\\sph/Senate_beta/")
  h_s = 'S'
}
source(file="F:\\Study_Backedup\\UCSC\\code_master/code_master/read_kh2.R")
source(file="F:\\Study_Backedup\\UCSC\\code_master/code_master/ymat_spit.R")
# hn_master = c(104:107,109,111:114)
# for(i in 1:length(hn_master)){
#   hn = hn_master[i]
#   out = ymat_spit(hn=hn,T)
#   ymat = out[[1]]
#   nr = nrow(ymat)
#   nc = ncol(ymat)
#   pol = out[[2]]
#   beta = t(read.csv(file=paste0(".\\House_beta\\beta_master_",hn,'.csv'),header=F))
# 
#   beta_master = matrix(beta,nr,length(beta)/nr)
#   save(file=paste0('.\\House_beta\\H',hn,"_betaChain.Rdata"),beta_master)
# 
# }

CI_mean<-function(what){
  return(quantile(what,probs = c(0.025,0.5,0.975)))
}
hn_master = c(100:116)
ln = length(hn_master)
col_s = c(rep('blue',4),rep('red',6),rep('blue',2),rep('red',4),'blue')
diff_master = v_pos_vec_gop =v_pos_vec_dem =v_pos_vec = matrix(0,3,ln)
var_mean_chain = diff_mean_master = rep(0,ln)
dev.new(width=15, height=15,noRStudioGD = T)
par(mfrow=c(5,5))
for(i in ln:length(hn_master)){
  hn = hn_master[i]
  if(hn==116 & house_check==T){
    load('h116.Rdata',verbose=T)
    load('h116_pol.Rdata',verbose= T)
  }else{
    out = ymat_spit(hn=hn,house_check)
    ymat = out[[1]]
    pol = out[[2]]
  }
  nr = nrow(ymat)
  nc = ncol(ymat)
  rm(ymat)
  rm(out)
  
  # load(file=paste0('.\\House_beta\\H',hn,"_betaChain.Rdata"),verbose=T)
  load(file=paste0(h_s,hn,"_beta_master.Rdata"),verbose=T)
  beta_mean = apply(beta_master,1,mean.circular)


  dem = grep("\\(D",pol) 
  gop = grep("\\(R",pol) 
  # if(hn==112){
  # gop = which(rank(-beta_mean)%in% c(1:length(gop)))
  # dem = which(rank(beta_mean)%in% c(1:length(dem)))
  # }
  
  var_chain = apply(beta_master,2,var.circular)
  var_chain_gop = apply(beta_master[gop,],2,var.circular)
  var_chain_dem = apply(beta_master[dem,],2,var.circular)
  diff_chain = apply(beta_master,2,function(x) abs(mean.circular(x[dem])-mean.circular(x[gop])))
  
  diff_mean_master[i] = mean.circular(diff_chain)
  var_mean_chain[i] = mean(var_chain)
  
  v_pos_vec[,i] = CI_mean(var_chain)
  v_pos_vec_gop[,i] = CI_mean(var_chain_gop)
  v_pos_vec_dem[,i] = CI_mean(var_chain_dem)
  diff_master[,i] = CI_mean(diff_chain)
  
  plot.circular(beta_mean[gop],axes=F,shrink=0.9,units="radians",zero=pi/2,col=rgb(1,0,0,0.5),cex=1,main=paste0('House',hn))
  par(new=T)
  plot.circular(beta_mean[dem],axes=F,shrink=0.9,units="radians",zero=pi/2,col=rgb(0,0,1,0.5),cex=1)

  # rm(beta_master,beta_mean,var_chain,var_chain_gop,var_chain_dem,dem,gop,nr,nc,pol,ymat,out)
  rm(beta_master)
}
#####################################################################################################################

dev.new(width=30, height=15,noRStudioGD = T)
par(mar = c(8, 10, 2, 2))

ylim_1 = c(min(v_pos_vec) - 0.02,max(v_pos_vec) + 0.02)
plot(0,xaxt = 'n',ylim=ylim_1,xlim=c(1,ln),las=1,xlab='',xaxt = 'n',cex.axis=2,ylab="Circular Variance",cex.axis=2.5,cex.lab=3,mgp=c(7,1,0))

if(house_check==T){
  rect(0, min(v_pos_vec) - 1, 4+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(4+0.5, min(v_pos_vec) - 1, 10+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect( 10+0.5, min(v_pos_vec) - 1, 12+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(12+0.5, min(v_pos_vec) - 1, 16+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect(16+0.5, min(v_pos_vec) - 1, 17+1, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  plotCI(v_pos_vec[2,],ui=v_pos_vec[3,],li=v_pos_vec[1,],main="",sfrac=0.005,ylim=ylim_1,slty = 1,pch=18,cex=3,mgp=c(6,5,0),
         yaxt = 'n',ylab="",xaxt = 'n',xlab='U.S House of Representatives',lwd=4,cex.lab=3,col='black')
  axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2.5,mgp=c(0,2,0))
}else{
  # c("D","D",'D','D','R','R','R','NA','R','R','D','D','D','D','R','R','R')
  
  rect(0, min(v_pos_vec) - 1, 4+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(4+0.5, min(v_pos_vec) - 1, 7+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect(8+0.5, min(v_pos_vec) - 1, 10+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(10+0.5, min(v_pos_vec) - 1, 14+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect(14+0.5, min(v_pos_vec) - 1, 17+1, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  plotCI(v_pos_vec[2,],ui=v_pos_vec[3,],li=v_pos_vec[1,],main="",sfrac=0.005,ylim=ylim_1,slty = 1,pch=18,cex=3,mgp=c(6,5,0),
         yaxt = 'n',ylab="",xaxt = 'n',xlab='U.S Senate',lwd=4,cex.lab=3,col='black')
  axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2.5,mgp=c(0,2,0))
  
}

# legend(0.5,0.45, legend=c("Republican Party", "Democratic Party"),
#        col=c("red", "blue"), lty=rep(1,2), cex=2,lwd=3)




#####################################################################################################################
# dev.new(width=30, height=15,noRStudioGD = T)
# par(mar=c(7,7,2,2))
# 
# plotCI(v_pos_vec_dem[2,],ui=v_pos_vec_dem[3,],li=v_pos_vec_dem[1,],main="",sfrac=0.005,
#        xaxt = 'n',ylab="",xlab='',lwd=2,cex.lab=2,col='blue',ylim=c(0,0.05))
# axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2,cex.main=3)
# #####################################################################################################################
# dev.new(width=30, height=15,noRStudioGD = T)
# par(mar=c(7,7,2,2))
# par(new=T)
# plotCI(v_pos_vec_gop[2,],ui=v_pos_vec_gop[3,],li=v_pos_vec_gop[1,],main="",sfrac=0.005,
#        xaxt = 'n',ylab="Circular Variance",xlab='U.S House of Representatives',lwd=2,cex.lab=2,col='red',ylim=c(0,0.05))
# axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2,cex.main=3)
# 
# legend(0.5,0.05, legend=c("Republican Party", "Democratic Party"),
#        col=c("red", "blue"), lty=rep(1,2), cex=2,lwd=3)

dev.new(width=30, height=15,noRStudioGD = T)
par(mar = c(8, 10, 2, 2))

ylim_1 = c(0,max(c(v_pos_vec_gop,v_pos_vec_dem)) + 0.02)
plot(0,xaxt = 'n',ylim=ylim_1,xlim=c(1,ln),las=1,xlab='',xaxt = 'n',cex.axis=2,ylab="Circular Variance",cex.axis=2.5,cex.lab=3,mgp=c(7,1,0))

if(house_check==T){
  rect(0, min(v_pos_vec) - 1, 4+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(4+0.5, min(v_pos_vec) - 1, 10+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect( 10+0.5, min(v_pos_vec) - 1, 12+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(12+0.5, min(v_pos_vec) - 1, 16+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect(16+0.5, min(v_pos_vec) - 1, 17+1, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  plotCI(v_pos_vec_gop[2,],ui=v_pos_vec_gop[3,],li=v_pos_vec_gop[1,],main="",sfrac=0.005,ylim=ylim_1,slty = 1,cex=2.5,mgp=c(6,5,0),
         yaxt = 'n',ylab="",xlab='U.S House of Representatives',xaxt = 'n',lwd=4,cex.lab=3,col='red',pch=15)
}else{
  rect(0, min(v_pos_vec) - 1, 4+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(4+0.5, min(v_pos_vec) - 1, 7+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect(8+0.5, min(v_pos_vec) - 1, 10+0.5, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  rect(10+0.5, min(v_pos_vec) - 1, 14+0.5, max(v_pos_vec) + 1, col = 'gray38', lty=0)
  par(new=T)
  rect(14+0.5, min(v_pos_vec) - 1, 17+1, max(v_pos_vec) + 1, col = 'gray88', lty=0)
  par(new=T)
  plotCI(v_pos_vec_gop[2,],ui=v_pos_vec_gop[3,],li=v_pos_vec_gop[1,],main="",sfrac=0.005,ylim=ylim_1,slty = 1,cex=2.5,mgp=c(6,5,0),
         yaxt = 'n',ylab="",xlab='U.S Senate',xaxt = 'n',lwd=4,cex.lab=3,col='red',pch=15)
}


plotCI(v_pos_vec_dem[2,],ui=v_pos_vec_dem[3,],li=v_pos_vec_dem[1,],main="",sfrac=0.005,ylim=ylim_1,slty = 1,cex=2.5,mgp=c(6,5,0),
       yaxt = 'n',ylab="",xaxt = 'n',lwd=4,cex.lab=3,col='blue',add=T,pch=18)
axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2.5,mgp=c(0,2,0))

legend(0.7, 0.07, c("Democratic Party","Republican Party"), pch=c(18,15), lwd=4,lty=c(1,1),cex=2.5, bty="n", col=c('blue', "red"))
#####################################################################################################################
# dev.new(width=30, height=15,noRStudioGD = T)
# par(mar=c(7,7,2,2))
# plotCI(diff_master[2,],ui=diff_master[3,],li=diff_master[1,],main="",sfrac=0.005,
#        xaxt = 'n',ylab="Absolute Difference of Party Mean",xlab='U.S House of Representatives',lwd=2,cex.lab=2,col=col_s)
# 
# axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2,cex.main=3)
# legend(0.5,1.7, legend=c("Republican Party", "Democratic Party"),
#        col=c("red", "blue"), lty=rep(1,2), cex=2,lwd=3)

dev.new(width=30, height=15,noRStudioGD = T)
par(mar = c(8, 10, 2, 2))

ylim_1 = c(min(diff_master) - 0.05,max(diff_master) + 0.05)
plot(0,xaxt = 'n',ylim=ylim_1,xlim=c(1,ln),las=1,xlab='',xaxt = 'n',cex.axis=2,ylab="Absolute Difference of Party Mean",cex.axis=2.5,cex.lab=3,mgp=c(7,1,0))

rect(0, min(diff_master) - 1, 4+0.5, max(diff_master) + 1, col = 'gray88', lty=0)
par(new=T)
rect(4+0.5, min(diff_master) - 1, 10+0.5, max(diff_master) + 1, col = 'gray38', lty=0)
par(new=T)
rect( 10+0.5, min(diff_master) - 1, 12+0.5, max(diff_master) + 1, col = 'gray88', lty=0)
par(new=T)
rect(12+0.5, min(diff_master) - 1, 16+0.5, max(diff_master) + 1, col = 'gray38', lty=0)
par(new=T)
rect(16+0.5, min(diff_master) - 1, 17+1, max(diff_master) + 1, col = 'gray88', lty=0)
par(new=T)
plotCI(diff_master[2,],ui=diff_master[3,],li=diff_master[1,],main="",sfrac=0.005,ylim=ylim_1,slty = 1,pch=18,cex=3,mgp=c(6,5,0),
       yaxt = 'n',ylab="",xaxt = 'n',xlab='U.S House of Representatives',lwd=4,cex.lab=3,col='black')
axis(1,at=1:ln,labels=paste0(rep(hn_master)), las = 1,cex.axis=2.5,mgp=c(0,2,0))



