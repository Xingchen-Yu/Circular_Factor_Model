library(circular)

##############################
z_to_x = function(z){
  x = (z+pi2)/pi22
  return(x)
}
pi2 = pi^2
pi22 = 2*pi2
log_pi22 = log(pi22)
log_pi = log(pi)
##############################
plot_theta = function(n,a,b,ccc_a, ccc_b){
  kappa_a = 1
  omega = rgamma(n,a,b)

  beta = sapply(1:n,function(x) rvonmises(1,pi,omega[x])) -pi
  tau_yes = runif(n,-pi,pi)
  tau_no = runif(n,-pi,pi)

  asset = acos(cos(tau_no-beta))^2-acos(cos(tau_yes-beta))^2

  x = z_to_x(asset)
  kappa = rgamma(n,kappa_a,rgamma(n,ccc_a,ccc_b))
  sup = pbeta(x,kappa,kappa)
  return(sup)
}


set.seed(2021)
n_draws = 10000
theta_ij_default = plot_theta(n=n_draws,a=1,b=1/10,ccc_a = 1,ccc_b=25)
theta_ij_alt_1 = plot_theta(n=n_draws,a=7,b=1,ccc_a = 1,ccc_b=100)
theta_ij_alt_2 = plot_theta(n=n_draws,a=1,b=1/2,ccc_a = 1,ccc_b=25)

# x11(width=30, height=10)
# par(mar=c(6,6,3,3),mgp=c(4,1,0),mfrow=c(1,3))
# hist(theta_ij_default,xlab=bquote(~theta[ij]),cex.lab=3,col='white',ylab='Density',main='Default Prior',cex.main=3)
# hist(theta_ij_alt_1,xlab=bquote(~theta[ij]),cex.lab=3,col='white',ylab='Density',main='Alternative Prior 1 (favors Euclidean)',cex.main=3)
# hist(theta_ij_alt_2,xlab=bquote(~theta[ij]),cex.lab=3,col='white',ylab='Density',main='Alternative Prior 2 (favors Circular)',cex.main=3)


# theta_ij_alt_3 = plot_theta(n=n_draws,a=1,b=1/10,ccc_a = 1,ccc_b=5)
theta_ij_alt_4 = plot_theta(n=n_draws,a=1,b=1/10,ccc_a = 1,ccc_b=0.5)

# x11(width=20, height=10)
# par(mar=c(6,6,3,3),mgp=c(4,1,0),mfrow=c(1,2))
# hist(theta_ij_alt_3,xlab=bquote(~theta[ij]),cex.lab=2,col='white',ylab='Density',main='Alternative Prior 3',cex.main=2)
# hist(theta_ij_alt_4,xlab=bquote(~theta[ij]),cex.lab=2,col='white',ylab='Density',main='Alternative Prior 4',cex.main=2)

setwd("F:/Study_Backedup/UCSC/depository/Circular_Factor_Model/new_plot")
pdf(file = paste0("default_prior.pdf"),width=10, height=10)
par(mar=c(10,10,1,1),mgp=c(6,2,0))
hist(theta_ij_default,xlab=bquote(~upsilon[ij]),main='',cex.axis =2,cex.lab=2,probability = T,col='white')
dev.off()

pdf(file = paste0("Alt_prior_1.pdf"),width=10, height=10)
par(mar=c(10,10,1,1),mgp=c(6,2,0))
hist(theta_ij_alt_1,xlab=bquote(~upsilon[ij]),main='',cex.axis =2,cex.lab=2,probability = T,col='white')
dev.off()

pdf(file = paste0("Alt_prior_2.pdf"),width=10, height=10)
par(mar=c(10,10,1,1),mgp=c(6,2,0))
hist(theta_ij_alt_2,xlab=bquote(~upsilon[ij]),main='',cex.axis =2,cex.lab=2,probability = T,col='white')
dev.off()

pdf(file = paste0("Alt_prior_3.pdf"),width=10, height=10)
par(mar=c(10,10,1,1),mgp=c(6,2,0))
hist(theta_ij_alt_4,xlab=bquote(~upsilon[ij]),main='',cex.axis =2,cex.lab=2,probability = T,col='white')
dev.off()
