library(mvnfast)
library(truncnorm)
library(wnominate)
library(coda)
library("plotrix")
#library(snowfall)
source(file="Euclidean_Functions.R")
####################################
hn = 116
seed = 2020
alpha_penalty = F
####################################
test = F
House = T
h_s = ifelse(House==T,'H','S')
####################################  
p = 2
iter = 40000
burnin = 20000
n_pos = iter - burnin 
likeli_chain = rep(0,n_pos)
####################################
p_mat<-diag((1:p)^2)

if(hn!=116){
  out = ymat_spit(hn=hn,House)
  ymat = out[[1]]
  pol = out[[2]]
}else{
  load('h116.Rdata',verbose=T)
  load('h116_pol.Rdata',verbose= T)
  squad = c(195,216,210,270)
  print(pol[squad])
  beta_squad_dim_2 = beta_squad_dim_1 = matrix(0,length(squad),n_pos)
}
print(dim(ymat))
nr<-nrow(ymat)
nc<-ncol(ymat)
a_mu<-rep(0,p+1)
if(alpha_penalty==T){
  sig_mu_alpha = diag(c(1,1/c(1:p)^2))
}else{
  sig_mu_alpha<-diag(p+1)
}
sig_mu_alpha_inv<-solve(sig_mu_alpha)

beta_rs<-matrix(0,nr,p)
mu_rs<-rep(0,nc)
alpha_rs<-matrix(0,nc,p)
zmat<-matrix(0,nr,nc)

if(test == T){
  na_ori_true = which(is.na(ymat)==T)
  na.position<-which(is.na(ymat)==T, arr.ind = T)
  set.seed(seed)
  for(j in 1:nc){
    tom = na.position[which(na.position[,2]==j),1]
    if(length(tom)>0){
      ymat[sample(c(1:nr)[-tom],ceiling(0.05 * (nr-length(tom)))),j] = NA
    }else{
      ymat[sample(c(1:nr),ceiling(0.05 * nr)),j] = NA
    }
  }
  temp_na = which(is.na(ymat)==T)
  na_test = temp_na[which(temp_na %!in% na_ori_true)]
  len_test = length(na_test)
  y_test = ymat_true[na_test]
}
pos_pred<-pos_pred2<-pos_pred3<-matrix(0,nr,nc)

yes<-which(ymat==1)
yes_l<-length(yes)
no<-which(ymat==0)
no_l<-length(no)
na<-which(is.na(ymat)==T)

no_na<-which(is.na(ymat)==F)
len_no_na = length(no_na)
na_l<-length(na)

yes_mean<-rep(0,yes_l)
no_mean<-rep(0,no_l)
na_mean<-rep(0,na_l)

beta<-matrix(0,nr,p)
cov_master =rep(NA,n_pos)
part2<-sig_mu_alpha_inv%*%a_mu
##
mu_alpha_mat<-matrix(0,nc,p+1)
douche<-function(kobe){
  sup<-rmvn(1,c(kobe),var_beta)
  return(sup)
}

jj = 1
for(i in 1:iter){
  
  ###latent###
  z_yes<-rtruncnorm(1,a=0,b=Inf,mean=yes_mean,sd=1)
  z_no<-rtruncnorm(1,a=-Inf,b=0,no_mean,sd=1)
  z_na<-rnorm(na_l,na_mean,sd=1)
  zmat[yes]<-z_yes
  zmat[no]<-z_no
  zmat[na]<-z_na
  
  ###alpha and mu###
  alpha_beta<-cbind(1,beta)
  po_sig_mual<-chol2inv(chol((crossprod(alpha_beta)+sig_mu_alpha_inv)))
  
  for(j in 1:nc){
    po_mean_mual<-po_sig_mual%*%(crossprod(alpha_beta,zmat[,j])+part2)
    up_mu_alpha<-rmvn(1,po_mean_mual,po_sig_mual)
    mu_alpha_mat[j,]<-up_mu_alpha
  }
  
  mu<-mu_alpha_mat[,1]
  alpha<-mu_alpha_mat[,-1]
  ###beta###
  tt<-t(zmat)-mu
  var_beta<-chol2inv(chol((crossprod(alpha)+p_mat)))
  po_mean_beta<-var_beta%*%crossprod(alpha,tt)
  beta<-t(apply(po_mean_beta,2,douche))
  
  mean_mat<-t(tcrossprod(alpha,beta) + mu)
  yes_mean<-mean_mat[yes]
  no_mean<-mean_mat[no]
  na_mean<-mean_mat[na]
  
  if(i %in% seq(1,iter,50)){
    cat("\rProgress: ",i,"/",iter)
    y_hat = predict_ymat(alpha,beta,mu)
    ymat[na] = y_hat[na]
    
    # print(paste('trainning accuracy is',length(which(y_hat[no_na]==ymat[no_na]))/len_no_na))
    # if(test==T){
    #   print(paste('test accuracy is',length(which(y_hat[na_test]==y_test))/len_test))
    # }
  }
  
  if(i>burnin){
    if(hn==116){
      if(mean(beta[squad[1:4],1])>0){
        beta_squad_dim_1[,jj] = rank(-beta[,1])[squad]
        beta_squad_dim_2[,jj] = rank(-beta[,2])[squad]
      }else{
        beta_squad_dim_1[,jj] = rank(beta[,1])[squad]
        beta_squad_dim_2[,jj] = rank(beta[,2])[squad]
      }
    }
    beta_cov = diag(cov(beta))
    cov_master[jj] =  beta_cov[1]>beta_cov[2]
    beta_rs = beta_rs + beta
    mu_rs = mu_rs + mu
    alpha_rs = alpha_rs + alpha
    kobe = dbinom(ymat,1,pnorm(mean_mat),log = T)
    likeli_chain[jj] = sum(kobe)
    pos_pred<-pos_pred+exp(kobe)
    pos_pred2<-pos_pred2+kobe
    pos_pred3<-pos_pred3+kobe^2
    jj = jj + 1
  }
  
}

sum(cov_master)==n_pos

##########################################################################################
waic_compute(n_pos,pos_red,pos_pred2,pos_pred3,no_na)
plot(likeli_chain)
beta_rs<-beta_rs/n_pos
mu_rs<-mu_rs/n_pos
alpha_rs<-alpha_rs/n_pos
y_hat = predict_ymat(alpha,beta,mu)

y_hat = predict_ymat(alpha_rs,beta_rs,mu_rs)
print(paste('trainning accuracy is',length(which(y_hat[no_na]==ymat[no_na]))/len_no_na))
print(paste('test accuracy is',length(which(y_hat[na_test]==y_test))/len_test))

beta_var = apply(beta_rs,2,var)
alpha_var = apply(alpha_rs,2,var)

pred_vec = pred_vec_test = rep(0,p)
for(ppp in 1:p){
    beta_temp = cbind(beta_rs[,1:ppp],matrix(0,nr,p-ppp))
    alpha_temp = cbind(alpha_rs[,1:ppp],matrix(0,nc,p-ppp))
    
    y_hat = predict_ymat(alpha_temp,beta_temp,mu_rs)
    pred_vec[ppp] = length(which(y_hat[no_na]==ymat[no_na]))/len_no_na
    pred_vec_test[ppp] = length(which(y_hat[na_test]==y_test))/len_test
}
if(alpha_penalty==T){
  save(file=paste0(h_s,hn,"_beta_nd_mean_scaled.Rdata"),beta_rs)
  save(file=paste0(h_s,hn,"_alpha_nd_mean_scaled.Rdata"),alpha_rs)
  save(file=paste0(h_s,hn,"_mu_mean_scaled.Rdata"),mu_rs)
}else{
  save(file=paste0(h_s,hn,"_beta_nd_mean.Rdata"),beta_rs)
  save(file=paste0(h_s,hn,"_alpha_nd_mean.Rdata"),alpha_rs)
  save(file=paste0(h_s,hn,"_mu_mean.Rdata"),mu_rs)
}

windows()
par(mfrow=c(1,2))
par(mar=c(7,7,5,5))
plot(1:p,pred_vec,type='l',ylim=c(0.8,1),cex.lab=2,ylab='training accuracy',xlab='dim')
abline(h=max(pred_vec),col='red')
plot(1:p,pred_vec_test,type='l',ylim=c(0.8,1),cex.lab=2,ylab='test accuracy',xlab='dim')
abline(h=max(pred_vec_test),col='red')
# pred_vec_test
# max(pred_vec_test)*0.99
# aa = abs(pred_vec_test-max(pred_vec_test)*0.99)
# which(aa==min(aa))
# windows()
# par(mfrow=c(3,3))
# 
# for(ppp in 1:p){
#   hist(beta_rs[,ppp],xlim=c(-3,3))
# # }
# windows()
# par(mfrow=c(1,2))
# plot(1:p,beta_var)
# plot(1:p,alpha_var)
e_val =  eigen(cov(beta_rs))$values
cumsum(e_val/sum(e_val))
#########################



rank_squad_dim1 = round(t(apply(beta_squad_dim_1,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rank_squad_dim2 = round(t(apply(beta_squad_dim_2,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)

rownames(rank_squad_dim1) = rownames(rank_squad_dim2) = pol[squad]
save(file='rank_squad_dim1',rank_squad_dim1)

save(file='rank_squad_dim2',rank_squad_dim2)



library(xtable)
xtable(rank_squad_dim1)

xtable(rank_squad_dim2)












