
args = commandArgs(trailingOnly = T)
la = length(args)
if(la==0){
  stop('No argument provided, please input 3 arguments')
}else if(la==3){
  iter = as.numeric(args[1])
  burnin = as.numeric(args[2])
  hn = as.numeric(args[3])
}else{
  stop('Argument length not correct, please input 3 arguments')
}

source(file="./source/read_kh2.R")
source(file="./source/ymat_spit.R")

required_package = c('mvnfast','truncnorm','wnominate','pscl')
check_package = sum(unlist(lapply(required_package, require, character.only = TRUE)))==4
if(check_package ==F){
  install.packages(required_package,repos = "http://cran.us.r-project.org")
  lapply(required_package, require, character.only = TRUE)
}

####################################
seed = 2020
####################################  
p = 2
n_pos = iter - burnin 
likeli_chain = rep(0,n_pos)
####################################
p_mat = diag((1:p)^2)

out = ymat_spit(hn=hn)
if(hn==116){
  squad = c(195,216,210,270)
  print(pol[squad])
  beta_squad_dim_2 = beta_squad_dim_1 = matrix(0,length(squad),n_pos)
  print('Only session 1 (first 700 votes) of the 116 House data will be analyzed')
}
ymat = out[[1]]
pol = out[[2]]
rm(out)

nr = nrow(ymat)
nc = ncol(ymat)
a_mu = rep(0,p+1)

sig_mu_alpha = diag(p+1)
sig_mu_alpha_inv = solve(sig_mu_alpha)

beta_rs = matrix(0,nr,p)
mu_rs = rep(0,nc)
alpha_rs = matrix(0,nc,p)
zmat = matrix(0,nr,nc)

yes = which(ymat==1)
yes_l = length(yes)
no = which(ymat==0)
no_l = length(no)
na = which(is.na(ymat)==T)

no_na = which(is.na(ymat)==F)
len_no_na = length(no_na)
na_l = length(na)

yes_mean = rep(0,yes_l)
no_mean = rep(0,no_l)
na_mean = rep(0,na_l)

beta = matrix(0,nr,p)
pos_pred = pos_pred2 = pos_pred3 = matrix(0,nr,nc)
cov_master =rep(NA,n_pos)
part2 = sig_mu_alpha_inv%*%a_mu
##
mu_alpha_mat = matrix(0,nc,p+1)
sample_rmvn = function(kobe){
  sup = rmvn(1,c(kobe),var_beta)
  return(sup)
}
y_hat = matrix(0,nr,nc)
jj = 1
for(i in 1:iter){
  
  ###latent###
  z_yes = rtruncnorm(1,a=0,b=Inf,mean=yes_mean,sd=1)
  z_no = rtruncnorm(1,a=-Inf,b=0,no_mean,sd=1)
  z_na = rnorm(na_l,na_mean,sd=1)
  zmat[yes] = z_yes
  zmat[no] = z_no
  zmat[na] = z_na
  
  ###alpha and mu###
  alpha_beta = cbind(1,beta)
  po_sig_mual = chol2inv(chol((crossprod(alpha_beta)+sig_mu_alpha_inv)))
  
  for(j in 1:nc){
    po_mean_mual = po_sig_mual%*%(crossprod(alpha_beta,zmat[,j])+part2)
    up_mu_alpha = rmvn(1,po_mean_mual,po_sig_mual)
    mu_alpha_mat[j,] = up_mu_alpha
  }
  
  mu = mu_alpha_mat[,1]
  alpha = mu_alpha_mat[,-1]
  ###beta###
  tt = t(zmat)-mu
  var_beta = chol2inv(chol((crossprod(alpha)+p_mat)))
  po_mean_beta = var_beta%*%crossprod(alpha,tt)
  beta = t(apply(po_mean_beta,2,sample_rmvn))
  
  mean_mat = t(tcrossprod(alpha,beta) + mu)
  yes_mean = mean_mat[yes]
  no_mean = mean_mat[no]
  na_mean = mean_mat[na]
  ###impute missing value####
  if(i %in% seq(1,iter,50)){
    cat("\rProgress: ",i,"/",iter)
    y_p = pnorm(mean_mat)
    index_yes = which(y_p>0.5)
    y_hat[index_yes] = 1
    y_hat[-index_yes] = 0
    ymat[na] = y_hat[na]
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
    pos_pred = pos_pred+exp(kobe)
    pos_pred2 = pos_pred2+kobe
    pos_pred3 = pos_pred3+kobe^2
    jj = jj + 1
  }
  
}
###compute waic with -1 scaling, recall that -2 corresponds to deviance scaling###
waic_spherical = -waic_compute(n_pos,pos_pred,pos_pred2,pos_pred3,no_na)
print(waic_spherical)

rank_squad_dim1 = round(t(apply(beta_squad_dim_1,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rank_squad_dim2 = round(t(apply(beta_squad_dim_2,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))),digits=0)
rownames(rank_squad_dim1) = rownames(rank_squad_dim2) = pol[squad]
save(file='rank_squad_dim1',rank_squad_dim1)
save(file='rank_squad_dim2',rank_squad_dim2)
print(rank_squad_dim1)
