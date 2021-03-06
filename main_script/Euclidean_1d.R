
args = commandArgs(trailingOnly = T)
la = length(args)
if(la==0){
  stop('No argument provided, please input 3 arguments')
}else if(la==3){
  iter = as.numeric(args[1]) ### iterations to run 
  burnin = as.numeric(args[2]) ### burnin period of mcmc
  hn = as.numeric(args[3]) ### the number of house data to be analyzed 
}else{
  stop('Argument length not correct, please input 3 arguments')
}
### loading functions
source(file="./source/read_kh2.R")
source(file="./source/ymat_spit.R")

required_package = c('mvnfast','truncnorm','wnominate','pscl')
check_package = sum(unlist(lapply(required_package, require, character.only = TRUE)))==4
if(check_package ==F){
  install.packages(required_package,repos = "http://cran.us.r-project.org")
  lapply(required_package, require, character.only = TRUE)
}
### set seed number
seed = 2016
set.seed(seed)
###load and preprocess data
out = ymat_spit(hn=hn)
if(hn==116){
  print('Only session 1 (first 700 votes) of the 116 House data will be analyzed')
}
ymat = out[[1]]
pol = out[[2]]
rm(out)

###Set hyperparameter and initialize variables
n_pos = iter - burnin 
no_na = which(is.na(ymat)==F)
len_no_na = length(no_na)
nr = nrow(ymat)
nc = ncol(ymat)
a_mu = a_alpha = rep(0,2)
sig_mu_alpha = diag(2)
sig_mu_alpha_inv = solve(sig_mu_alpha)
sig_beta = diag(nr)
sig_beta_inv = solve(sig_beta)
beta_master = matrix(0,nr,n_pos)

zmat = matrix(0,nr,nc)
yes = which(ymat==1)
yes_l = length(yes)
no = which(ymat==0)
no_l = length(no)
na = which(is.na(ymat)==T)
na_l = length(na)

yes_mean = rep(0,yes_l)
no_mean = rep(0,no_l)
na_mean = rep(0,na_l)

beta = rep(0,nr)
beta[grep('\\(D',pol)] = -1
beta[grep('\\(R',pol)] = 1

part2 = sig_mu_alpha_inv%*%a_mu

mu_alpha_mat = matrix(0,nc,2)
pos_pred = pos_pred2 = pos_pred3 = matrix(0,nr,nc)
p_index = 1
y_hat = matrix(0,nr,nc)
likeli_chain = rep(0,n_pos)

### run starts
for(i in 1:iter){
  if(i %in% seq(1,iter,100)){
    cat("\rProgress: ",i,"/",iter)
  }
  ###sample latent Z###
  z_yes = rtruncnorm(1,a=0,b=Inf,mean=yes_mean,sd=1)
  z_no = rtruncnorm(1,a=-Inf,b=0,no_mean,sd=1)
  z_na = rnorm(na_l,na_mean,sd=1)
  zmat[yes] = z_yes
  zmat[no] = z_no
  zmat[na] = z_na
  
  ###sample \alpha and \mu###
  alpha_beta = cbind(1,beta)
  po_sig_mual = chol2inv(chol((crossprod(alpha_beta)+sig_mu_alpha_inv)))
  for(j in 1:nc){
    po_mean_mual = po_sig_mual%*%(crossprod(alpha_beta,zmat[,j])+part2)
    up_mu_alpha = rmvn(1,po_mean_mual,po_sig_mual)
    mu_alpha_mat[j,] = up_mu_alpha
  }
  mu = mu_alpha_mat[,1]

  alpha = mu_alpha_mat[,2]
  ### sample \beta
  tt = t(zmat)-mu
  var_beta = c(chol2inv(chol((crossprod(alpha)+1))))
  po_mean_beta = var_beta%*%crossprod(alpha,tt)
  beta = c(rmvn(1,po_mean_beta,diag(var_beta,nr)))
  
  mean_mat = t(tcrossprod(alpha,beta) + mu)
  
  yes_mean = mean_mat[yes]
  no_mean = mean_mat[no]
  na_mean = mean_mat[na]
  ###impute missing value
  if(i %in% seq(1,iter,25)){
    y_p = pnorm(mean_mat)
    index_yes = which(y_p>0.5)
    y_hat[index_yes] = 1
    y_hat[-index_yes] = 0
    ymat[na] = y_hat[na]
  }
  ### record paratmer after burnin and compute waic using running sums
  if(i>burnin){
    beta_master[,p_index] = beta
    p_s1 = pnorm(t(tcrossprod(alpha,beta) + mu))
    kobe = dbinom(ymat,1,p_s1,log = T)
    likeli_chain[p_index] = sum(kobe)
    pos_pred = pos_pred + exp(kobe)
    pos_pred2 = pos_pred2 + kobe
    pos_pred3 = pos_pred3 + kobe^2
    p_index = p_index + 1
  }
  
}
###compute waic with -1 scaling, recall that -2 corresponds to deviance scaling###
waic = -waic_compute(n_pos,pos_pred,pos_pred2,pos_pred3,no_na)
print(waic)

save(file=paste0('H',hn,"_beta_master_eu.Rdata"),beta_master)

