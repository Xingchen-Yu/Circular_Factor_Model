args = commandArgs(trailingOnly = T)
la = length(args)
if(la==0){
  stop('No argument provide, please provide argument')
}else if(la==4){
  iter = as.numeric(args[1])
  burnin = as.numeric(args[2])
  core = as.numeric(args[3])
  hn = as.numeric(args[4])
}else{
  stop('Arguments length not correct, accept 4 arguments')
}

required_package = c('Rcpp','snowfall','wnominate','rlecuyer','RcppArmadillo')
lapply(required_package, require, character.only = TRUE)

######################
source(file="../source/read_kh2.R")
source(file="../source/geodesic_snowfall_rcpp_wrapper_v3_delta_jitter_bracket.R")
source(file="../source/ymat_spit.R")

if(hn==116){
  cluster_seed = 888
  seed = 8888
  b_range = c(0.01,0.04)
}else if(hn==112){
  cluster_seed = 1990
  seed = 2080
  b_range = c(0.005,0.04)
}else{
  cluster_seed = 1234
  seed = 4321
  b_range = c(0.01,0.04)
}
##################
sourceCpp(code = RcppCode)
sfInit(parallel=TRUE,cpus=core)
sfClusterSetupRNG( type="RNGstream",seed=cluster_seed)
sfLibrary("Rcpp", character.only=TRUE)
sfLibrary("RcppArmadillo", character.only=TRUE)
##################
set.seed(seed)
n_pos = iter - burnin 
start_jitter = 51
skip = 50
yn_range = c(0.01,0.105)
l_range = c(1,10)
###hyperparameter######
a = 1
b = 1/10
mu = 0
ccc_a = 1
ccc_b = 25
kappa_a = 1
omega_sd = 0.1
#######data###########
out = ymat_spit(hn=hn)
ymat = out[[1]]
pol = out[[2]]
rm(out)

nr = nrow(ymat)
nc = ncol(ymat)
dem = grep("\\(D",pol)
gop = grep("\\(R",pol)
ind = grep("\\(I",pol)
##################################
t_sig = rep(0.5,nc)
omega = a/b
ccc = ccc_a/ccc_b

kappa = rep(0.1,nc)
omega_ini = 0

beta = rep(0,nr)
beta[dem] = runif(length(dem),-pi/2,0)
beta[gop] = runif(length(gop),0,pi/2)
tau_no = runif(nc,-pi,pi)
tau_yes = runif(nc,-pi,pi)

leap = sample(l_range[1]:l_range[2],nr,replace=T)
leap_tau = sample(l_range[1]:l_range[2],nc,replace=T)

delta = runif(nr,b_range[1],b_range[2])
delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])

delta2 = delta/2
delta2_yes = delta_yes/2
delta2_no = delta_no/2
################################
lol = c(sin(mu),cos(mu))
cbeta_prior = lol * omega
c_alpha = kappa_a*nc+ccc_a
##################################
nr_par = round(seq(0,nr,length.out = core+1))
nc_par = round(seq(0,nc,length.out = core+1))
sfExportAll(except=c("dem","gop",'ind','omega_ini'))
##################################
na = which(is.na(ymat==T)) - 1
len_na = length(na)
na.position = which(is.na(ymat)==T, arr.ind = T)
i_index = as.numeric(na.position[,1]) - 1
j_index = as.numeric(na.position[,2]) - 1
##################################
col_stream = rep('red',nr)
col_stream[dem] = 'blue'
col_stream[ind] = 'green'
##################################
if(n_pos>0){
  pos_pred = pos_pred2 = pos_pred3 = matrix(0,nr,nc)
  no_na = which(is.na(ymat)==F)
  kappa_master = matrix(0,nc,n_pos)
  beta_master = matrix(0,nr,n_pos)
  omega_master = ccc_master = rep(0,n_pos)
}
posterior_chain = likeli_chain = rep(0,iter)

wrapper_beta = function(t){
  update_beta(t,nr_par,delta,delta2,leap,nc,omega,cbeta_prior,beta,tau_yes,tau_no,kappa,ymat)
}

wrapper_yes = function(t){
  update_tau_yes(t,nc_par,delta_yes,delta2_yes,leap_tau,nr,beta,tau_yes,tau_no,kappa,ymat)
}

wrapper_no = function(t){
  update_tau_no(t,nc_par,delta_no,delta2_no,leap_tau,nr,beta,tau_yes,tau_no,kappa,ymat)
}

wrapper_kappa = function(t){
  update_kappa(t,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat,kappa_a,ccc,t_sig)
}

wrapper_waic = function(t){
  waic_cpp(t,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat)
}
wrapper_predict = function(t){
  predict(t,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat)
}
wrapper_predict_mean = function(t){
  predict(t,nc_par,nr,beta_mean,yes_mean,no_mean,kappa_mean,ymat)
}
sfClusterEval(sourceCpp(code = RcppCode))
core_1 = core - 1
node = 0:core_1
j = 1
omega_ratio = 0
beta_accept_rs_all = rep(0,nr)
kappa_accept_rs = rep(0,nc)
kappa_accept_rs_all = no_accept_rs_all = yes_accept_rs_all = rep(0,nc)

for(i in 1:iter){
  if(i %in% seq(start_jitter,iter,skip)){
    
    leap = sample(l_range[1]:l_range[2],nr,replace=T)
    leap_tau = sample(l_range[1]:l_range[2],nc,replace=T)
    
    delta = runif(nr,b_range[1],b_range[2])
    delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])
    #######################################
    delta2 = delta/2
    delta2_yes = delta_yes/2
    delta2_no = delta_no/2
    
    if(i<2000 & continue==F){
      ks = kappa_accept_rs/skip
      kappa_skip = min(ks)
      
      out = update_tsig(0.6,0.3,t_sig,ks,nc)
      t_sig = out[[1]]
      kappa_mod = out[[2]]
      
      kappa_accept_rs = rep(0,nc)
      print(paste0('percent kappa changed ',kappa_mod))
      print(paste0('min kappa acceptance ',kappa_skip))
    }
    
    
    sfExport("delta",'delta2','delta_yes','delta2_yes','delta_no','delta2_no','t_sig','leap','leap_tau')
  }
  
  
  if(i %in% seq(1,iter,20)){
    ymat = impute_NA(na, i_index, j_index, ymat, tau_yes, tau_no, beta, kappa, len_na)
    sfExport('ymat')
  }
  ################################################################################

  out = sfLapply(node,wrapper_kappa)
  kappa = unlist(lapply(out,"[[",1))
  sfExport("kappa")
  
  ccc = rgamma(1,c_alpha,ccc_b+sum(kappa))
  sfExport("ccc")
  
  kappa_ratio = sum(sapply(out,"[[",2))/nc
  haha = unlist(lapply(out,"[[",3))
  kappa_accept_rs = kappa_accept_rs  + haha
  kappa_accept_rs_all = kappa_accept_rs_all  + haha
  ################################################################################
  out = sfLapply(node,wrapper_beta)
  beta = unlist(lapply(out,"[[",1))
  beta_ratio = sum(sapply(out,"[[",2))/nr
  
  haha = unlist(lapply(out,"[[",3))
  beta_accept_rs_all = beta_accept_rs_all  + haha
  sfExport("beta")
  ################################################################################
  out = sfLapply(node,wrapper_yes)
  tau_yes = unlist(lapply(out,"[[",1))
  yes_ratio = sum(sapply(out,"[[",2))/nc
  
  haha = unlist(lapply(out,"[[",3))
  yes_accept_rs_all = yes_accept_rs_all  + haha
  sfExport("tau_yes")
  ################################################################################
  out = sfLapply(node,wrapper_no)
  tau_no = unlist(lapply(out,"[[",1))
  no_ratio = sum(sapply(out,"[[",2))/nc
  
  haha = unlist(lapply(out,"[[",3))
  no_accept_rs_all = no_accept_rs_all  + haha
  sfExport("tau_no")
  ################################################################################
  out = update_omega(omega,beta,nr,a,b,omega_sd)
  omega = out[[1]]
  omega_ratio = omega_ratio + out[[2]]
  cbeta_prior = lol * omega
  sfExport('omega',"cbeta_prior")
  
  waic_out = sfLapply(node,wrapper_waic)
  temp = do.call("cbind",lapply(waic_out,"[[",1))
  sum_temp = sum(temp)
  likeli_chain[i] = sum_temp

  if(i %in% seq(1,iter,50)){
    cat("\rProgress: ",i,"/",iter)
    print(paste('beta ar is',round(beta_ratio,digits = 2)))
    print(paste('yes ar is',round(yes_ratio,digits = 2)))
    print(paste('no ar is',round(no_ratio,digits = 2)))
    print(paste('kappa ar is',round(kappa_ratio,digits = 2)))
    print(paste('omega ar is',round(omega_ratio/i,digits = 2)))

    # print(paste('min kappa ar is',min(kappa_accept_rs_all/i)))
    # print(paste('min beta ar is',min(beta_accept_rs_all/i)))
    # print(paste('min yes ar is',min(yes_accept_rs_all/i)))
    # print(paste('min no ar is',min(no_accept_rs_all/i)))

    print(paste('loglikelihood is ',round(sum_temp,digits = 0)))
  }
  if(i>burnin){
    beta_master[,j] = beta
    omega_master[j] = omega
    kappa_master[,j] = kappa

    pos_pred = pos_pred + do.call("cbind",lapply(waic_out,"[[",2))
    pos_pred2 = pos_pred2 + temp
    pos_pred3 = pos_pred3 + do.call("cbind",lapply(waic_out,"[[",3))
    j = j + 1
  }
}
waic_spherical = waic_compute(n_pos,pos_pred,pos_pred2,pos_pred3,no_na)

nnn = n_pos
pos_pred_master = pos_pred[no_na]/nnn
lpd = sum(log(pos_pred_master))
va = sum((pos_pred3[no_na]/nnn-(pos_pred2[no_na]/nnn)^2))*nnn/(nnn-1)
waic_spherical = lpd-va
print(waic_spherical)

# save.image(file=paste0('H',hn,"_workspace.Rdata"))
