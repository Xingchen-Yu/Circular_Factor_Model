# args = commandArgs(trailingOnly = T)
# la = length(args)
# if(la==0){
#   stop('No argument provided, please input 4 arguments')
# }else if(la==4){
#   iter = as.numeric(args[1]) ### total iterations to run
#   burnin = as.numeric(args[2]) ### burnin period of mcmc
#   core = as.numeric(args[3]) ### number of cores
#   hn = as.numeric(args[4]) ### the number of house data to be analyzed 
# }else{
#   stop('Argument length not correct, please input 4 arguments')
# }

iter = 21000
burnin = 1000
core = 10
hn = 101
house = T
h_s = 'H'
#### checking and installing required packages###
required_package = c('Rcpp','snowfall','wnominate','rlecuyer','RcppArmadillo','pscl','matrixStats')
check_package = sum(unlist(lapply(required_package, require, character.only = TRUE)))==7
if(check_package ==F){
  install.packages(required_package,repos = "http://cran.us.r-project.org")
  lapply(required_package, require, character.only = TRUE)
}
### loading functions
source(file="./source/read_kh2.R")
source(file="./source/ymat_spit.R")
source(file="./source/circular_factor_model_functions.R")
#### setting both cluster seed and local seed for reproducibility 

WAIC_group = T ###grouped waic calculation
continue = T ###
if(hn==116){
  cluster_seed = 888
  seed = 8888
  b_range = c(0.01,0.04) 
}else if(hn==112){
  cluster_seed = 1990
  seed = 2080
  b_range = c(0.005,0.04)
}else{
  cluster_seed = 12345
  seed = 43215
  b_range = c(0.005,0.04)
}
if(continue==T){
  cluster_seed = cluster_seed + as.numeric(continue)
  seed = seed + as.numeric(continue)
}
###Set up parallel environment 
sourceCpp(code = RcppCode)
sfInit(parallel=TRUE,cpus=core)
sfClusterSetupRNG( type="RNGstream",seed=cluster_seed)
sfLibrary("Rcpp", character.only=TRUE)
sfLibrary("RcppArmadillo", character.only=TRUE)
###load and preprocess data
out = ymat_spit(hn=hn)
if(hn==116){
  print('Only session 1 (first 700 votes) of the 116 House data will be analyzed')
}
ymat = out[[1]]
pol = out[[2]]
rm(out)
###Set hyperparameter and initialize variables
set.seed(seed)
n_pos = iter - burnin 
start_jitter = 51
skip = 50
yn_range = c(0.01,0.105)
l_range = c(1,10)

a = 1
b = 1/10
mu = 0
ccc_a = 1
ccc_b = 25
kappa_a = 1
omega_sd = 0.1

nr = nrow(ymat)
nc = ncol(ymat)
dem = grep("\\(D",pol)
gop = grep("\\(R",pol)
ind = grep("\\(I",pol)

if(continue==F){
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
}else{
  load(file=paste0('./continue/',h_s,hn,"_beta_start.Rdata"),verbose = T)
  load(file=paste0('./continue/',h_s,hn,"_tau_yes_start.Rdata"),verbose = T)
  load(file=paste0('./continue/',h_s,hn,"_tau_no_start.Rdata"),verbose = T)
  load(file=paste0('./continue/',h_s,hn,"_kappa_start.Rdata"),verbose = T)
  load(file=paste0('./continue/',h_s,hn,"_ccc_start.Rdata"),verbose = T)
  load(file=paste0('./continue/',h_s,hn,"_omega_start.Rdata"),verbose = T)
  load(file=paste0('./continue/',h_s,hn,"_tsig_start.Rdata"),verbose = T)
  
  # load(file=paste0(h_s,hn,"_beta_start.Rdata"),verbose = T)
  # load(file=paste0(h_s,hn,"_tau_yes_start.Rdata"),verbose = T)
  # load(file=paste0(h_s,hn,"_tau_no_start.Rdata"),verbose = T)
  # load(file=paste0(h_s,hn,"_kappa_start.Rdata"),verbose = T)
  # load(file=paste0(h_s,hn,"_ccc_start.Rdata"),verbose = T)
  # load(file=paste0(h_s,hn,"_omega_start.Rdata"),verbose = T)
  # load(file=paste0(h_s,hn,"_tsig_start.Rdata"),verbose = T)
}

leap = sample(l_range[1]:l_range[2],nr,replace=T)
leap_tau = sample(l_range[1]:l_range[2],nc,replace=T)

delta = runif(nr,b_range[1],b_range[2])
delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])

delta2 = delta/2
delta2_yes = delta_yes/2
delta2_no = delta_no/2

lol = c(sin(mu),cos(mu))
cbeta_prior = lol * omega
c_alpha = kappa_a*nc+ccc_a
######parallel miscellaneous############################
nr_par = round(seq(0,nr,length.out = core+1))
nc_par = round(seq(0,nc,length.out = core+1))
sfExportAll(except=c("dem","gop",'ind','omega_ini'))

na = which(is.na(ymat==T)) - 1
len_na = length(na)
na.position = which(is.na(ymat)==T, arr.ind = T)
i_index = as.numeric(na.position[,1]) - 1
j_index = as.numeric(na.position[,2]) - 1
### initialize output vecotr/matrix
if(n_pos>0){
  if(WAIC_group==T){
    pos_pred_group = matrix(0,nr,n_pos)
  }
  pos_pred = pos_pred2 = pos_pred3 = matrix(0,nr,nc)
  
  no_na = which(is.na(ymat)==F)
  kappa_master = matrix(0,nc,n_pos)
  beta_master = matrix(0,nr,n_pos)
  omega_master = ccc_master = rep(0,n_pos)
}
posterior_chain = likeli_chain = rep(0,iter)

### Define wrapper function for parallel
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
### Run starts
for(i in 1:iter){
  
  if(i %in% seq(start_jitter,iter,skip)){
    #### step size and leap steps jittering
    leap = sample(l_range[1]:l_range[2],nr,replace=T)
    leap_tau = sample(l_range[1]:l_range[2],nc,replace=T)
    
    delta = runif(nr,b_range[1],b_range[2])
    delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])
    delta2 = delta/2
    delta2_yes = delta_yes/2
    delta2_no = delta_no/2
    #######################################
    ### tuning the proposal variance of the scale parameter \kappa during the initial 2000 iterations
    if(i<burnin){
      ks = kappa_accept_rs/skip
      kappa_skip = min(ks)
      ### targeting acceptance ratio between 0.3 and 0.6#####
      out = update_tsig(0.6,0.3,t_sig,ks,nc)
      t_sig = out[[1]]
      kappa_mod = out[[2]]
      
      kappa_accept_rs = rep(0,nc)
      print(paste0('percent kappa changed ',kappa_mod))
      print(paste0('min kappa acceptance ',kappa_skip))
    }
    
    
    sfExport("delta",'delta2','delta_yes','delta2_yes','delta_no','delta2_no','t_sig','leap','leap_tau')
  }
  
  ### impute the missing values
  if(i %in% seq(1,iter,20)){
    ymat = impute_NA(na, i_index, j_index, ymat, tau_yes, tau_no, beta, kappa, len_na)
    sfExport('ymat')
  }
  ################################################################################
  ### Update scale parameter \kappa
  out = sfLapply(node,wrapper_kappa)
  kappa = unlist(lapply(out,"[[",1))
  sfExport("kappa")
  ###update hyperprior parameter for \kappa
  ccc = rgamma(1,c_alpha,ccc_b+sum(kappa))
  sfExport("ccc")
  
  kappa_ratio = sum(sapply(out,"[[",2))/nc
  haha = unlist(lapply(out,"[[",3))
  kappa_accept_rs = kappa_accept_rs  + haha
  kappa_accept_rs_all = kappa_accept_rs_all  + haha
  ################################################################################
  ### update ideal points \beta_i's
  out = sfLapply(node,wrapper_beta)
  beta = unlist(lapply(out,"[[",1))
  beta_ratio = sum(sapply(out,"[[",2))/nr
  
  haha = unlist(lapply(out,"[[",3))
  beta_accept_rs_all = beta_accept_rs_all  + haha
  sfExport("beta")
  ################################################################################
  ### update \psi_j's
  out = sfLapply(node,wrapper_yes)
  tau_yes = unlist(lapply(out,"[[",1))
  yes_ratio = sum(sapply(out,"[[",2))/nc

  haha = unlist(lapply(out,"[[",3))
  yes_accept_rs_all = yes_accept_rs_all  + haha
  sfExport("tau_yes")
  ################################################################################
  ### update \zeta_j's
  out = sfLapply(node,wrapper_no)
  tau_no = unlist(lapply(out,"[[",1))
  no_ratio = sum(sapply(out,"[[",2))/nc
  
  haha = unlist(lapply(out,"[[",3))
  no_accept_rs_all = no_accept_rs_all  + haha
  sfExport("tau_no")
  ################################################################################
  ### update hyperparameter \omega for \beta_i's
  out = update_omega(omega,beta,nr,a,b,omega_sd)
  omega = out[[1]]
  omega_ratio = omega_ratio + out[[2]]
  cbeta_prior = lol * omega
  sfExport('omega',"cbeta_prior")
  
  ### compute joint loglikelihood
  waic_out = sfLapply(node,wrapper_waic)
  temp = do.call("cbind",lapply(waic_out,"[[",1))
  sum_temp = sum(temp[no_na]) ###exclude NA for DIC calculation
  likeli_chain[i] = sum_temp
  
  ## output average acceptance ratio and joint loglikelihood every 50 iterations
  if(i %in% seq(1,iter,50)){
    cat("\rProgress: ",i,"/",iter)
    print(paste('beta ar is',round(beta_ratio,digits = 2)))
    print(paste('yes ar is',round(yes_ratio,digits = 2)))
    print(paste('no ar is',round(no_ratio,digits = 2)))
    print(paste('kappa ar is',round(kappa_ratio,digits = 2)))
    print(paste('omega ar is',round(omega_ratio/i,digits = 2)))
    print(paste('loglikelihood is ',round(sum_temp,digits = 0)))
  }
  ### record paratmer after burnin and compute waic using running sums
  if(i>burnin){
    beta_master[,j] = beta
    omega_master[j] = omega
    
    if(WAIC_group==T){
      # temp[na+1] = 0 needs to include NA for WAIC computation
      temp2 = apply(temp,1,sum)
      pos_pred_group[,j] = temp2
    }
    
    pos_pred = pos_pred + do.call("cbind",lapply(waic_out,"[[",2))
    pos_pred2 = pos_pred2 + temp
    pos_pred3 = pos_pred3 + do.call("cbind",lapply(waic_out,"[[",3))
    
    j = j + 1
  }
}
waic_compute_new3 = function(nnn,pos_pred2){
  lpd = sum(apply(pos_pred2,1,function(x) logSumExp(x) - log(nnn)))
  penalty = 2 * (lpd - sum(apply(pos_pred2,1,mean)))
  penalty_va = sum(apply(pos_pred2,1,var))
  return(c(lpd,penalty, penalty_va))
}
waic_group = waic_compute_new3(n_pos,pos_pred_group)
rm(pos_pred_group)
### compute waic with -1 scaling, recall that -2 corresponds to deviance scaling
# waic_spherical = -waic_compute(n_pos,pos_pred,pos_pred2,pos_pred3,no_na)
# print(waic_spherical)

### save paramters for further analysis
save.image(file=paste0('H',hn,"_workspace_grouped_waic_dic.Rdata"))
# save(file=paste0('H',hn,"_beta_master_sph.Rdata"),beta_master)
# save(file=paste0('H',hn,"_pol.Rdata"),pol)

