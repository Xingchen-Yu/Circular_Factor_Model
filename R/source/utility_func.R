
## group WAIC
waic_compute_group = function(nnn,pos_pred2){
  lpd = sum(apply(pos_pred2,1,function(x) logSumExp(x) - log(nnn)))
  penalty = 2 * (lpd - sum(apply(pos_pred2,1,mean)))
  penalty_va = sum(apply(pos_pred2,1,var))
  return(c(lpd,penalty, penalty_va))
}
## WAIC
waic_compute = function(n_pos,pos_pred,pos_pred2,pos_pred3,no_na){
  pos_pred_master = pos_pred[no_na]/n_pos
  lpd = sum(log(pos_pred_master))
  va = sum((pos_pred3[no_na]/n_pos-(pos_pred2[no_na]/n_pos)^2))*n_pos/(n_pos-1)
  waic_spherical = lpd-va
}
update_tsig = function(upper,lower,t_sig,check,nnn){
  l_s = which(check<=lower)
  u_s = which(check>=upper)
  
  t_sig[l_s] = t_sig[l_s] * 0.95
  t_sig[u_s] = t_sig[u_s] * 1.05
  return(list(t_sig,length(c(l_s,u_s))/nnn))
}

z_to_x<-function(z){
  pi2 = pi^2
  pi22 = 2 * pi2
  x<-(z+pi2)/pi22
  return(x)
}
ymat_sim<-function(nr,nc,kappa_a,ccc,omega_a,omega_b){
  omega_master = rgamma(nr,omega_a,omega_b)
  beta = sapply(omega_master,function(x) circular::rvonmises(1,pi,x)-pi)
  
  tau_no = runif(nc,-pi,pi)
  tau_yes = runif(nc,-pi,pi)
  
  kappa = rgamma(nc,kappa_a,rgamma(nc,1,ccc))
  kappa_mat = matrix(kappa,nr,nc,byrow=T)
  asset<-acos(cos(t(outer(tau_no,beta,FUN="-"))))^2-acos(cos(t(outer(tau_yes,beta,FUN="-"))))^2
  x<-z_to_x(asset)
  sup = pbeta(x,kappa_mat,kappa_mat)
  sup_log_lower = pbeta(x,kappa_mat,kappa_mat,log.p = T,lower.tail = F)
  ymat = matrix(0,nr,nc)
  for(i in 1:nr){
    for(j in 1:nc){
      p_success = sup[i,j]
      p_success_c = exp(sup_log_lower[i,j])
      ymat[i,j] = sample(c(1,0),1,prob=c(p_success,p_success_c))
    }
  }
  
  return(list(ymat,beta,tau_no,tau_yes,kappa))
}

ymat_sim_eu = function(nr,nc,model_index){
  y_hat = matrix(0,nr,nc)
  mu = rnorm(nc)
  if(model_index==1){
    alpha = rnorm(nc)
    beta = rnorm(nr)
  }else{
    alpha = rmvn(nc,rep(0,model_index),diag(model_index))
    beta = rmvn(nr,rep(0,model_index),diag(model_index))
  }
  mean_mat = t(tcrossprod(alpha,beta) + mu)
  y_p = pnorm(mean_mat)
  y_p_log_lower = pnorm(mean_mat,lower.tail = F,log.p = T)
  for(i in 1:nr){
    for(j in 1:nc){
      yp = y_p[i,j]
      yp_lower = exp(y_p_log_lower[i,j])
      y_hat[i,j] = sample(c(1,0),1,prob=c(yp,yp_lower))
    }
  }
  return(y_hat)
}
# set.seed(2021)###toy data set in package
# out = ymat_sim(20,10,1,50,1,1/10)