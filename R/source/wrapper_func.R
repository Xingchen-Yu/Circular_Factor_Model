wrapper_kappa = function(t){
  update_kappa(t,nc_par,nr,beta_i,tau_yes,tau_no,kappa_j,ymat,kappa_a,ccc,t_sig)
}

wrapper_beta = function(t){
  update_beta(t,nr_par,delta,delta2,leap,nc,omega,cbeta_prior,beta_i,mu,tau_yes,tau_no,kappa_j,ymat)
}

wrapper_yes = function(t){
  update_tau_yes(t,nc_par,delta_yes,delta2_yes,leap_tau,nr,beta_i,tau_yes,tau_no,kappa_j,ymat)
}

wrapper_no = function(t){
  update_tau_no(t,nc_par,delta_no,delta2_no,leap_tau,nr,beta_i,tau_yes,tau_no,kappa_j,ymat)
}

wrapper_waic = function(t){
  waic_cpp(t,nc_par,nr,beta_i,tau_yes,tau_no,kappa_j,ymat)
}
