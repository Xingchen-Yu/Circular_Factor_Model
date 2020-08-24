RcppCode = '

#include <RcppArmadillo.h>
#include <Rmath.h>
#include<functional>  
#include <boost/math/special_functions/bessel.hpp>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]

const double pi2 = pow(M_PI,2);
const double pi22 = 2 * pow(M_PI,2);
const double pi22_i = 1/pi22;
const double log_pi = log(M_PI);
const int MAXIT = 1000;
const double EPS = 3e-12;
const double FPMIN = 1e-30;
const double p = 2;
const arma::mat dp = arma::eye(p,p);
const double mu = 0;
const double log2_p_pi = log(2) + log(M_PI);

arma::vec s_to_e(double angle, arma::vec x){
  x(0) = sin(angle);
  x(1) = cos(angle);
  return x;
}

double betaContFrac(double a, double b, double x) {

  double qab = a + b;
  double qap = a + 1;
  double qam = a - 1;
  double c = 1;
  double d = 1 - qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  double h = d;
  int m;
  for (m = 1; m <= MAXIT; m++) {
    int m2 = 2 * m;
    double aa = m * (b-m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= (d * c);
    aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    double del = d*c;
    h *= del;
    if (fabs(del - 1) < EPS) break;
  }
  return h;
}

double betaInc(double a, double b, double x) {
  if (x == 0)
    return 0;
  else if (x == 1)
    return 1;
  else {
    double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
    if (x < (a+1) / (a+b+2))
      return exp(logBeta) * betaContFrac(a, b, x) / a;
    else
      return 1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b;
  }
}

double betaInc_log(double a, double b, double x) {
  if (x == 1)
    return 0;
  else {
    double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
    if (x < (a+1) / (a+b+2))
      return logBeta + log(betaContFrac(a,b,x))-log(a);
    else
      return log(1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b);
  }
}

double betaInc_log_lower(double a, double b, double x) {
  if (x == 0)
    return 0;
  else {
    double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
    if (1-x < (b+1) / (a+b+2))
      return logBeta + log(betaContFrac(b,a,1-x))-log(b);
    else
      return log(1 - exp(logBeta) * betaContFrac(a, b, x) / a);
  }
}

double dgamma_log (double x, double a, double b){
  return a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;
}

double likeli_omega(double omega, arma::vec beta,int nr,double a, double b){
  return  - nr * (log2_p_pi + log(boost::math::cyl_bessel_i(0,omega))) +
    omega * sum(cos(beta-mu)) + dgamma_log (omega, a, b);
}

double sb_c(double x){
  return (x + pi2)/pi22;
}
arma::vec sb_c_vec(arma::vec x){
  return (x + pi2)/pi22;
}


arma::vec Rf_rbinom_vec(arma::vec x,arma::vec shape1,const int n){
  arma::vec kobe = arma::zeros(n);
  double temp,shape;
  for(int i=0 ;i < n; i++){
    shape = shape1(i);
    temp = betaInc(shape,shape,x(i));
    //kobe(i) = Rf_rbinom(1,temp);
    if(as_scalar(arma::randu(1))<= temp){
      kobe(i) = 1;
    }else{
      kobe(i) = 0;
    }
  }
  return kobe;
}
// [[Rcpp::export]]
arma::mat impute_NA(arma::uvec na, arma::uvec i_index, arma::uvec j_index, arma::mat ymat,
                    arma::vec tau_yes, arma::vec tau_no, arma::vec beta, arma::vec kappa, const int n){

  arma::vec beta_na = beta.rows(i_index);
  arma::vec yes_na = tau_yes.rows(j_index);
  arma::vec no_na = tau_no.rows(j_index);
  arma::vec kappa_na = kappa.rows(j_index);

  arma::vec impute = sb_c_vec(square(acos(cos(no_na-beta_na))) - square(acos(cos(yes_na-beta_na)))) ;
  ymat.elem(na) = Rf_rbinom_vec(impute, kappa_na,n);

  return ymat;
}


double likeli_beta (double omega, double beta, arma::vec tau_yes,
                    arma::vec tau_no,arma::vec kappa ,int nc,arma::rowvec ymat_row){

  const arma::vec asset = sb_c_vec(square(acos(cos(tau_no-beta))) - square(acos(cos(tau_yes-beta))));
  arma::vec avec = arma::zeros(nc);
  double shape;

  for(int i = 0 ;i < nc; i++){
    shape = kappa(i);
    if(ymat_row(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec) + omega * cos(beta-mu);
}


arma::vec gradient_beta(arma::vec cbeta_prior, double beta, arma::vec tau_yes, arma::vec tau_no,arma::vec kappa,
                        const int nc, arma::rowvec ymat_row){

 arma::vec temp_x = arma::zeros(p);
 arma::vec please = arma::zeros(p);

 double temp,shape,shape2,logbeta,check,beta_tau_yes,beta_tau_no,avec,buffer_yes,buffer_no,yes_temp,no_temp;
  for(int i = 0 ;i < nc; i++){
    shape = kappa(i);
    yes_temp = tau_yes(i);
    no_temp = tau_no(i);
    beta_tau_yes = yes_temp - beta;
    beta_tau_no = no_temp - beta;
    buffer_yes = acos(cos(beta_tau_yes));
    buffer_no = acos(cos(beta_tau_no));
    temp = sb_c(pow(buffer_no,2)-pow(buffer_yes,2));


    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_row(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += 2 * avec * (buffer_yes/fabs(sin(beta_tau_yes)) * s_to_e(yes_temp,temp_x) -
      buffer_no/fabs(sin(beta_tau_no)) * s_to_e(no_temp,temp_x));
  }

  return please + cbeta_prior;
}


double likeli_tau_yes(double tau_yes,arma::vec beta,double shape,const int nr,
                      arma::vec ymat_col,arma::vec temp_no){
  
  const arma::vec asset = sb_c_vec(temp_no - square(acos(cos(tau_yes-beta))));
  arma::vec avec = arma::zeros(nr);
  
  for(int i = 0 ;i < nr; i++){
    if(ymat_col(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec);
}


arma::vec gradient_tau_yes(double tau_yes,arma::vec beta,double shape,const int nr,
                           arma::vec ymat_col,arma::vec temp_no){
  
  arma::vec temp_x = arma::zeros(p);
  arma::vec please = arma::zeros(p);
  
  double temp,shape2,logbeta,check,beta_tau_yes,avec,buffer_yes,beta_temp;
  for(int i = 0 ;i < nr; i++){
    beta_temp = beta(i);
    beta_tau_yes = tau_yes - beta_temp;
    buffer_yes = acos(cos(beta_tau_yes));
    temp = sb_c(temp_no(i) - pow(buffer_yes,2));
    
    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_col(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += 2 * avec * buffer_yes/fabs(sin(beta_tau_yes)) * s_to_e(beta_temp,temp_x);
  }
  
  return please;
}

double likeli_kappa(int nr, arma::vec asset,double shape,double kappa_a, double ccc,arma::vec ymat_col){
  arma::vec avec = arma::zeros(nr);
  
  for(int i = 0 ;i < nr; i++){
    if(ymat_col(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec) + dgamma_log(shape,kappa_a,ccc);
}

double likeli_tau_no(double tau_no,arma::vec beta,double shape,const int nr,
                      arma::vec ymat_col,arma::vec temp_yes){
  
  const arma::vec asset = sb_c_vec(square(acos(cos(tau_no-beta))) - temp_yes);
  arma::vec avec = arma::zeros(nr);
  
  for(int i = 0 ;i < nr; i++){
    if(ymat_col(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec);
}


arma::vec gradient_tau_no(double tau_no,arma::vec beta,double shape,const int nr,
                           arma::vec ymat_col,arma::vec temp_yes){
  
  arma::vec temp_x = arma::zeros(p);
  arma::vec please = arma::zeros(p);
  
  double temp,shape2,logbeta,check,beta_tau_no,avec,buffer_no,beta_temp;
  for(int i = 0 ;i < nr; i++){
    beta_temp = beta(i);
    beta_tau_no = tau_no - beta_temp;
    buffer_no = acos(cos(beta_tau_no));
    temp = sb_c(pow(buffer_no,2) - temp_yes(i));
    
    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_col(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += - 2 * avec * buffer_no/fabs(sin(beta_tau_no)) * s_to_e(beta_temp,temp_x);
  }
  
  return please;
}

// [[Rcpp::export]]
List update_beta(int t, arma::vec nr_par, arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec,int nc, double omega,
                      arma::vec cbeta_prior, arma::vec beta,arma::vec tau_yes,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
  arma::vec nu = arma::zeros(p);
  arma::vec x_temp = arma::zeros(p);
  arma::vec x = arma::zeros(p);
  double count = 0;
  double alpha, ae, cosat, sinat, h, h_new, accept,beta_new,beta_prev;
  int start = nr_par(t);
  int end = nr_par(t+1);
  int e_s = end - start;
  arma::vec beta_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);
  
for(int i = start; i< end; i++,count++ ){
    double epsilon = epsilon_vec(i);
    double epsilon2 = epsilon2_vec(i);
    int leap = leap_vec(i);
    arma::rowvec ymat_row = ymat.row(i);
    beta_prev = beta_new = beta(i);
    x = s_to_e(beta_prev,x);
    nu = arma::randn(p);
    nu = (dp - x * x.t()) * nu;
  
    h = likeli_beta(omega,beta_prev,tau_yes,tau_no,kappa,nc,ymat_row) - as_scalar(0.5 * nu.t() * nu);
    for(int j = 0; j < leap; j++) {
      nu = nu + epsilon2 * gradient_beta(cbeta_prior, beta_new, tau_yes, tau_no,kappa, nc, ymat_row);
      nu = (dp - x * x.t()) * nu;

      alpha = norm(nu);
      ae = alpha * epsilon;
      cosat = cos(ae);
      sinat = sin(ae);
      x_temp = x;
      //geodesic flow//
      x = x * cosat + nu/alpha * sinat;
      nu = nu * cosat - alpha * x_temp * sinat;
      //
      beta_new = atan2(x(0),x(1));
      //
      nu = nu + epsilon2 * gradient_beta(cbeta_prior, beta_new, tau_yes, tau_no,kappa ,nc,ymat_row);
      nu = nu = (dp - x * x.t()) * nu;

    }
    h_new = likeli_beta(omega,beta_new,tau_yes,tau_no,kappa,nc,ymat_row) - as_scalar(0.5 * nu.t() * nu);
    accept = exp(h_new - h);
    if(accept < as_scalar(arma::randu(1))){
      beta_new = beta_prev;
      e_s = e_s - 1;
      accept_chain(count) = 0;
    }
    beta_out(count) = beta_new;
  }   
  return List::create(beta_out,e_s,accept_chain);
}
// [[Rcpp::export]]
List update_tau_yes(int t, arma::vec nc_par,arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec, int nr, 
                    arma::vec beta,arma::vec tau_yes,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);                  
  int e_s = end - start;
  arma::vec yes_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);
  
  arma::vec nu = arma::zeros(p);
  arma::vec x_temp = arma::zeros(p);
  arma::vec x = arma::zeros(p);
  double alpha, ae, cosat, sinat, h, h_new, accept,tau_new,tau_prev;
  double count = 0;
  
for(int j = start; j< end; j++,count++ ){ 
  double kappa_j = kappa(j);
  double epsilon = epsilon_vec(j);
  double epsilon2 = epsilon2_vec(j);
  int leap = leap_vec(j);
  //reuse
  arma::vec temp_no = pow(acos(cos(tau_no(j)-beta)),2);
  //
  arma::vec ymat_col = ymat.col(j);
  tau_prev = tau_new = tau_yes(j);
  x = s_to_e(tau_prev,x);
  nu = arma::randn(p);
  nu = (dp - x * x.t()) * nu;
  
  h = likeli_tau_yes(tau_prev,beta,kappa_j,nr,ymat_col,temp_no) - as_scalar(0.5 * nu.t() * nu);
  for(int i = 0; i < leap; i++) {
    nu = nu + epsilon2 * gradient_tau_yes(tau_new, beta,kappa_j,nr,ymat_col,temp_no);
    nu = (dp - x * x.t()) * nu;
    
    alpha = norm(nu);
    ae = alpha * epsilon;
    cosat = cos(ae);
    sinat = sin(ae);
    x_temp = x;
    //geodesic flow//
    x = x * cosat + nu/alpha * sinat;
    nu = nu * cosat - alpha * x_temp * sinat;
    //
    tau_new = atan2(x(0),x(1));
    //
    nu = nu + epsilon2 * gradient_tau_yes(tau_new,beta,kappa_j,nr,ymat_col,temp_no);
    nu = nu = (dp - x * x.t()) * nu;
    
  }
  h_new = likeli_tau_yes(tau_new,beta,kappa_j,nr,ymat_col,temp_no) - as_scalar(0.5 * nu.t() * nu);
  accept = exp(h_new - h);
  if(accept < as_scalar(arma::randu(1))){
    e_s -= 1;
    tau_new = tau_prev;
    accept_chain(count) = 0;
  }
  yes_out(count) = tau_new;
}
  return List::create(yes_out,e_s,accept_chain);
}
// [[Rcpp::export]]
List update_tau_no(int t,arma::vec nc_par, arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec, int nr, 
                      arma::vec beta,arma::vec yes_out,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);    
  int e_s = end - start;
  arma::vec no_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);
    
  arma::vec nu = arma::zeros(p);
  arma::vec x_temp = arma::zeros(p);
  arma::vec x = arma::zeros(p);
  double alpha, ae, cosat, sinat, h, h_new, accept,tau_new,tau_prev;
  double count = 0;
  
for(int j = start; j< end; j++,count++ ){   
  double kappa_j = kappa(j);
  double epsilon = epsilon_vec(j);
  double epsilon2 = epsilon2_vec(j);
  int leap = leap_vec(j);
  //reuse
  arma::vec temp_yes = pow(acos(cos(yes_out(j)-beta)),2);
  //
  arma::vec ymat_col = ymat.col(j);
  tau_prev = tau_new = tau_no(j);
  x = s_to_e(tau_prev,x);
  nu = arma::randn(p);
  nu = (dp - x * x.t()) * nu;
  
  h = likeli_tau_no(tau_prev,beta,kappa_j,nr,ymat_col,temp_yes) - as_scalar(0.5 * nu.t() * nu);
  for(int i = 0; i < leap; i++) {
    nu = nu + epsilon2 * gradient_tau_no(tau_new, beta,kappa_j,nr,ymat_col,temp_yes);
    nu = (dp - x * x.t()) * nu;
    
    alpha = norm(nu);
    ae = alpha * epsilon;
    cosat = cos(ae);
    sinat = sin(ae);
    x_temp = x;
    //geodesic flow//
    x = x * cosat + nu/alpha * sinat;
    nu = nu * cosat - alpha * x_temp * sinat;
    //
    tau_new = atan2(x(0),x(1));
    //
    nu = nu + epsilon2 * gradient_tau_no(tau_new,beta,kappa_j,nr,ymat_col,temp_yes);
    nu = nu = (dp - x * x.t()) * nu;
    
  }
  h_new = likeli_tau_no(tau_new,beta,kappa_j,nr,ymat_col,temp_yes) - as_scalar(0.5 * nu.t() * nu);
  accept = exp(h_new - h);
  if(accept < as_scalar(arma::randu(1))){
    tau_new = tau_prev;
    e_s -= 1;
    accept_chain(count) = 0;
  }
  no_out(count) = tau_new;
}
  return List::create(no_out,e_s,accept_chain);
}

// [[Rcpp::export]]
List update_omega(double omega, arma::vec beta, int nr, double a, double b,double omega_sd){
  double eta = log(omega);
  double eta_new = omega_sd * as_scalar(arma::randn(1)) + eta;
  double omega_new = exp(eta_new);
  double accept;
  int accept_out = 1;
  
  accept = exp(likeli_omega(omega_new,beta,nr,a,b) + eta_new - likeli_omega(omega,beta,nr,a,b) - eta);
  if(accept < as_scalar(arma::randu(1))){
    omega_new = omega;
    accept_out = 0;
  }
  return List::create(omega_new,accept_out);
}
// [[Rcpp::export]]
List update_kappa(int t,arma::vec nc_par, int nr,arma::vec beta, arma::vec yes_out, arma::vec no_out,
                      arma::vec kappa,arma::mat ymat, double kappa_a, double ccc, arma::vec t_sig_vec){
                      
  int start = nc_par(t);
  int end = nc_par(t+1);    
  int e_s = end - start;
  arma::vec kappa_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);
  
  double count = 0;
  
for(int j = start; j< end; j++,count++ ){   
  arma::vec asset = sb_c_vec(pow(acos(cos(no_out(j)-beta)),2) - pow(acos(cos(yes_out(j)-beta)),2));
  double kappa_last = kappa(j);
  double t_sig = t_sig_vec(j);
  double tt = log(kappa_last);
  double tt_new = t_sig * as_scalar(arma::randn(1)) + tt;
  double kappa_new = exp(tt_new);
  double accept;
  arma::vec ymat_col = ymat.col(j);

  accept = exp(likeli_kappa(nr,asset,kappa_new,kappa_a,ccc,ymat_col) + tt_new -
    likeli_kappa(nr,asset,kappa_last,kappa_a,ccc,ymat_col) - tt);

  if(accept < as_scalar(arma::randu(1))){
    kappa_new = kappa_last;
    e_s -= 1;
    accept_chain(count) = 0;
  }
  kappa_out(count) = kappa_new;
}
  return List::create(kappa_out,e_s,accept_chain);
}

// [[Rcpp::export]]
List waic_cpp(int t, arma::vec nc_par,int nr,arma::vec beta, arma::vec tau_yes, arma::vec tau_no, arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);
  arma::mat amat = arma::zeros(nr,end-start);
  arma::mat amat_exp = amat;
  arma::mat amat_2 = amat;
  double count = 0;
  double temp = 0;
  for(int j = start; j< end; j++,count++){
    double kappa_j = kappa(j);
    arma::vec ymat_col = ymat.col(j);
    arma::vec asset = sb_c_vec(square(acos(cos(tau_no(j)-beta))) - square(acos(cos(tau_yes(j)-beta))));
    
    for(int i = 0 ;i < nr; i++){
      if(ymat_col(i) == 1){
        amat(i,count) = temp = betaInc_log(kappa_j,kappa_j,asset(i));
      }else{
        amat(i,count) = temp = betaInc_log_lower(kappa_j,kappa_j,asset(i));
      }
        amat_exp(i,count) = exp(temp);
        amat_2(i,count) = pow(temp,2);
    }  
  }
    return List::create(amat,amat_exp,amat_2);
}

// [[Rcpp::export]]
arma::mat predict(int t, arma::vec nc_par,int nr,arma::vec beta, arma::vec tau_yes, arma::vec tau_no, arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);
  arma::mat amat = arma::zeros(nr,end-start);
  arma::mat amat_exp = amat;
  arma::mat amat_2 = amat;
  double count = 0;
  double temp = 0;
  for(int j = start; j< end; j++,count++){
    double kappa_j = kappa(j);
    arma::vec ymat_col = ymat.col(j);
    arma::vec asset = sb_c_vec(square(acos(cos(tau_no(j)-beta))) - square(acos(cos(tau_yes(j)-beta))));
    
    for(int i = 0 ;i < nr; i++){
      temp = betaInc(kappa_j,kappa_j,asset(i));
      if(temp >= 0.5){
        amat(i,count) = 1;
      }  
    }
  }
    return amat;
}

'