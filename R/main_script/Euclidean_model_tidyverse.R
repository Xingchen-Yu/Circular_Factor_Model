
LFM = function(out, n_pos=1000,burnin=500,thin = 1, hyperparams=)

data = 2
iter = 30000
burnin = 10000
House = F
WAIC_group = T
n_pos = iter - burnin

  out = ymat_spit(hn=hn)
  ymat = out[[1]]
  pol = out[[2]]

no_na<-which(is.na(ymat)==F)
len_no_na = length(no_na)
nr<-nrow(ymat)
nc<-ncol(ymat)
a_mu<-a_alpha<-rep(0,2)
#a_beta<-rep(0,nr)
sig_mu_alpha<-diag(2)
sig_mu_alpha_inv<-solve(sig_mu_alpha)
sig_beta<-diag(nr)
sig_beta_inv<-solve(sig_beta)

beta_master<-matrix(0,nr,n_pos)
#mu_master<-alpha_master<-matrix(0,iter,nc)
zmat<-matrix(0,nr,nc)
yes<-which(ymat==1)
yes_l<-length(yes)
no<-which(ymat==0)
no_l<-length(no)
na<-which(is.na(ymat)==T)
na_l<-length(na)
# na_master<-matrix(0,iter,na_l)
yes_mean<-rep(0,yes_l)
no_mean<-rep(0,no_l)
na_mean<-rep(0,na_l)

beta<-rep(0,nr)
if(data==1){
  beta[grep('\\(D',pol)] = -1
  beta[grep('\\(R',pol)] = 1
}
if(WAIC_group==T){
  pos_pred_group = matrix(0,nr,n_pos)
}
part2<-sig_mu_alpha_inv%*%a_mu
##
mu_alpha_mat<-matrix(0,nc,2)
pos_pred<-pos_pred2<-pos_pred3<-matrix(0,nr,nc)
p_index = 1
y_hat = matrix(0,nr,nc)
likeli_chain = rep(0,n_pos)
for(i in 1:iter){
  if(i %in% seq(1,iter,100)){
    cat("\rProgress: ",i,"/",iter)
  }
  ###latent###
  z_yes<-rtruncnorm(1,a=0,b=Inf,mean=yes_mean,sd=1)
  z_no<-rtruncnorm(1,a=-Inf,b=0,no_mean,sd=1)
  z_na<-rnorm(na_l,na_mean,sd=1)
  zmat[yes]<-z_yes
  zmat[no]<-z_no
  zmat[na]<-z_na
  # na_master[i,]<-ifelse(z_na>0,1,0)
  
  ###alpha and mu###
  alpha_beta<-cbind(1,beta)
  po_sig_mual<-chol2inv(chol((crossprod(alpha_beta)+sig_mu_alpha_inv)))
  for(j in 1:nc){
    po_mean_mual<-po_sig_mual%*%(crossprod(alpha_beta,zmat[,j])+part2)
    up_mu_alpha<-rmvn(1,po_mean_mual,po_sig_mual)
    mu_alpha_mat[j,]<-up_mu_alpha
  }
  mu<-mu_alpha_mat[,1]
  #mu_master[i,]<-mu
  alpha<-mu_alpha_mat[,2]
  #alpha_master[i,]<-alpha
  ###beta###
  tt<-t(zmat)-mu
  var_beta<-c(chol2inv(chol((crossprod(alpha)+1))))
  po_mean_beta<-var_beta%*%crossprod(alpha,tt)
  beta<-c(rmvn(1,po_mean_beta,diag(var_beta,nr)))
  
  mean_mat<-t(tcrossprod(alpha,beta) + mu)
  
  yes_mean<-mean_mat[yes]
  no_mean<-mean_mat[no]
  na_mean<-mean_mat[na]
  if(i %in% seq(1,iter,25)){
    y_p = pnorm(mean_mat)
    index_yes = which(y_p>0.5)
    y_hat[index_yes] = 1
    y_hat[-index_yes] = 0
    ymat[na] = y_hat[na]
    # print(paste('trainning accuracy is',length(which(y_hat[no_na]==ymat[no_na]))/len_no_na))
  }
  
  if(i>burnin){
    beta_master[,p_index]<-beta
    p_s1<-pnorm(t(tcrossprod(alpha,beta) + mu))
    kobe<-dbinom(ymat,1,p_s1,log = T)
    likeli_chain[p_index] = sum(kobe[no_na])
    
    if(WAIC_group==T){
      kobe2 = apply(kobe,1,sum)
      pos_pred_group[,p_index] = kobe2
    }
    
    pos_pred<-pos_pred + exp(kobe)
    pos_pred2<-pos_pred2 + kobe
    pos_pred3<-pos_pred3 + kobe^2
    p_index = p_index + 1
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
# pos_pred_master<-pos_pred[no_na]/n_pos
# lpd<-sum(log(pos_pred_master))
# va<-sum((pos_pred3[no_na]/n_pos-(pos_pred2[no_na]/n_pos)^2))*n_pos/(n_pos-1)
# waic_eucliwdean<-lpd-va
# waic_master[hh]=waic_euclidean
#waic_master[hh,3] = waic_euclidean
#waic_master[hh,1:2]=c(nr,nc)
# if(House==T){
#   save(file=paste0("H",hn,"_betaChain_eu.Rdata"),beta_master)
# }else{
#   save(file=paste0("S",hn,"_betaChain_eu.Rdata"),beta_master)
# }
# save.image(file=paste0(h_s,hn,"_workspace.Rdata"))
# save.image(file=paste0('../paper1_house_eu/',"H",hn,"_workspace.Rdata"))
# rm(list=setdiff(ls(), c('hh')))
# print(waic_master)
# }
# save(file='../paper1_house_eu/waic_master.Rdata',waic_master)

# save.image(file=paste0(h_s,hn,'_',1,"_nd_eu_workspace.Rdata"))
# plot(likeli_chain)
# # }
# save(file = "Senate_WAIC.Rdata",waic_master)
#
# waic_master_eu = waic_master
# save(file = "Senate_WAIC",waic_master)
#
# senate_all = cbind(waic_master,waic_master_eu[,3])
# colnames(senate_all) = c("I","J","Spherical","Euclidean")
#
# save(file="Senate_WAIC_master.Rdata",senate_all)
# write.csv(file="Senate_WAIC_master.csv",senate_all)
# save(file='h112_beta_master_eu.Rdata',beta_master)
