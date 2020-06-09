ymat_spit<-function(hn){
  vote2<-readKH2(file=paste0("../data/H",hn,"_votes.ord"))
  if(hn==116){
    vote<-as.matrix(vote2$votes)[,1:700]
    ####Amash switches party at 430, we combine the record for Amash together####
    vote[202,430:700] = vote[203,430:700]
    vote = vote[-203,]
  }else{
    vote<-as.matrix(vote2$votes)
  }
  pol<-rownames(vote)
  ####0,7,8,9 represents NA###
  ind3<-which(vote==0 | vote==7 | vote==8| vote==9)
  vote[ind3]<-NA
  ####1,2,3 represents Yea####
  ind1<-which(vote==1 | vote==2 | vote==3)
  vote[ind1]<-1
  ####4,5,6 represents Nah####
  ind2<-which(vote==4 | vote==5| vote==6)
  vote[ind2]<-0
  #####legislators who miss more than 40% of the vote is excluded####
  abs_percent<-as.numeric(apply(vote,1,function(xx) length(which(is.na(xx)==T))))/ncol(vote)
  abs_ind<-which(abs_percent>=0.4)
  if(length(abs_ind)>0){
    pol<-pol[-abs_ind]
    ymat<-vote[-abs_ind,]
  }else{
    ymat<-vote
  }
  colnames(ymat)<-NULL
  return(list(ymat,pol))
}
waic_compute = function(n_pos,pos_pred,pos_pred2,pos_pred3,no_na){
  pos_pred_master = pos_pred[no_na]/n_pos
  lpd = sum(log(pos_pred_master))
  va = sum((pos_pred3[no_na]/n_pos-(pos_pred2[no_na]/n_pos)^2))*n_pos/(n_pos-1)
  waic_spherical = lpd-va
  print(waic_spherical)
}
update_tsig = function(upper,lower,t_sig,check,nnn){
  l_s = which(check<=lower)
  u_s = which(check>=upper)
  
  t_sig[l_s] = t_sig[l_s] * 0.95
  t_sig[u_s] = t_sig[u_s] * 1.05
  return(list(t_sig,length(c(l_s,u_s))/nnn))
}
