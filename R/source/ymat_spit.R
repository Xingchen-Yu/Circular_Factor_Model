ymat_spit<-function(hn,house){
  if(house==T){
    vote2<-readKH2(file=paste0("./data/H",hn,"_votes.ord"))
  }else{
    vote2<-readKH2(file=paste0("./data/S",hn,"_votes.ord"))
  }
  # vote2<-readKH2(file=paste0("./data/H",hn,"_votes.ord"))
  if(hn==116 & house==T){
    # vote<-as.matrix(vote2$votes)[,1:700]
    vote<-as.matrix(vote2$votes)
    ####Amash switches party at 430, we combine the record for Amash together####
    vote[202,430:ncol(vote)] = vote[203,430:ncol(vote)]
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
source(file="./source/utility_func.R")
