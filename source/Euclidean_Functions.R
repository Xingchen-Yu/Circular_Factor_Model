'%!in%' <- function(x,y)!('%in%'(x,y))
# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)
readKH2<-function (file, dtl = NULL, yea = c(1, 2, 3), nay = c(4, 5, 
                                                               6), missing = c(7, 8, 9), notInLegis = 0, desc = NULL, debug = FALSE) 
{
  cat("Attempting to read file in Keith Poole/Howard Rosenthal (KH) format.\n")
  warnLevel <- options()$warn
  options(warn = -1)
  data <- try(readLines(con = file), silent = TRUE)
  if (inherits(data, "try-error")) {
    cat(paste("Could not read", file, "\n"))
    return(invisible(NULL))
  }
  options(warn = warnLevel)
  cat("Attempting to create roll call object\n")
  voteData <- substring(data, 37)
  n <- length(voteData)
  m <- nchar(voteData)[1]
  rollCallMatrix <- matrix(NA, n, m)
  for (i in 1:n) {
    asdf<-gsub(" ","",voteData[i])
    rollCallMatrix[i, ] <- as.numeric(unlist(strsplit(asdf, 
                                                      split = character(0))))
  }
  rm(voteData)
  if (!is.null(desc)) 
    cat(paste(desc, "\n"))
  cat(paste(n, "legislators and", m, "roll calls\n"))
  cat("Frequency counts for vote types:\n")
  tab <- table(rollCallMatrix, exclude = NULL)
  print(tab)
  icpsrLegis <- as.numeric(substring(data, 4, 8))
  party <- as.numeric(substring(data, 21, 23))
  partyfunc <- function(x) {
    party <- partycodes$party[match(x, partycodes$code)]
    party[party == "Democrat"] <- "D"
    party[party == "Republican"] <- "R"
    party[party == "Independent"] <- "Indep"
    party
  }
  partyName <- partyfunc(party)
  statename <- function(x) {
    state.info$state[match(x, state.info$icpsr)]
  }
  state <- as.numeric(substring(data, 9, 10))
  KHstateName <- substring(data, 13, 20)
  stateName <- statename(state)
  stateAbb <- datasets::state.abb[match(stateName, datasets::state.name)]
  stateAbb[grep(KHstateName, pattern = "^USA")] <- "USA"
  cd <- as.numeric(substring(data, 11, 12))
  cdChar <- as.character(cd)
  cdChar[cd == 0] <- ""
  lnames <- substring(data, 26, 36)
  for (i in 1:n) {
    #lnames[i] <- strip.trailing.space(lnames[i])
    lnames[i] <- trim.leading(lnames[i])
    lnames[i] <- trim.trailing(lnames[i])
    #lnames[i] <- strip.after.comma(lnames[i])
  }
  legisId <- paste(lnames, " (", partyName, " ", stateAbb, 
                   "-", cdChar, ")", sep = "")
  legisId <- gsub(x = legisId, pattern = "-)", replacement = ")")
  if (any(duplicated(legisId))) {
    dups <- duplicated(legisId)
    legisId[dups] <- paste(legisId[dups], icpsrLegis[dups])
  }
  legis.data <- data.frame(state = stateAbb, icpsrState = state, 
                           cd = cd, icpsrLegis = icpsrLegis, party = partyName, 
                           partyCode = party)
  dimnames(legis.data)[[1]] <- legisId
  vote.data <- NULL
  if (!is.null(dtl)) {
    vote.data <- dtlParser(dtl, debug = debug)
  }
  rc <- rollcall(data = rollCallMatrix, yea = yea, nay = nay, 
                 missing = missing, notInLegis = notInLegis, legis.names = legisId, 
                 legis.data = legis.data, vote.data = vote.data, desc = desc, 
                 source = file)
  rc
}
ymat_spit<-function(hn,house){
  if(house==T){
    vote2<-readKH2(file=paste0("H",hn,"_votes.ord"))
  }else{
    vote2<-readKH2(file=paste0("S",hn,"_votes.ord"))
  }
  vote<-as.matrix(vote2$votes)
  pol<-rownames(vote)
  ind3<-which(vote==0 | vote==7 | vote==8| vote==9)
  vote[ind3]<-NA
  ind1<-which(vote==1 | vote==2 | vote==3)
  vote[ind1]<-1
  ind2<-which(vote==4 | vote==5| vote==6)
  vote[ind2]<-0
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

waic_compute = function(nnn,pos_red,pos_pred2,pos_pred3,no_na){
  pos_pred_master<-pos_pred[no_na]/nnn
  lpd<-sum(log(pos_pred_master))
  va<-sum((pos_pred3[no_na]/nnn-(pos_pred2[no_na]/nnn)^2))*nnn/(nnn-1)
  waic_euclidean_pd<-lpd-va
  return(waic_euclidean_pd)
}
predict_ymat = function(aa,bb,mm){
  y_hat = matrix(0,nr,nc)
  mean_mat<-t(tcrossprod(aa,bb) + mm)
  y_p = pnorm(mean_mat)
  index_yes = which(y_p>0.5)
  y_hat[index_yes] = 1
  return(y_hat)
}