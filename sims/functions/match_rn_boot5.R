## write a little function to (1) select treatment matches and (2) among treatment matches, select confounder matches ##
match_rn_boot5<-function(bootlist,mdqtl,calpE){
  TSboot<-bootlist[1:4]
  USboot<-bootlist[5:8]
  
  ## components of TSboot and USboot are bootstrapped:
  #1. Yobs
  #2. Eobs
  #3. Ecf
  #4. X
  
  #trtdiff<-TSboot[[3]]*calpE
  trtdiff<-calpE
  
  ## get inverse covariance matrix of counfounders to compute mahalanobis distance ##
  S_inv<-solve(cov(USboot[[4]]))
  
  all_trt_match_ind<-mapply(exmatch_trt_v2,split(TSboot[[3]],1:nrow(TSboot[[3]])),split(trtdiff,1:nrow(trtdiff)),MoreArgs=list(trtobs=USboot[[2]]))
  
  m.ind<-mapply(mdmatch_conf_v2,all_trt_match_ind,split(TSboot[[4]],1:nrow(TSboot[[4]])),MoreArgs=list(confounders=USboot[[4]],S_inv=S_inv,mdqtl=mdqtl))
  
  #message(length(which(is.na(m.ind)==T)))
  
  ## for units with no matches, require them to have one match ##
  ## to find this match, get 10 closest units on treatment and choose the one that is closest in terms of confounders ##
  S_inv_E<-solve(cov(USboot[[2]]))
  S_inv_X<-solve(cov(USboot[[4]]))
  for (i in which(is.na(m.ind)==T)){
    temp<-matrix(as.numeric(TSboot[[3]][i,]),nrow=nrow(USboot[[2]]),ncol=ncol(USboot[[2]]),byrow=T)
    MD_E<-sqrt(rowSums(as.matrix(temp-USboot[[2]])%*%S_inv_E*(as.matrix(temp-USboot[[2]]))))
    ## 10 smallest MD's ##
    small_mde<-order(MD_E)[1:20]
    
    temp<-matrix(as.numeric(TSboot[[4]][i,]),nrow=length(small_mde),ncol=ncol(USboot[[4]]),byrow=T)
    MD_i<-sqrt(rowSums(as.matrix(temp-USboot[[4]][small_mde,])%*%S_inv_X*(as.matrix(temp-USboot[[4]][small_mde,]))))
    m.ind[[i]]<-small_mde[which.min(MD_i)]
  }

  # construct a dataset with treatment indices in one column and their matched control indices in the other ##
  matchdat<-NULL
  for (i in 1:length(m.ind)){
    matchdat<-rbind(matchdat,cbind(rep(i,length(m.ind[[i]])),m.ind[[i]]))
  }

# matchdat<-NULL
# for (i in 1:length(m.ind)){
#   if (sum(is.na(m.ind[[i]]))==length(m.ind[[i]])){}
#   else{
#     matchdat<-rbind(matchdat,cbind(rep(i,length(m.ind[[i]])),m.ind[[i]]))
#   }
# }


  matchdat<-as.data.frame(matchdat)
  names(matchdat)<-c('trtind','ctlind')

  return(matchdat)
  
}
