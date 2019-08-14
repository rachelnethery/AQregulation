## write a little function to (1) select treatment matches and (2) among treatment matches, select confounder matches ##
match_rn_app_v2<-function(trtobs,trtcf,confounders,trtdiff,mdqtl){
  
  ## get inverse covariance matrix of counfounders to compute mahalanobis distance ##
  S_inv<-solve(cov(confounders))
  sds<-apply(trtobs,2,sd)  
  
  all_trt_match_ind<-mapply(exmatch_trt_v2,split(trtcf,1:nrow(trtcf)),split(trtdiff,1:nrow(trtdiff)),MoreArgs=list(trtobs=trtobs))
  
  m.ind<-mapply(mdmatch_conf_v2,all_trt_match_ind,split(confounders,1:nrow(confounders)),MoreArgs=list(confounders=confounders,S_inv=S_inv,mdqtl=mdqtl))


  ## construct a dataset with treatment indices in one column and their matched control indices in the other ##
  matchdat<-NULL
  for (i in 1:length(m.ind)){
    if (sum(is.na(m.ind[[i]]))==length(m.ind[[i]])){}
    else{
      matchdat<-rbind(matchdat,cbind(rep(i,length(m.ind[[i]])),m.ind[[i]]))
    }
  }
  
  colnames(matchdat)<-c('trtind','ctlind')
  
  return(matchdat)
  
}
