mdmatch_conf_v2<-function(foo2,foo3,confounders,S_inv,mdqtl){
  if (sum(is.na(foo2))==length(foo2)){
    m.ind<-NA
  } else{
    temp<-matrix(foo3,nrow=nrow(confounders),ncol=ncol(confounders),byrow=T)
    
    MD_i<-sqrt(rowSums(as.matrix(temp-confounders)%*%S_inv*(as.matrix(temp-confounders))))
    #MD_cut<-quantile(MD_i,probs=mdqtl)
    
    if (sum(foo2 %in% which(MD_i<=mdqtl))==0){
      m.ind<-NA
    } else {
      m.ind<-foo2[which(foo2 %in% which(MD_i<=mdqtl))]
    }
  }
  return(m.ind)
}