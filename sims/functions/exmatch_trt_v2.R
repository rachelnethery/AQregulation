exmatch_trt_v2<-function(foo0,foo1,trtobs){
  lwr<-foo0-foo1
  upr<-foo0+foo1
  
  #trt_matches<-apply(trtobs,1,function(x) as.numeric(sum(x>lwr & x<upr)==length(x)))
  blah<-matrix(NA,nrow=nrow(trtobs),ncol=length(lwr))
  for (v in 1:length(lwr)){
    blah[,v]<-as.numeric(between(trtobs[,v],lwr[v],upr[v]))
  }
  trt_matches<-apply(blah,1,function(x) sum(x)==length(x))
  
  if (sum(trt_matches)==0){
    all_trt_match_ind<-NA
  } else{
    all_trt_match_ind<-which(trt_matches>0)
  }
  return(all_trt_match_ind)
}