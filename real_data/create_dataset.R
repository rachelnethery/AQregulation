N<-10000

set.seed(2)
  
## simulate 2 factual (observed) pollutants similar to pm2.5 annual average and o3 summer average ##
Eobs<-rmvnorm(n=N,mean=c(10.5,46.9),sigma=diag(2)*c(7,37.8))

## simulate 2 counterfactual pollutants with larger mean and larger spread ##
Ecf<-Eobs+cbind(runif(n=N,min=0,max=7.08),runif(n=N,min=0,max=15))

## simulate confounders ##
betaE<-matrix(c(0.1,0.05,-0.1,-0.05,0.02,.012,-.005,-.0035,.08,-.05),nrow=5,byrow=F)
sdE<-c(0.5,0.85,1,1,1)
for (i in 1:5){
  confi<-rnorm(n=N,mean=Eobs%*%t(betaE[i,,drop=F]),sd=sdE[i])
  if (i==1){
    X<-matrix(confi)
  } else{
    X<-cbind(X,confi)
  }
}

mdat<-list()

for (yr in 1:2){
  
  ## beta values and predictor forms ##
  betaconf<-c(0.001400702, -0.0024289523,  0.0049607745, -0.0037199380,  0.00023785807,  -0.0033575231,  0.0023112471, -0.0034723793)
  
  pform.obs<-cbind(1,X,X[,1]*X[,2],X[,3]^2,exp(X[,4])/(1+exp(X[,4])),Eobs,Eobs[,1]^2,Eobs[,1]*Eobs[,2],Eobs[,1]*Eobs[,2]*X[,5],Eobs[,1]*Eobs[,2]*X[,4])
  pform.cf<-cbind(1,X,X[,1]*X[,2],X[,3]^2,exp(X[,4])/(1+exp(X[,4])),Ecf,Ecf[,1]^2,Ecf[,1]*Ecf[,2],Ecf[,1]*Ecf[,2]*X[,5],Ecf[,1]*Ecf[,2]*X[,4])
  
  for (health in 1:3){
    ## simulate observed counts of mortality (health=1), dementia (health=2), and cvd (health=3), associated with pollutants and confounders ##
    if (health==1){
      beta<-matrix(c(3.75,betaconf,.008,.0005,.0017,.00001,.0005,.0003),ncol=1)
      Yobs<-matrix(rpois(n=N,lambda=exp(pform.obs%*%beta)))
    } else if (health==2){
      beta<-matrix(c(2.5,betaconf,.008,.0005,.0017,.00001,.0005,.0003),ncol=1)
      Yobs<-cbind(Yobs,rpois(n=N,lambda=exp(pform.obs%*%beta)))
    } else{
      beta<-matrix(c(4,betaconf,.008,.0005,.0017,.00001,.0005,.0003),ncol=1)
      Yobs<-cbind(Yobs,rpois(n=N,lambda=exp(pform.obs%*%beta)))
    }
  }
  
  alldat<-data.frame(1:N,Yobs,Eobs,Ecf,X)
  
  names(alldat)<-c('id','mort','dementia','cvd','pmWith','ozWith','pmNo','ozNo','X1','X2','X3','X4','X5')
  
  mdat<-c(mdat,list(alldat))
  
}

save(mdat,file='analysis_data.RData')
