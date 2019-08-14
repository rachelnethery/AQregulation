set.seed(2)

nboot<-50
calp<-0.1

source('functions/exmatch_trt_v2.R')
source('functions/mdmatch_conf_v2.R')
source('functions/match_rn_app_v2.R')
source('functions/match_rn_boot5.R')

## load the data for analysis ##
load('analysis_data.RData')

## for both years of data ##
for (yr in 1:2){
  
  ## pick off that year's dataset ##
  adat<-mdat[[yr]]
  
  N<-nrow(adat)

  ## create separate matricies of data for the observed exposures, counterfactual exposures, and confounders ##
  Eobs<-as.matrix(adat[,c('pmWith','ozWith')])
  Ecf<-as.matrix(adat[,c('pmNo','ozNo')])
  X<-as.matrix(adat[,grep('X',names(adat))])
  
  ## create an external text file to save results ##
  sink(paste0('results_yr',yr,'.txt'))
  
  cat('=====================================\n')
  cat("Data Features\n")
  cat('=====================================\n')
  cat(paste0('Units in original dataset=',N,'\n'))
  
  ##############
  ## MATCHING ##
  ##############
  
  ## omega values, tolerances for exact matching on E ##
  ecut<-apply(Ecf,2,sd)*calp
  xx<-matrix(ecut,nrow=nrow(Ecf),ncol=ncol(Ecf),byrow=T)
  
  ## caliper for mahalanobis distance matching on confounders ##
  quan<-NULL
  S_inv<-solve(cov(X))
  for (j in 1:N){
    temp<-matrix(rep(X[j,],N-1),nrow=N-1,byrow=T)
    quan<-c(quan,quantile(sqrt(rowSums((temp-X[-j,])%*%S_inv*(temp-X[-j,]))),.1))
  }
  
  ## carry out matching ##
  ## the matching function outputs a dataframe with columns
  ## trtind: the indices of the units for which we found matches for their counterfactual pollutants and their confounders
  ## ctlind: the indices of the matched units for the unit given in the trtind column in the same row
  ## note that units for which multiple matches are found will be repeated in the trtind column,
  ## with each row having a different ctlind value corresponding to a different matched unit
  orig.matchout.1<-match_rn_app_v2(trtobs=Eobs,trtcf=Ecf,confounders=X,trtdiff=xx,mdqtl=mean(quan))
  orig.matchind.1<-as.data.frame(orig.matchout.1)
  names(orig.matchind.1)<-c('trtind','ctlind')
  keep.1<-unique(orig.matchind.1$trtind)
  if (yr==1){
    save(orig.matchind.1,keep.1,file='match_output.RData')
  }
  
  cat(paste0('Units in trimmed dataset=',length(keep.1),'\n'))
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("Matching Point Estimates\n")
  cat('=====================================\n')
  
  pname<-c('Mortality Events','Neuro Events','CVD Events')
  
  ## bias correction and point estimates ##
  allest<-NULL
  
  ## for each of the 3 health outcomes ##
  for (i in 2:4){
    ## pick off the health outcome ##
    Yobs<-as.matrix(adat[,i])
    
    ## fit regression for bias correction ##
    regdat<-data.frame(Yobs,Eobs,X)
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
    fit.1<-glm(y~pm25+o3+conf1+conf2+conf3+conf4+conf5, family="poisson", data=regdat)
    ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
    predXi<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$trtind,])
    predXj<-data.frame(Ecf[orig.matchind.1$trtind,],X[orig.matchind.1$ctlind,])
    names(predXi)<-c('pm25','o3',paste0('conf',1:5))
    names(predXj)<-c('pm25','o3',paste0('conf',1:5))
    ## take difference to be used for bias correction ##
    diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
    ## impute counterfactuals ##
    EYcf<-tapply(Yobs[orig.matchind.1$ctlind]+diffmeans,orig.matchind.1$trtind,mean)
    ## difference in imputed counterfactuals and observed number of deaths in each area ## 
    orig.est.1<-sum(EYcf-Yobs[as.numeric(names(EYcf))])
    
    ## write results to the text file ##
    cat(paste0('Number of ',pname[i-1]," Prevented: ",round(orig.est.1,4),'\n'))
    allest<-c(allest,orig.est.1)
  }
  
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("All Matching Results\n")
  cat('=====================================\n')

  ######################
  ## do bootstrapping ##
  ######################
  
  Yobs<-as.matrix(adat[,2])
  
  ## initiate a list to save the bootstrapped data ##
  allboot<-list()
  for (j in 1:nboot){
    
    ## bootstrap from the trimmed sample ##
    bootind1<-sample(x=keep.1,size=length(keep.1),replace=T)
    Yobs.TSboot<-Yobs[bootind1]
    Eobs.TSboot<-Eobs[bootind1,]
    Ecf.TSboot<-Ecf[bootind1,]
    X.TSboot<-X[bootind1,]
    
    ## bootstrap from the untrimmed sample ##
    bootind2<-sample(x=1:N,size=N,replace=T)
    Yobs.USboot<-Yobs[bootind2]
    Eobs.USboot<-Eobs[bootind2,]
    Ecf.USboot<-Ecf[bootind2,]
    X.USboot<-X[bootind2,]
    
    ## save the bootstrapped data in the list ##
    allboot[[j]]<-list(Yobs.TSboot,Eobs.TSboot,Ecf.TSboot,X.TSboot,Yobs.USboot,Eobs.USboot,Ecf.USboot,X.USboot,bootind1,bootind2)
  }
  
  ## run the matching on the bootstrapped data as described in Section 3 of the paper ##
  matchout.1<-lapply(allboot,match_rn_boot5,mdqtl=mean(quan),calpE=xx[1:length(keep.1),])
  
  ## resulting indices are in terms of the bootstrapped datasets, correct them so that they point to the index in the original data ##
  for (j in 1:nboot){
    matchout.1[[j]]$trtind.orig<-allboot[[j]][[9]][matchout.1[[j]]$trtind]
    matchout.1[[j]]$ctlind.orig<-allboot[[j]][[10]][matchout.1[[j]]$ctlind]
  }
  

  ##########################################################
  ## compute matching estimates for each bootstrap sample ##
  ##########################################################
  
  ## create an object to save results ##
  boot_results<-NULL
  
  ## for each health outcome ##
  for (i in 2:4){
    Yobs<-as.matrix(adat[,i])
    
    est.method1.boot<-NULL
    for (j in 1:nboot){
      matchind.1<-matchout.1[[j]]
      ## fit regression with bootstrapped data ##
      regdat<-data.frame(Yobs[allboot[[j]][[10]]],Eobs[allboot[[j]][[10]],],X[allboot[[j]][[10]],])
      names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
      fit.1<-glm(y~pm25+o3+conf1+conf2+conf3+conf4+conf5, family="poisson", data=regdat)
      ## for each CF level in matched data, use regression function to predict outcome with true confounders and with confounders of matches ##
      predXi<-data.frame(Ecf[matchind.1$trtind.orig,],X[matchind.1$trtind.orig,])
      predXj<-data.frame(Ecf[matchind.1$trtind.orig,],X[matchind.1$ctlind.orig,])
      names(predXi)<-c('pm25','o3',paste0('conf',1:5))
      names(predXj)<-c('pm25','o3',paste0('conf',1:5))
      ## take difference to be used for bias correction ##
      diffmeans<-predict(fit.1, predXi, type="response")-predict(fit.1, predXj, type="response")
      ## impute counterfactuals ##
      EYcf<-tapply(Yobs[matchind.1$ctlind.orig]+diffmeans,matchind.1$trtind,mean)
      
      unqtrt<-NULL
      for (k in 1:max(matchind.1$trtind)){
        unqtrt<-c(unqtrt,matchind.1$trtind.orig[which(matchind.1$trtind==k)][1])
      }
      
      ## difference in imputed counterfactuals and observed number of deaths in each area ## 
      ice_hat.1<-EYcf-Yobs[unqtrt]
      ## sum over all areas to get total number of deaths prevented ##
      est.method1.boot<-c(est.method1.boot,sum(ice_hat.1))
    }
    boot_results<-rbind(boot_results,est.method1.boot)
  }

  ## write results to the text file ##
  for (i in 1:3){
    cat(paste0(pname[i],": ",round(allest[i],4)," (",round(quantile(boot_results[i,],.025),4),',',round(quantile(boot_results[i,],.975),4),')\n'))
  }
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("BART Results\n")
  cat('=====================================\n')
  
  ##########
  ## BART ##
  ##########
  
  ## for each health outcome ##
  for (i in 2:4){
    ## pick off the health outcome ##
    Yobs<-as.matrix(adat[,i])
    
    ## fit the BART model to the observed data ##
    ## use BART to predict counterfactuals ##
    bartfit<-bart(x.train=cbind(Eobs,X),y.train=as.numeric(Yobs),x.test=cbind(Ecf,X)[keep.1,],nskip=1000,ndpost=800,verbose=F)
    
    ## use BART to estimate the TEA (note that we are only using the units retained after trimming for estimation) ##
    bartest_tr<-sum(bartfit$yhat.test.mean-Yobs[keep.1])
    
    ## compute the credible interval ##
    tau_pd_tr<-NULL
    for (j in 1:800){
      tau_pd_tr<-c(tau_pd_tr,sum(rnorm(n=length(keep.1),mean=bartfit$yhat.test[j,],sd=bartfit$sigma[j])-Yobs[keep.1]))
    }
    bart_ci<-quantile(tau_pd_tr,c(.025,.975))
    
    ## write results to the text file ##
    cat(paste0(pname[i-1],': ',round(bartest_tr,4)," (",round(bart_ci[1],4),',',round(bart_ci[2],4),')\n'))
    
  }
  
  cat('\n')
  cat('\n')
  cat('\n')
  
  cat('=====================================\n')
  cat("Linear Model Results\n")
  cat('=====================================\n')
  
  ##################
  ## LINEAR MODEL ##
  ##################
  
  ## for each health outcome ##
  for (i in 2:4){
    ## pick off the health outcome ##
    Yobs<-as.matrix(adat[,i])
    
    ## fit linear model including all variables from unmatched data ##
    regdat<-data.frame(Yobs,Eobs,X)
    names(regdat)<-c('y','pm25','o3',paste0('conf',1:5))
    fit.3<-glm(y~pm25+o3+conf1+conf2+conf3+conf4+conf5, family="poisson", data=regdat)
    ## predict at counterfactual pollution levels ##
    predXi<-data.frame(Ecf[keep.1,],X[keep.1,])
    names(predXi)<-c('pm25','o3',paste0('conf',1:5))
    EYcf<-predict(fit.3, predXi, type="response")
    ice_hat.3<-EYcf-Yobs[keep.1]
    lmest<-sum(ice_hat.3)
    ## compute var-cov matrix for the predictions so that we can get the variance of the their sum ##
    vdat<-data.frame(Yobs[keep.1],Ecf[keep.1,],X[keep.1,])
    names(vdat)<-c('y','pm25','o3',paste0('conf',1:5))
    foo<-model.matrix(y~pm25+o3+conf1+conf2+conf3+conf4+conf5,data=vdat)
    pred_vcov<-(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))%*%vcov(fit.3)%*%
      t(foo*matrix(rep(EYcf,ncol(foo)),nrow=length(keep.1),ncol=ncol(foo),byrow=F))
    
    ## compute the confidence limits ##
    lmlb<-sum(ice_hat.3)-1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))
    lmub<-sum(ice_hat.3)+1.96*sqrt(matrix(1,ncol=length(keep.1),nrow=1)%*%pred_vcov%*%matrix(1,ncol=1,nrow=length(keep.1)))
    
    ## write results to the text file ##
    cat(paste0(pname[i-1],': ',round(lmest,4)," (",round(lmlb,4),',',round(lmub,4),')\n'))
  }
  
  ## stop writing to the text file ##
  sink()
  
}