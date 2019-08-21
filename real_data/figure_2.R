allest<-NULL
allub<-NULL
alllb<-NULL

for (i in 1:2){
  
  yr1<-readLines(paste0('results_yr',i,'.txt'))
  
  ## clean the data ##
  # get rid of header stuff
  yr1<-yr1[-(1:(which(yr1=="All Matching Results")-1))]
  
  # only keep elements with the word "events"
  yr1<-yr1[grep('Events:',yr1)]
  
  temp<-strsplit(x=yr1,split="[, ()]")
  
  yr1num<-as.numeric(unlist(lapply(temp,function(x) x[-c(1,2,4)])))
  
  ## extract point estimates ##
  est<-yr1num[seq(1,length(yr1num),3)]
  
  ## extract LBs ##
  lb<-yr1num[seq(1,length(yr1num),3)+1]
  
  ## extract UBs ##
  ub<-yr1num[seq(1,length(yr1num),3)+2]
  
  allest<-c(allest,est)
  allub<-c(allub,ub)
  alllb<-c(alllb,lb)

}

disease<-factor(rep(rep(c('Mortality','Dementia','CVD'),3),2),levels=c('Mortality','Dementia','CVD'))
Method<-factor(rep(rep(c('Match','BART','PR'),each=3),2),levels=c('Match','BART','PR'))
yr<-factor(rep(c('PA-2000','PA-2001'),each=9),levels=c('PA-2000','PA-2001'))

all_results<-data.frame(allest,alllb,allub,disease,Method,yr)

## make a grid of plots from these results ##
plot_results<-ggplot(all_results,aes(x=Method,y=allest,colour=Method))+
  geom_errorbar(aes(ymin=alllb,ymax=allub))+
  geom_point()+
  geom_hline(yintercept=0)+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text=element_text(size=18))+
  facet_grid(rows=vars(disease),cols=vars(yr))

pdf('figure_2.pdf')
print(plot_results)
dev.off()
