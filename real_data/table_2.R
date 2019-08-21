load('analysis_data.RData')
load('match_output.RData')

adat<-mdat[[1]]
N<-nrow(adat)

adat<-adat[,-1]

meansall<-round(apply(adat,2,mean),2)
sdsall<-round(apply(adat,2,sd),2)

meanstrim<-round(apply(adat[-keep.1,],2,mean),2)
sdstrim<-round(apply(adat[-keep.1,],2,sd),2)

meansuntrim<-round(apply(adat[keep.1,],2,mean),2)
sdsuntrim<-round(apply(adat[keep.1,],2,sd),2)

tab1<-data.frame(paste0(meansall,' (',sdsall,')'),paste0(meanstrim,' (',sdstrim,')'),paste0(meansuntrim,' (',sdsuntrim,')'))
rownames(tab1)<-names(adat)
names(tab1)<-c('Original Sample','Discarded Units','Trimmed Sample')
sink('table_2.txt')
print(tab1)
sink()

