setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/biobank/clinicaldata")
library(Rtsne)  
library("ggplot2")  
data<-read.table("GHCC.RA.CCP.uni.R",head=T,sep="\t",row.names=1)
data=bdata(data)
d <- dist(data) # euclidean distances between the rows
rlt<-Rtsne(d, is_distance=T)
pdf("tsne_RA_clinical.pdf")
tsne = as.data.frame(rlt$Y)  
plot(tsne$Y,pch=16,cex=0.5)
dev.off()
load("tsne.RData")
cols <- rainbow(13)
label=kmeans(tsne,13)
plot(tsne, t='n')
text(tsne, labels=label$cluster, col=data$CCP,cex=0.8)


par(mfrow=c(2,2))

plot(tsne, pch=16,col=varible2color(data$CCP),cex=0.8)
text(tsne, labels=label$cluster, col=varible2color(data$CCP),cex=0.8)

plot(tsne, t='n')
text(tsne, labels=label$cluster, col=varible2color(data$RFIGM),cex=0.8)

plot(tsne, t='n')
text(tsne, labels=label$cluster, col=varible2color(data$ESR),cex=0.8)

plot(tsne, t='n')
text(tsne, labels=label$cluster, col=varible2color(data$CRP),cex=0.8)

plot(tsne, t='n')
text(tsne, labels=label$cluster, col=varible2color(data$RFIGG),cex=0.8)

plot(tsne, t='n')
text(tsne, labels=label$cluster, col=varible2color(data$RFIGA),cex=0.8)



g9=which(label$cluster==9)
c9=which(label$cluster !=9)

rlt<-apply(data,2,function(x) t.test(x[g9],x[c9],na.rm=T))
rbind(lapply(rlt,function(x) x$p.value),lapply(rlt,function(x) x$statistic))



varible2color<-function(x){
  x<-x+1
  x[is.na(x)]<-0
  x
}


bdata<-function(temp){
temp$CCP[temp$CCP<=25]<-1
temp$CCP[temp$CCP>25]<-2
temp$CRP[temp$CRP<=10]<-1
temp$CRP[temp$CRP>10]<-2
temp$ESR[temp$ESR<=30]<-1
temp$ESR[temp$ESR>30]<-2
temp$RFIGM[temp$RFIGM<=30]<-1
temp$RFIGM[temp$RFIGM>30]<-2
temp$RFIGG[temp$RFIGG<=18]<-1
temp$RFIGG[temp$RFIGG>18]<-2
temp$RFIGA[temp$RFIGA<=18]<-1
temp$RFIGA[temp$RFIGA>18]<-2
temp$WBC[temp$WBC<=9.16]<-1
temp$WBC[temp$WBC>9.16]<-2
return(temp)
}

