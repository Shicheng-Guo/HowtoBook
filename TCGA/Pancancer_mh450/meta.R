install.packages("metafor")
library("metafor")
rm(list=ls())
load("CD4_RGSS_data_beforecombat.RData")
phen<-read.csv("CD4_RGSS_clinical.csv")
methdata<-CD4_RGSS_data_beforecombat
dim(phen)
dim(methdata)
head(phen)

i=500
Seq<-paste(phen[,3],phen[,2],sep="_")
mean<-tapply(as.numeric(methdata[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(methdata[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(methdata[i,]),Seq,function(x) length(x))

m1i=mean[seq(1,8,by=2)]
m2i=mean[seq(2,8,by=2)]
sd1i=sd[seq(1,8,by=2)]
sd2i=sd[seq(2,8,by=2)]
n1i=num[seq(1,8,by=2)]
n2i=num[seq(2,8,by=2)]
Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
data<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
data$source=Source
data
es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=data)
md <- rma(es,slab=source,method = "REML", measure = "MD",data=data)
plot(md)
