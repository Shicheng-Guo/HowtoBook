
load("MCRI.GBM.brain.RData")
beta<-GBM$beta
phen<-GBM$phen
beta[1:5,1:5]
phen[1:5,1:3]
phen$disease=as.character(phen$disease)
phen$tissue=as.character(phen$tissue)

load("~/hpc/methylation/Normal.PBMC.GEO.HM450K.Beta.RData")
BUR<-subset(normalpbmc450beta,mean<0.1 & Quantile.75.<0.1)
beta<-beta[match(rownames(BUR),rownames(beta)),]
beta<-data.matrix(beta)
lapply(beta,function(x) aov(x ~ phen$disease))

GBM<-list()
GBM$beta<-beta
GBM$phen<-phen
save(GBM,file="MCRI.GBM.MachineLearning.RData")


