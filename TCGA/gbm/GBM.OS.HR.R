
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
methdata[1:5,1:5]
phen4<-id2phen4(colnames(methdata))
phen3<-id2phen3(colnames(methdata))
bin<-id2bin(colnames(methdata))
pid<-id2pid(colnames(methdata))
phen<-data.frame(phen4=phen4,phen3=phen3,pid=pid,bin=bin)
exclude<-which(c(phen$bin !=1 & phen$bin !=11))
phen<-phen[-exclude,]
input<-methdata[,-exclude]
Seq<-paste(phen$pid,phen$bin,sep="-")
head(phen)
input[1:5,1:5]

BRCA<-grep("BRCA",colnames(input))
newinput<-input[,BRCA]
newphen<-phen[BRCA,]
newinput[1:5,1:5]
library("survival")
library("survminer")
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-newinput[,which(id2bin(colnames(newinput))==1)]
newdata<-data[,id2phen3(colnames(data)) %in% OS$submitter_id]
colnames(newdata)<-id2phen3(colnames(newdata))
newdata<-RawNARemove(newdata)
phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$month=phen$time/30
head(phen)

HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=0.3]<-0
  dat$Rna[dat$Rna>0.3]<-1
  hr<-summary(coxph(Surv(month,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
rownames(HR)<-rownames(newdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
rlt<-data.frame(HR,map[match(rownames(HR),map[,4]),])
write.table(rlt,file="~/hpc/methylation/TCGA_HM450_GBM_Suvival_HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)


