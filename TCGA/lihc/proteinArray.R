#####################################################################################
id2phen4<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
}

id2phen3<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*"))
}

id2bin<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
  as.numeric(lapply(strsplit(filename,"-"),function(x) x[4]))
}

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}

RawZeroRemove<-function(data,missratio=0.5){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) sum(x==0)>=threshold))
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}
###################################################################################################################
##################### ENSG, ENST, ENSP, Symbol in hg19 ############################################################
###################################################################################################################
setwd("~/hpc/methylation/Pancancer/RNA-seq")
load("rnaseqdata.pancancer.RData")
phen1=read.table("/home/guosa/hpc/methylation/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("/home/guosa/hpc/methylation/Pancancer/RNA-seq/File_metadata2.txt",header = T,sep="\t")
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
head(phen1)
head(phen2)
head(phen)

colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]

phen1<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen2<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen3<-id2bin(phen$cases.0.samples.0.submitter_id)

phen$bin=phen3

include<-which(c(phen$bin==1 | phen$bin==11))
phen<-phen[include,]
input<-rnaseqdata[,include]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]
colnames(input)<-phen$id
Seq<-paste(phen$project_id,phen$bin,sep="-")

library("metafor")
data<-input
data[1:5,1:5]
i=500
Seq<-paste(phen$project_id,phen$bin,sep="-")
mean<-tapply(as.numeric(data[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(data[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(data[i,]),Seq,function(x) length(x))

exclude<-names(which(table(unlist(lapply(strsplit(names(mean),"-"),function(x) x[2])))<2))
exclude <-grep(paste(exclude,collapse="|"),phen$project_id)

phen<-phen[-exclude,]
input<-input[,-exclude]
phen$id=id2phen4(phen$cases.0.samples.0.submitter_id)
colnames(input)<-phen$id
dim(phen)
dim(input)
input[1:5,1:5]
phen[1:5,1:5]

input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$bin,sep="-")

rlt<-c()
coll<-c()
for(i in 1:nrow(input)){
  print(i)
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
}
rownames(rlt)<-rownames(input)[coll]
colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)

write.table(rlt,file="TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)
save.image("TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.RData")


SEL<-grep("LIHC",Seq)
input<-input[,SEL]
Seq<-Seq[SEL]
yphen<-abs(as.numeric(as.factor(Seq))-2)
rlt<-c()
for(i in 1:nrow(input)){
  fit<-summary(glm(yphen~input[i,]))$coefficients[2,]
  rlt<-rbind(rlt,fit)
}
rownames(rlt)<-rownames(input)
rlt2<-read.table(file="TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.txt",sep="\t",head=T,row.names=1)

ENST2Symbol<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
ENSG<-unlist(lapply(strsplit(as.character(rownames(newrlt)),split="[.]"),function(x) x[1]))
Symbol<-ENST2Symbol[match(ENSG,ENST2Symbol$V7),5]
newrlt<-cbind(rlt,rlt2,Symbol)
head(newrlt[order(newrlt[,4]),],n=50)

theratarget<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/therapeutic_targets_603.txt")
NewRlt<-subset(newrlt,newrlt[,4]<10^-6 & newrlt[,7]<10^-6 & newrlt[,1]*newrlt[,6]>0)
NewRlt<-NewRlt[na.omit(match(theratarget[,1],NewRlt$Symbol)),]
write.table(NewRlt,file="~/hpc/methylation/TCGA-LIHC-FPKM-UQ.DGE.txt",sep="\t",quote=F,col.names=NA,row.names=T)


library("survival")
library("survminer")
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
input<-input[,which(id2bin(colnames(input))==1)]
newdata<-input[,na.omit(match(OS$submitter_id,id2phen3(colnames(input))))]
colnames(newdata)<-id2phen3(colnames(newdata))
phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$week=phen$time/7
head(phen)
ENST2Symbol<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
ENSG<-unlist(lapply(strsplit(as.character(rownames(newdata)),split="[.]"),function(x) x[1]))
Symbol<-ENST2Symbol[match(ENSG,ENST2Symbol$V7),5]
HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
rownames(HR)=rownames(newdata)
NewRltHR<-cbind(rlt,rlt2,Symbol,HR)


theratarget<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/therapeutic_targets_603.txt")
New<-subset(NewRltHR,NewRltHR[,4]<10^-9 & NewRltHR[,7]<10^-9 & NewRltHR[,17]<10^-3 & NewRltHR[,1]*NewRltHR[,6]>0 & NewRltHR[,1]*NewRltHR[,13]>0)
New<-New[New$Symbol %in% theratarget[,1],]
write.table(New,file="~/hpc/methylation/TCGA-LIHC-FPKM-UQ.DGE-Overall_82.txt",sep="\t",quote=F,col.names=NA,row.names=T)
