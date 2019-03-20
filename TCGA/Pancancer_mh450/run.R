##############################################################################################
##############################################################################################
##############################################################################################
manifest="gdc_manifest.2019-03-09.txt"
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > file_metadata.txt")

##############################################################################################
##############################################################################################
##############################################################################################

setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
files=list.files(pattern="*gdc_hg38.txt$",recursive = T)
methdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=T,sep="\t",row.names = 1)
  methdata<-cbind(methdata,temp[,1])
  print(i)
}
colnames(methdata)<-files
rownames(methdata)<-rownames(temp)
methdata[1:5,1:5]
save(methdata,file="methdata.pancancer.RData")

##############################################################################################
##############################################################################################
##############################################################################################

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

id2pid<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"edu_...."))
  unlist(lapply(filename,function(x) unlist(strsplit(x,"[_]"))[2]))
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
##############################################################################################
##############################################################################################
##############################################################################################
library("metafor")
data<-input
i=500
Seq<-paste(phen$pid,phen$bin,sep="-")
mean<-tapply(as.numeric(data[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(data[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(data[i,]),Seq,function(x) length(x))
exclude<-names(which(table(unlist(lapply(strsplit(names(mean),"-"),function(x) x[1])))<2))
if(length(exclude)>0){
exclude <-grep(paste(exclude,collapse="|"),phen$pid)
length(exclude)
head(exclude)
phen<-phen[-exclude,]
input<-input[,-exclude]
colnames(input)<-phen$phen4
}
dim(phen)
dim(input)
input[1:3,1:3]
head(phen)
#############################################################################################
#############################################################################################
#############################################################################################
newinput<-RawNARemove(input)
Seq<-paste(phen$pid,phen$bin,sep="-")
newinput[1:3,1:3]
head(phen)
rlt<-c()
coll<-c()
for(i in 1:nrow(newinput)){
  print(i)
  mean<-tapply(as.numeric(newinput[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(newinput[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(newinput[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
  output<-data.frame(cbind(n1i,m1i,sd1i ,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  rlt<-rbind(rlt,c(i,C=mean(m1i,na.rm=T),N=mean(m2i,na.rm=T),md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
  coll<-c(coll,i)
}
rownames(rlt)<-rownames(newinput)[coll]
colnames(rlt)<-c("idx","C","N","beta","pval","cilb","ciub","i2","tau2")
rlt<-data.frame(rlt)
write.table(rlt,file="TCGA-Pancancer-MH450.Meta.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)


load("/mnt/bigdata/Genetic/Projects/shg047/methylation/GEO/Normal.PBMC.GEO.HM450K.Beta.RData")

write.table(tcga,file="TCGA-Pancancer-MH450.Meta.PBMC.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)
save(tcga,file="TCGA-Pancancer-MH450.Meta.PBMC.diff.RData")

load("TCGA-Pancancer-MH450.Meta.PBMC.diff.RData")
hypemarker<-subset(tcgapbmc,beta>0.15 & pval<10^-15 & Quantile.75. <0.3) 
dim(hypemarker)
head(hypemarker)
hypomarker<-subset(tcgapbmc,beta< -0.15 & pval<10^-15 & Quantile.75.>0.6) 
dim(hypomarker)
marker<-rbind(hypemarker,hypomarker)
# temp<-read.table("ffc283aa-7924-4b1d-b9e7-d156510ae5f9/jhu-usc.edu_COAD.HumanMethylation450.14.lvl-3.TCGA-NH-A5IV-01A-42D-A36Y-05.gdc_hg38.txt",head=T,sep="\t",row.names = 1)
output<-data.frame(temp[match(rownames(marker),rownames(temp)),2:4],marker)
output[,2]<-output[,2]-250
output[,3]<-output[,3]+250
write.table(output,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)

map<-read.table("~/hpc/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newoutput<-data.frame(map[match(rownames(output),map[,4]),],output)
write.table(newoutput,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)
dim(newoutput)
head(newoutput)

bedtools intersect -wo -a TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed > TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.anno.bed

pan5methMarkerTCGA<-read.table("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.anno.bed")
colnames(pan5methMarkerTCGA)[1:ncol(newoutput)]<-colnames(newoutput)
panRnaMarker<-read.table("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",head=T,sep="\t",row.names = 1)
head(pan5methMarkerTCGA)
head(panRnaMarker)

pan5methRnaMarker<-merge(pan5methMarkerTCGA,panRnaMarker,by.x="V34",by.y="hgnc_symbol")
marker<-subset(pan5methRnaMarker,pval.x<10^-6 & pval.y<10^-6 & beta.x*beta.y<0)
dim(marker)
write.table(marker,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.FinallMarker_E_6.bed",sep="\t",quote=F,col.names=NA,row.names=T)
unique(marker[,1])


################################################################################################################################
################################################################################################################################
library("Haplin")
pdf("TCGA-Pancancer-mh450k.meta.qqplot.pdf")
pQQ(rlt$pval, nlabs =nrow(rlt), conf = 0.95) 
dev.off()
################################################################################################################################
################################################################################################################################
################################################################################################################################
Sig<-head(rlt[order(rlt$pval),],n=100)
Sig<-subset(Sig,abs(beta)>0.15)
new<-data.frame(temp[match(rownames(Sig),rownames(temp)),2:4],Sig)
head(new)
write.table(new,file="TCGA-Sig.Pancancer-MH450.Meta.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)

for(i in match(rownames(marker),rownames(newinput))){
  print(i)
  mean<-tapply(as.numeric(newinput[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(newinput[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(newinput[i,]),Seq,function(x) length(x))
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
  pdf(paste(rownames(newinput)[i],".pdf",sep=""))
  plot(md)
  dev.off()
}
#############################################################################################
################### Parallel Solution for DMR analysis ################################
#############################################################################################
metaDMR<-function(input,Seq){
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
    output<-data.frame(cbind(n1i,m1i,sd1i ,n2i,m2i,sd2i))
    output$source=Source
    output<-na.omit(output)
    es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
    md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
    rlt<-rbind(rlt,c(i,C=mean(m1i,na.rm=T),N=mean(m2i,na.rm=T),md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
    coll<-c(coll,i)
  }
  rownames(rlt)<-rownames(input)[coll]
  colnames(rlt)<-c("idx","beta","pval","cilb","ciub","i2","tau2")
  rlt<-data.frame(rlt)
}

library(foreach)
library(doParallel)
registerDoParallel(20) 

rlt<-foreach (i=seq(1,nrow(newinput,by=45000)),.combine=rbind) %do% {
  metaDMR(newinput,Seq)
}


