#!/usr/bin/R
# TCGA Methylation Pan-cancer Analysis
# Contact: Shicheng Guo
# Version 1.3
# Update: 1/4/2017
# Input: the ts-MHL counts for each samples in each reference given specific MHL positive threshold (>0.13)
# Failed: cannot predict the tissue with hypermethyalted regions in its own tissue.

# Download all the HM450K file and save them in the working directory
library("stringr")

CIMP<-function(dataframe){
  ## Functions related to CIMP (CpG island methylator phenotype)
  Q<-apply(dataframe,1,function(x) quantile(x,c(0,0.5),na.rm = T))
  names(Q)=rownames(dataframe)
  QP<-which(apply(Q,2,function(x) x[1]>0.5))
  return(QP)
}

samplesize<-function(file){
  return(table(unlist(lapply(file,function(x) unlist(strsplit(x,"[_.]"))[3]))))
}

## Start Here

file=list.files(pattern="jhu-usc.edu_*")
idv<-unique(as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*")))
idv1<-as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(any(grepl(t1,file)),any(grepl(t2,file)))){
    pairidv<-c(pairidv,t1,t2)
  }
}

# calculate the sample size for each cancer type
pairfile<-file[sapply(pairidv,function(x){x<-grep(x,file);x[length(x)]})]
samplesize(pairfile)
write.table(samplesize(pairfile),file="paired.samples.txt",col.names=NA,row.names=T,sep="\t",quote=F)

allfile<-file[sapply(idv1,function(x){x<-grep(x,file);x[length(x)]})]
samplesize(allfile)
write.table(samplesize(allfile),file="all.samples.txt",col.names=NA,row.names=T,sep="\t",quote=F)


# build cancer CIMP
file=list.files(pattern="jhu-usc.edu_LUAD")
idv<-unique(as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*")))
idv1<-as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(any(grepl(t1,file)),any(grepl(t2,file)))){
    pairidv<-c(pairidv,t1,t2)
  }
}
pairfile<-file[sapply(pairidv,function(x){x<-grep(x,file);x[length(x)]})]
samplesize(pairfile)
allfile<-file[sapply(idv1,function(x){x<-grep(x,file);x[length(x)]})]
samplesize(allfile)
data1<-c()
for(i in 1:length(pairfile)){
  tmp<-read.table(pairfile[i],head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data1<-cbind(data1,tmp[,2])
  print(c(i,pairfile[i]))
  rownames(data1)<-tmp[,1]
}

data2<-c()
for(i in 1:length(allfile)){
  tmp<-read.table(allfile[i],head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data2<-cbind(data2,tmp[,2])
  rownames(data2)<-tmp[,1]
  print(c(i,allfile[i]))
}

# CIMP: Q0>0.3 and Q50>0.6
cimp.luad<-CIMP(data1)
positive.ratio<-apply(data2[match(names(cimp.luad),rownames(data2)),],2,function(x) sum(x>0.3,na.rm=T)/length(na.omit(x)))

# validate lung and luad CIMP in all other paired samples
for(i in names(samplesize(pairfile))){
  ratio<-c()
  data3<-c()
  file=pairfile[grep(i, pairfile)]
  if(length(file)>=10){
    file<-file[1:10]
  }
  for(j in file){
  tmp<-read.table(j,head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data3<-cbind(data3,tmp[,2])
  rownames(data3)<-tmp[,1]
  #print(paste(i,j))
  }
  positive.ratio<-apply(data3[match(names(cimp.luad),rownames(data3)),],2,function(x) sum(x>0.6,na.rm=T)/length(na.omit(x)))
  ratio<-c(ratio,positive.ratio)
  print(paste(i,positive.ratio))
}
names(ratio)<-names(samplesize(pairfile))
write.table(ratio,file="luad.cimp.postive.ratio.in.paired.samples.txt",col.names=NA,row.names=T,sep="\t",quote=F)

# validate lung and luad CIMP in all other paired and non-paired samples
ratio<-c()
for(i in names(samplesize(allfile))){
  data4<-c()
  file=allfile[grep(i, allfile)]
  for(j in file){
    tmp<-read.table(j,head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
    data4<-cbind(data4,tmp[,2])
    rownames(data4)<-tmp[,1]
    print(paste(i,j))
  }
  positive.ratio<-apply(data4[match(names(cimp.luad),rownames(data4)),],2,function(x) sum(x>0.3,na.rm=T)/length(na.omit(x)))
  ratio<-c(ratio,positive.ratio)
  print(paste(i,positive.ratio))
}
names(ratio)<-names(samplesize(pairfile))
write.table(ratio,file="luad.cimp.postive.ratio.in.all.samples.txt",col.names=NA,row.names=T,sep="\t",quote=F)


