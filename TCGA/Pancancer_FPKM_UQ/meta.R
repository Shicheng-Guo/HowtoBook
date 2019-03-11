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

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq")

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

##############################################################################################
##############################################################################################
##############################################################################################
library("metafor")

data<-input
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
head(rlt[order(rlt$pval),])

i=8761
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
print(md)
dev.off()


install.packages('biomaRt')
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- df$genes
df<-df[,-4]
G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")




