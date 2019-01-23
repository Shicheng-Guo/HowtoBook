#!/usr/bin/R
# TCGA Methylation Pan-cancer Analysis
# Contact: Shicheng Guo
# Version 1.3
# Update: 01/20/2018
# Input: the ts-MHL counts for each samples in each reference given specific MHL positive threshold (>0.13)
# Failed: cannot predict the tissue with hypermethyalted regions in its own tissue.
# Download all the HM450K file and save them in the working directory
setwd("/home/local/MFLDCLIN/guosa/hpc/tcga/xiong/luad/miRNA")
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

manifest="gdc_manifest.2019-01-20.txt"
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > barcode.txt")


file=list.files(pattern="*mirnas.quantification.txt")
barcode<-read.table("barcode.txt",sep="\t",head=T)
tsz=samplesize(file)

data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=F,sep="\t",as.is=F)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data<-cbind(data,tmp[,2])
  print(paste(i,"in",tsz,file[i],sep=" "))
  rownames(data)<-tmp[,1]
}

barcodename<-barcode[match(file,barcode$file_name),ncol(barcode)]
colnames(data)<-barcodename
save(data,file="TCGA-LUAD.miRNAseq.RData")
