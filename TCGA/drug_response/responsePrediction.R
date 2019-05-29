setwd("/home/guosa/hpc/project/TCGA")

phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
barcode<-read.table("/home/guosa/hpc/project/TCGA/pancancer/FPKM/barcode.txt",head=T,sep="\t")
load("rnaseqdata.pancancer.env.RData")

ncn<-barcode[match(unlist(lapply(colnames(rnaseqdata),function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name<-gsub(".gz","",barcode$file_name)),]
ncol<-match("cases.0.samples.0.submitter_id",colnames(ncn))
colnames(rnaseq)<-ncn[,ncol]
phen$ID<-paste(phen$bcr_patient_barcode,"-01",sep="")

