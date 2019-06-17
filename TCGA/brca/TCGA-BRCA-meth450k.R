##############################################################################################
##############################################################################################
##############################################################################################
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
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
########################## Read mh450K level 3 data  ########################################
##############################################################################################
setwd("/home/guosa/hpc/methylation/BRCA/450K/beta")
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
save(methdata,file="methdata.BRCA.RData")

##############################################################################################
###################### assemble methylation and phenotype data   #############################
##############################################################################################
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/meth450Pancancer.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/meth450Pancancer.R")
##################################################################################################### 
####################### Step 1: Read methylation 450K files ######################################### 
##################################################################################################### 
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
##################################################################################################### 
#######################  Step 2: Map cpg probe to GeneSymbol    ##################################### 
##################################################################################################### 
rltquantile<-apply(input,1,function(x) quantile(x,na.rm = T))
newrltquantile<-t(rltquantile[,which(rltquantile[2,]>0.6)])

map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(newrltquantile),map[,4]),],tcga_quantile=newrltquantile)
head(newrlt)
write.table(newrlt,file="/home/guosa/hpc/methylation/methbase/tcga.hyper.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

# system('wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')
# system("gunzip GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz")

data<-read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct",head=T,sep="\t",skip=2)
rnaseq<-data.matrix(data[,3:ncol(data)])
ratio<-rnaseq[,ncol(rnaseq)]/(rowMeans(rnaseq[,1:(ncol(rnaseq)-1)]))
newdata<-data.frame(data,ratio)
BOG<-subset(newdata,ratio>20 & Whole.Blood>1)
BOGHyerMeth<-merge(newrlt,BOG,by.x="V5",by.y="Description")
newBOGHyperMeth<-BOGHyerMeth[,c(2:4,1,5:ncol(BOGHyerMeth))]

write.table(newBOGHyperMeth,file="newBOGHyperMeth.hg19.bed",sep="\t",quote=F,row.names = F,col.names = F)

# bedtools intersect -wao -a newBOGHyperMeth.hg19.bed -b ../methbase/Human_PBMC.hmr.bed > newBOGHyperMeth.BUR.hg19.bed
###########################################################################################################
######################### Step 3A: Meta-analysis to differential methylation dataset #######################
###########################################################################################################
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
###########################################################################################################
######################### Step 3B: Meta-analysis to differential methylation dataset #######################
###########################################################################################################
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
write.table(normalpbmc450beta,file="TCGA-Pancancer-MH450.Meta.PBMC.diff.txt",sep="\t",quote=F,col.names=NA,row.names=T)
save(normalpbmc450beta,file="TCGA-Pancancer-MH450.Meta.PBMC.diff.RData")
###########################################################################################################
######################### Step 3B: DMR biomarker selection filtering with normal PBMC  ###################
###########################################################################################################
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("TCGA-Pancancer-MH450.Meta.PBMC.diff.RData")
hypemarker<-subset(tcgapbmc,beta>0.15 & pval<10^-7 & Quantile.75. <0.3) 
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newoutput<-data.frame(map[match(rownames(hypemarker),map[,4]),],hypemarker)
write.table(newoutput,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)
panRnaMarker<-read.table("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",head=T,sep="\t",row.names = 1)
pan5methRnaMarker<-merge(newoutput,panRnaMarker,by.x="V5",by.y="hgnc_symbol")
marker<-subset(pan5methRnaMarker,pval.x<10^-6 & pval.y<10^-6 & beta.x*beta.y<0)
dim(marker)
write.table(marker,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.FinallMarker.hg19.method1.bed",sep="\t",quote=F,col.names=NA,row.names=T)
GENE1<-unique(marker[,1])
###########################################################################################################
######################### Step 3B: DMR biomarker selection filtering with Cancer Expression ###############
###########################################################################################################
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("TCGA-Pancancer-MH450.Meta.PBMC.diff.RData")
hypemarker<-subset(tcgapbmc,beta>0.15 & pval<10^-6 & Quantile.75. <0.35) 
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newoutput<-data.frame(map[match(rownames(hypemarker),map[,4]),],hypemarker)
write.table(newoutput,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

system("bedtools intersect -wo -a TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed > TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.anno.bed")
pan5methMarkerTCGA<-read.table("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.anno.bed")
colnames(pan5methMarkerTCGA)[1:ncol(newoutput)]<-colnames(newoutput)
head(pan5methMarkerTCGA)
GENE<-unique(pan5methMarkerTCGA$V31)
panRnaMarker<-read.table("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",head=T,sep="\t",row.names = 1)
head(pan5methMarkerTCGA)
head(panRnaMarker)
pan5methRnaMarker<-merge(pan5methMarkerTCGA,panRnaMarker,by.x="V31",by.y="hgnc_symbol")
marker<-subset(pan5methRnaMarker,pval.x<10^-5 & pval.y<10^-5 & beta.x*beta.y<0)
dim(marker)
write.table(marker,file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.FinallMarker.hg19.method2.bed",sep="\t",quote=F,col.names=NA,row.names=T)
GENE2<-unique(marker[,1])
###########################################################################################################
######################### Step 3B: DMR biomarker selection filtering with other condictions   ###############
###########################################################################################################
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
MethMarker<-read.table(file="TCGA-Sig.Pancancer-MH450.Meta.PBMC.diff.hg19.FinallMarker.hg19.method2.bed",sep="\t",row.names = 1,header = T)

tsg1<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/1218.tumor.suppressor.gene.txt")
tsg2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/COSMIC.TSG.hg19.bed",sep="\t")
lung_GWAS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/lung_GWAS.Genelist.txt",sep="\t")
cellcycle<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/Gene/1180.cellcycle.gene.txt",sep="\t",head=T)
lungene<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/gene.txt")
dmrgene<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/CandidateGene.hg19.bed")
mutation<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/mostfrequentMutated.Genelist.txt")
tf<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Encode_ChIP_Seq_TF.txt")

GENE1[GENE1 %in% tsg1[,1]]
GENE2[GENE2 %in% tsg1[,1]]
GENE1[GENE1 %in% tsg2[,4]]
data.frame(GENE2[GENE2 %in% tsg2[,4]])
GENE1[GENE1 %in% lung_GWAS[,1]]
data.frame(GENE2[GENE2 %in% lung_GWAS[,1]])
data.frame(GENE2[GENE2 %in% cellcycle[,2]])
data.frame(GENE2[which(! GENE2 %in% lungene[,1])])
data.frame(GENE[GENE %in% mutation[,1]])

LEHMTF<-unique(MethMarker[MethMarker[,1]%in% tf[,1],c(1:5,11,12,14,22,37,38)])
write.table(LEHMTF,file="../TCGA_HyperMeth_LowExp_TF.txt",sep="\t",quote=F,col.names=NA,row.names=T)

###########################################################################################################
######################### Step 4: QQ-plot to shown P-value distribution   #################################
###########################################################################################################
library("Haplin")
pdf("TCGA-Pancancer-mh450k.meta.qqplot.pdf")
pQQ(rlt$pval, nlabs =nrow(rlt), conf = 0.95) 
dev.off()
###########################################################################################################
######################### Step 5: forest plot to significant DMRs   #######################################
###########################################################################################################
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
######################################################################################################
################### Step 6:Suvival analysis based on DNA methylation dataset  ########################
#####################################################################################################
library("survival")
library("survminer")
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-input[,which(id2bin(colnames(input))==1)]
newdata<-data[,na.omit(match(OS$submitter_id,id2phen3(colnames(data))))]
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
HR1<-HR
HR<-HR1
NewRowName<-paste(Symbol,rownames(newdata),sep="_")
HRV=cbind(HR,Symbol,rownames(newdata))
write.table(HRV,file="../../RNAseq_Suvival.HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)
HRVV<-read.table(file="../../RNAseq_Suvival.HR.txt",sep="\t",head=T)
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

#############################################################################################
################### Biomarker screening for BRCA  ################################
#############################################################################################
library("pROC")
library("ggplot2")

Table2Generator = function(methdata){
  seq.case = which(methdata[,1] ==1)
  seq.control = which(methdata[,1] == 11)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methdata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case],na.rm=T)$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  #Logistic regression analysis
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(ncol(methdata)-1 )){
    temp = methdata[,c(1,i+1 )]
    temp[,1] = ifelse(temp[,1] ==1,1,0)
    temp[,1] = as.factor(temp[,1])
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i]=log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
    print(i)
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}

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

methdata=data.frame(newphen$bin,t(newinput))
methdata[1:5,1:5]
newmethdata=ColNARemove(methdata)
newmethdata[1:5,1:5]
methdata<-newmethdata
rlt<-Table2Generator(newmethdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(rlt),map[,4]),],rlt)
head(newrlt)
write.table(newrlt,file="TCGA_BRCA_meth450_marker.txt",col.names=NA,row.names = T,quote=F,sep="\t")

library("survival")
library("survminer")
OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-input[,which(id2bin(colnames(input))==1)]
newdata<-data[,na.omit(match(OS$submitter_id,id2phen3(colnames(data))))]
colnames(newdata)<-id2phen3(colnames(newdata))
newdata<-RawNARemove(newdata)
phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$month=phen$time/30
head(phen)

BRCAdata<-newdata[,grep("BRCA",colnames(newdata))]
BRCAphen<-phen[grep("BRCA",colnames(newdata)),]

HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
HR1<-HR
HR<-HR1
NewRowName<-paste(Symbol,rownames(newdata),sep="_")
HRV=cbind(HR,Symbol,rownames(newdata))
write.table(HRV,file="../../RNAseq_Suvival.HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)
HRVV<-read.table(file="../../RNAseq_Suvival.HR.txt",sep="\t",head=T)


#############################################################################################
################### Biomarker screening for LUAD  ################################
#############################################################################################

setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")

install.packages("pROC")
install.packages("ggplot2")
library("pROC")
library("ggplot2")
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


Table2Generator = function(methdata){
  seq.case = which(methdata[,1] ==1)
  seq.control = which(methdata[,1] == 11)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methdata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case],na.rm=T)$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  #Logistic regression analysis
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(ncol(methdata)-1 )){
    temp = methdata[,c(1,i+1 )]
    temp[,1] = ifelse(temp[,1] ==1,1,0)
    temp[,1] = as.factor(temp[,1])
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i]=log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
    print(i)
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}


ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}

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

LUAD<-grep("LUAD",colnames(input))
newinput<-input[,LUAD]
newphen<-phen[LUAD,]
newinput[1:5,1:5]
head(newphen)
methdata=data.frame(newphen$bin,t(newinput))
methdata[1:5,1:5]
newmethdata=ColNARemove(methdata)
newmethdata[1:5,1:5]
rlt<-Table2Generator(newmethdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(rlt),map[,4]),c(1:5)],rlt)
write.table(newrlt,file="../TCGA_BRCA_meth450_marker.txt",col.names=T,row.names = F,quote=F,sep="\t")
BRCAsig<-subset(newrlt,Pvalue<10^-15)
dim(BRCAsig)


# grep -v rs ~/hpc/db/hg19/GPL13534_450K_hg19.bed > ~/hpc/db/hg19/GPL13534_450K_hg19_V2.bed
# bedtools intersect -wa -a ~/hpc/db/hg19/GPL13534_450K_hg19_V2.bed -b ~/hpc/methylation/methbase/Human_PBMC.hmr.bed > ~/hpc/db/hg19/GPL13534_450K_hg19_PBMC_BUR.bed
BRCAsig<-read.table(file="~/hpc/methylation/TCGA_BRCA_meth450_marker.txt",head=T,sep="\t")
map<-read.table("~/hpc/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
BURCpG<-read.table("~/hpc/db/hg19/GPL13534_450K_hg19_PBMC_BUR.bed",sep="\t")

BRCAsig<-subset(BRCAsig,Pvalue<10^-10 & AUC>0.85 & OR>0)
dim(BRCAsig)
BRCAsigBUR<-merge(BRCAsig,BURCpG,by="V4")
dim(BRCAsigBUR)

UKBiobankBRCA<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/brca/UKbiobank.pheWas.txt",sep="\t")
COSMICTSG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/COSMIC.TSG.hg19.bed",sep="\t")
FDAdrugGene<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/FDA_approved_drugtarget.txt",sep="\t",head=T)
PanCancerDMG<-subset(read.table("~/hpc/methylation/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",sep="\t",head=T,row.names = 1),pval<10^-8 & beta<0)
TSG2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/TSGene2.0.txt",sep="\t",head=T,row.names = 1)
OverallSuvival<-subset(read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_FPKM_UQ/TCGA_RNAseq_Suvival.HR.txt",sep="\t",head=T),exp.coef.<1 & Pr...z..<10^-5)
BRCAOSFavarible<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/brca/TCGA_BRCA_OverallSuvival_favorible_genes.txt",head=T,sep="\t")

SA<-unique(as.character(BRCAsigBUR[BRCAsigBUR[,5] %in% UKBiobankBRCA$V2,5]))
SB<-unique(as.character(BRCAsigBUR[BRCAsigBUR[,5] %in% COSMICTSG$V4,5]))
SC<-unique(as.character(BRCAsigBUR[BRCAsigBUR[,5] %in% FDAdrugGene[,1],5]))
SD<-unique(as.character(BRCAsigBUR[BRCAsigBUR[,5] %in% PanCancerDMG$hgnc_symbol,5]))
SE<-unique(as.character(BRCAsigBUR[BRCAsigBUR[,5] %in% TSG2$GeneSymbol,5]))
SF<-unique(as.character(BRCAsigBUR[BRCAsigBUR[,5] %in% BRCAOSFavarible$Gene,5]))

sort(table(c(SA,SB,SC,SD,SE)))
which(table(c(SA,SB,SC,SD,SE))>=2)
names(which(table(c(SA,SB,SC,SD,SE,SF))>=2))
length(which(table(c(SA,SB,SC,SD,SE))>=2))
SF[SF %in% names(which(table(c(SA,SB,SC,SD,SE))>=2))]
BRCAsigBURBiomarker<-BRCAsigBUR[match(names(which(table(c(SA,SB,SC,SD,SE,SF))>=2)),BRCAsigBUR[,5]),]
write.table(BRCAsigBURBiomarker,file="../TCGA_BRCAsigBUR_meth450_marker.txt",col.names=T,row.names = F,quote=F,sep="\t")



newrlt<-data.frame(map[match(rownames(rlt),map[,4]),c(1:5)],rlt)
head(newrlt)
write.table(newrlt,file="../TCGA_LUAD_meth450_marker.txt",col.names=T,row.names = F,quote=F,sep="\t")

map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(rlt),map[,4]),c(1:5)],rlt)
head(newrlt)
write.table(newrlt,file="../TCGA_LUSC_meth450_marker.txt",col.names=T,row.names = F,quote=F,sep="\t")

TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt

#############################################################################################
################### Biomarker screening for LGG/GBM  ################################
#############################################################################################
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

Seq<-paste(phen$pid,phen$bin,sep="-")
head(phen)
input[1:5,1:5]
BC<-grep("LGG|GBM",colnames(input))
newinput<-input[,BC]
newphen<-phen[BC,]
newinput[1:5,1:5]
head(newphen)

load("~/GSE41826.RData")
load("~/GSE89707.RData")
load("~/GSE74486.RData")
load("~/GSE103659.RData")
load("~/GSE114534.RData")
load("~/GSE66351.RData")

GSE41826$beta<-GSE41826$beta[match(rownames(newinput),rownames(GSE41826$beta)),]
GSE89707$beta<-GSE89707$beta[match(rownames(newinput),rownames(GSE89707$beta)),]
GSE74486$beta<-GSE74486$beta[match(rownames(newinput),rownames(GSE74486$beta)),]
GSE103659$beta<-GSE103659$beta[match(rownames(newinput),rownames(GSE103659$beta)),]
GSE114534$beta<-GSE114534$beta[match(rownames(newinput),rownames(GSE114534$beta)),]
GSE66351$beta<-GSE66351$beta[match(rownames(newinput),rownames(GSE66351$beta)),]

dim(GSE41826$beta)
dim(GSE89707$beta)
dim(GSE74486$beta)
dim(GSE103659$beta)
dim(GSE114534$beta)
dim(GSE66351$beta)

GSE41826$beta
GSE89707$beta
GSE74486$beta
GSE103659$beta
GSE114534$beta
GSE66351$beta

Name=c("disease","tissue","gender","cell type")
GSE41826$phen[,c(38,41,36,37)]
GSE89707$phen[,c(35,36,37)]
GSE74486$phen[,c(43:47)]
GSE103659$phen[,c(40:47)]
GSE114534$phen[,c(34:37)]
GSE66351$phen[,c(41:49)]

input<-data.frame(newinput,GSE41826$beta,GSE89707$beta,GSE74486$beta,GSE103659$beta,GSE66351$beta,GSE114534$beta)


data=data.frame(newphen$bin,t(newinput))
data[1:5,1:5]
newdata=ColNARemove(data)
newdata[1:5,1:5]
rlt<-Table2Generator(newdata)
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(rlt),map[,4]),c(1:5)],rlt)
write.table(newrlt,file="../TCGA_BRCA_meth450_marker.txt",col.names=T,row.names = F,quote=F,sep="\t")
BRCAsig<-subset(newrlt,Pvalue<10^-15)
dim(BRCAsig)

