setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/pmrp/phase2/RA/C2")
library("Haplin")
data<-read.table("MCRI_RA_C2_Exom2.hg19.bed",head=F,sep="\t")
data<-data[order(data[,10],decreasing = F),]
head(data)
input<-unique(data[,c(1:4,10)])
head(input)
par(cex.lab=1.5,cex.axis=1.5,cex=1)
pQQ(input[,5], nlabs =nrow(input), conf = 0.95,lim=c(0,4.5))
write.table(data,file="MCRI_RA_C2_Exom2_PvalueSort.hg19.bed",sep="\t",quote=F,col.names =F,row.names = F)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/pmrp/phase1/plink/RA")
library("Haplin")
data<-read.table("MCRI_RA_C2_Exom2.hg19.bed",head=F,sep="\t")
data<-data[order(data[,10],decreasing = F),]
head(data)
input<-unique(data[,c(1:3,10)])
head(input)
par(cex.lab=1.5,cex.axis=1.5,cex=1)
pQQ(input[,4], nlabs =nrow(input), conf = 0.95,lim=c(0,9))
write.table(data,file="MCRI_RA_C2_Exom1_PvalueSort.hg19.bed",sep="\t",quote=F,col.names =F,row.names = F)



setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/pmrp/phase2/RA/C2")
library("Haplin")
data<-read.table("MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta",head=T,sep="")
data<-data[order(data[,8],decreasing = F),]
head(data)
input<-unique(data[,c(1:3,7,8)])
head(input)
par(cex.lab=1.5,cex.axis=1.5,cex=1)
pQQ(input[,4], nlabs =nrow(input), conf = 0.95,lim=c(0,9))
pQQ(input[,5], nlabs =nrow(input), conf = 0.95,lim=c(0,9))
write.table(data,file="MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.bed",sep="\t",quote=F,col.names =T,row.names = F)

meta<-read.table("MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta",head=T,sep="")
head(meta)
ManhattanPlot(meta)
mylimma=meta


#####################################################################################
#########################    Psoriatic Arthritis    ################################
#####################################################################################
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/pmrp/phase1/plink")
phen<-read.table("FinalRelease_QC_20140311_Team1_Marshfield.phen",sep="\t",head=T)
head(phen)
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/pmrp/phase2/PA/")
meta<-read.table("MCRI.PheTyp1_PA_C2_Exom1_Exom2.meta",head=T)
meta<-meta[order(meta[,8],decreasing = F),]
write.table(meta,file="MCRI.PheTyp1_PA_C2_Exom1_Exom2.meta.Sortpvalue.txt",sep="\t",col.names=T,row.names = F,quote=F)
head(meta)
mylimma<-meta
library(qqman)
res <- mylimma
SNP=res$SNP
CHR=res$CHR
if(length(grep("X",CHR))>0){
  CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
  CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
}
CHR<-as.numeric(CHR)
BP=res$BP
P=res$P
manhattaninput=data.frame(SNP,CHR,BP,P)
max<-max(2-log(manhattaninput$P,10))
genomewideline=0.05/nrow(manhattaninput)
par(cex.lab=1.25,cex.axis=1.25)
manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),genomewideline=F,lwd=2, suggestiveline=F)

ManhattanPlot<-function(mylimma){
  library(qqman)
  res <- mylimma
  SNP=res$SNP
  CHR=res$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=res$BP
  P=res$P
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  pdf("manhattan.pdf")
  par(cex.lab=1.25,cex.axis=1.25)
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),genomewideline=F,lwd=2, suggestiveline=F)
  dev.off()
}

getwd()
