library("Haplin")
library("qqman")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/SSc/FGF6")


data<-read.table("Iron.Gene.in.SSc.Trim.Pvalue.txt",head=T,row.names = 1)
bed<-read.table("Iron.Gene.in.SSc.hg19.bed",head=F)

input<-data.frame(data,bed[match(data[,6],bed[,6]),])
mylimma<-data.frame(probeID=input$V6,CHR=input$V1,MAPINFO=input$V2,P=input$P.Value)
ManhattanPlot(mylimma,annotatePval=0.00001,cex=2,snpsOfInterest="HAMP")

pQQ(data[,2], nlabs =nrow(data), conf = 0.95) 
head(data)
library("calibrate")

ManhattanPlot<-function(mylimma,annotatePval,cex,snpsOfInterest){
  library(qqman)
  res <- mylimma
  SNP=res$probeID
  CHR=res$CHR
  CHR=gsub("chr","",CHR)
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=res$MAPINFO
  P=res$P
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  manhattan(manhattaninput,col = c("blue4", "orange3"),highlight = snpsOfInterest,cex=cex,suggestiveline=-log(0.00001,10),genomewideline=F,ylim = c(0,12),lwd=1, annotatePval=annotatePval,annotateTop = FALSE)
}

annotatePval=0.00001;
cex=1.5;
snpsOfInterest="HAMP";

pdf("manhattan.pdf",width=864, height=631)
manhattan(manhattaninput,col = c("blue4", "orange3"),highlight = snpsOfInterest,cex=3,suggestiveline=-log(0.00001,10),genomewideline=F,ylim = c(0,12),lwd=1, annotatePval=annotatePval,annotateTop = FALSE)
dev.off()

