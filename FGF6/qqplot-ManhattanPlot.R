##############################################################
install.packages("Haplin")
library("Haplin")
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

#Generate some fake data that deviates from the null

set.seed(42)
par(mfrow=c(2,2))
pvalues=runif(10000)

# NULL UNIF distribution
pvalue=data$P
data<-read.table("FG",head=T)
pvalue=data$P

Gene<-read.table("~/hpc/db/hg19/refGene.hg19.bed")
out<-data.frame(data,Gene[match(data$GENE,Gene$V5),])
input<-data.frame(CHR=out$V1,POS=out$V2,P=out$P)
input<-na.omit(input)
input$CHR=unlist(lapply(strsplit(as.character(input$CHR),split="chr"),function(x) unlist(x)[2]))
ManhattanPlot(input)

pdf("FGF6-qqplot.pdf")
pQQ(na.omit(pvalue), nlabs =200, conf = 0.95, mark = F) 
dev.off()

##############################################################

setwd("/mnt/bigdata/Genetic/Projects/shg047/hemochromatosis/Manuscript")
Gene<-read.table("~/hpc/db/hg19/refGene.hg19.bed")

ManhattanPlot<-function(mylimma){
  library(qqman)
  res <- mylimma
  SNP=res$probeID
  CHR=res$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=res$MAPINFO
  P=res$P.Value
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  pdf("manhattan.pdf")
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,9),genomewideline=-log10(genomewideline),lwd=1.5, suggestiveline=F)
  dev.off()
}

pvalue=data$P
data<-read.table("qqplot.txt",head=T)
pvalue=data$P
out<-data.frame(data,Gene[match(data$GENE,Gene$V5),])
input<-data.frame(probeID=out$GENE,CHR=out$V1,MAPINFO=out$V2,P.Value=out$P)
input<-na.omit(input)
input$CHR=unlist(lapply(strsplit(as.character(input$CHR),split="chr"),function(x) unlist(x)[2]))
head(input)
ManhattanPlot(input)


###################################################################
