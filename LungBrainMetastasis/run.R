install.packages("openxlsx")
library(openxlsx)
library("devtools")
library("GSEABase")

devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)

ci95<-function(x){
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  m<-round(mean(x),2)
  d<-round(mean(x)-error,2)
  u<-round(mean(x)+error,2)
  paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
}

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis")

data<-read.xlsx("Result.xlsx",sheet=6,rowNames=T)

ci95(data[,9])
ci95(data[seq(1,26,by=2),9])
ci95(data[seq(2,26,by=2),9])

data<-read.xlsx("Result.xlsx",sheet=7,rowNames=T)
ci95(data[1:14,8])

data<-read.xlsx("Result.xlsx",sheet=4,rowNames=F)
db<-read.xlsx("Result.xlsx",sheet=9,rowNames=F)

temp<-db[match(data[,1],as.character(db[,1])),]
out=na.omit(data.frame(data,temp))
out[1:5,1:5]
mean(apply(out[,2:27],2,function(x) sum(x>0)))
write.table(out,file="LungBrain-Pancancer.txt",sep="\t",quote=F)


data<-read.xlsx("Result.xlsx",sheet=4,rowNames=F)
head(data,n=30)
par(las=1,cex.axis=0.6)
barplot(data$Sum2[1:30]/26,col="red",names.arg=data[1:30,1],horiz = T,xlim=c(0,0.5),xlab="Freqency of Mutation")

# grep -v '#' *.vcf | perl -lane '{print "$1>$2" if $_=~/(\w)>(\w)/}' > sub.txt
data<-read.table("sub.txt")
head(data)

pdf("Figure2.pdf")
par(las=1,cex.axis=1)
barplot(table(data[,1]),col="red",horiz = T,xlab="Counts of different Mutation type")
dev.off()
write.table(table(data[,1]),file="Transition-transversion.txt",sep="\t",quote=F)


data<-read.xlsx("Result.xlsx",sheet=10,rowNames=F)
head(data,n=30)
pdf("pathway.enrichment.pdf")
par(las=2,cex.axis=0.5)
barplot(data$Fold.Enrichment,col="red",names.arg=data[,2],ylab="Fold of Enrichment")
dev.off()

data<-read.xlsx("Result.xlsx",sheet=11,rowNames=F)
head(data,n=30)
pdf("Kewwords.enrichment.pdf")
par(las=2,cex.axis=0.75)
barplot(data$Fold.Enrichment,col="red",names.arg=data[,2],ylab="Fold of Enrichment",ylim=c(0,3))
dev.off()




