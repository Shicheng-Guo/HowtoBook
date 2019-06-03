setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/biobank/clinicaldata")
file=list.files(pattern="*.txt.uni")
data<-c()
for(i in file){
  print(i)
  temp<-read.table(i,head=F,sep="\t")
  data<-rbind(data,temp)
}

library("readxl")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/biobank/gwas")
RA<-data.frame(readxl::read_xlsx("RA-2019-GWAS.xlsx",sheet=1))
CN<-data.frame(readxl::read_xlsx("Normal-2019-GWAS.xlsx",sheet=1))
warnings()
head(RA)

setwd("/home/guosa/hpc/rheumatology/biobank/clinicaldata")
file=list.files(pattern="*.txt.uni")
data<-c()
for(i in file){
  print(i)
  temp<-read.table(i,head=F,sep="\t")
  data<-rbind(data,temp)
}

library("readxl")
setwd("/home/guosa/hpc/rheumatology/biobank/gwas")
RA<-data.frame(readxl::read_xlsx("RA-2019-GWAS.xlsx",sheet=1))
CN<-data.frame(readxl::read_xlsx("Normal-2019-GWAS.xlsx",sheet=1))
head(RA)
head(CN)

quantile(as.numeric(subset(RA,Gender=="F")[,3]))
quantile(as.numeric(subset(CN,Gender=="F")[,3]))

RA1<-subset(RA,Gender=="F")
CN1<-subset(CN,Gender=="F")
RA2<-subset(RA,Gender=="M")
CN2<-subset(CN,Gender=="M")
CN2<-CN2[head(order(CN2$Age,decreasing = T),n=183),]

data<-rbind(RA1,RA2,CN1,CN2)
write.table(data,file="GWAS-GH-RA.txt",sep="\t",quote=F,col.names = T,row.names=F)

par(cex.lab=2,cex.axis=2)
plot(density(as.numeric(subset(RA,Gender=="F")[,3])),col="blue",main="Age Distribution",lwd=2)
lines(density(as.numeric(subset(CN,Gender=="F")[,3])),col="red",lwd=2)

catlog<-read.table("https://raw.githubusercontent.com/CNAID/GWAS/master/GWAS-Catlog-Disease-Gene-Pair.20190101.txt",sep="\t")
RA<-read.table("https://raw.githubusercontent.com/CNAID/GWAS/master/RA/RA-GWAS-Gene.txt")
RAN<-catlog[catlog[,2] %in% RA[,1],]
table(RAN[,1])
head(sort(table(RAN[,1]),decreasing = T))
xx<-head(RA,n=1100)
