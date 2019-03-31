GSE112658: 10 RA vs 10 OA
```
setwd("/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/GSE112658/RNAseq")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication/GSE112658/RNAseq")
ENST2SymbolRef<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
head(ENST2SymbolRef)
data<-read.table("GSE112656_Counts_matix.txt",head=T,row.names=1,sep="\t")
head(rownames(data))
ENSG<-unlist(lapply(strsplit(as.character(rownames(data)),split="[.]"),function(x) x[1]))
head(ENSG)
Rownames<-ENST2SymbolRef[match(ENSG,ENST2SymbolRef$V7),5]
newdata=data/colMeans(data)
head(newdata)

P=apply(newdata,1,function(x) t.test(x[1:10],x[11:20],na.rm=T)$p.value)
head(P)
output<-data.frame(Rownames,P)
head(output)
head(rlt,50)

target<-c("FSTL1","FSTL3","IQUB","CDON","FGF10","IL4","TNF","MYD88","IL6")
FSTL<-data.frame(t(newdata[match(target,output$Rownames),]),type=c(rep("RA",10),rep("OA",10)))
colnames(FSTL)<-c(target,"type")
head(FSTL)
par(mfrow=c(3,3))
for(i in 1:length(target)){
  print(i)
  boxplot(FSTL[,i]~type,FSTL,col=c("red","blue"),ylab=target[i])
}
```
