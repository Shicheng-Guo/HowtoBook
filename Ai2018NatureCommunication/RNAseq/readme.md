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
head(data)
P=apply(data,1,function(x) t.test(x[1:10],x[11:20],na.rm=T)$p.value)
head(P)
output<-data.frame(Rownames,P)
head(output)

FSTL1<-data.frame(fstl1=t(data[match("FSTL1",output$Rownames),]),type=c(rep("RA",10),rep("OA",10)))
FSTL1
boxplot(ENSG00000163430.5~type,FSTL1)
```
