```
setwd("~/hpc/methylation/Pancancer/RNA-seq")
DEG<-read.table(file="TCGA-BRCA-RNAseq-FPKM-UQ.DEG.txt",sep="\t",head=T,row.names=1)
DMR<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/brca/meth450/TCGA-BRCA-BUR-PAN-TopVar.txt",head=T)
```
