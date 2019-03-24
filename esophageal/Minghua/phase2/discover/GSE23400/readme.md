```
GSE23400<-read.table("~/hpc/methylation/esophageal_carcinoma/GEO/GSE23400.txt",head=T,row.names = 1,sep="\t")
GSE23400$B=-GSE23400$B
GSE23400$logFC=-GSE23400$logFC
output<-GSE23400[GSE23400$Gene.symbol %in% FULL,]
write.table(output,file="~/hpc/methylation/PanCancerRNAseq/ESCS.EXP.diff.txt",sep="\t",quote=F,row.names = T,col.names = NA)
```

