
```
miRdata<-read.table("/home/guosa/hpc/db/mirTarget/Predicted_Target_Locations.default_predictions.hg19.bed")
RAGene<-read.table("/home/guosa/hpc/rheumatology/RA/KEGG_90_RA_GeneList.txt")
miR<-miRdata[miRdata[,4] %in% RAGene[,1],]
write.table(miR,file="miRNA_KEGG_90Gene.txt",sep="\t",quote=F,col.names = F,row.names=F)
```
