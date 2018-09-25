### A Gene-Based Recessive Diplotype Exome Scan Discovers Novel Type II diabetes and Obesity Genes
```
cd /gpfs/home/guosa/hpc/T2D
list.files()
data1<-read.table("T2D.Candidate.Gene.txt")
data2<-read.table("T2D.GWAS.Gene.txt")
rlt<-data.frame(data1,data1[,1] %in% data2[,1])
write.table(rlt,file="T2D.Candidate.Gene.GWAS.txt",sep="\t",quote=F)
data3<-read.table("Obesity.Candidate.Gene.txt")
data4<-read.table("Obesity.GWAS.Gene.txt")
rlt<-data.frame(data3,data3[,1] %in% data4[,1])
write.table(rlt,file="Obesity.Candidate.Gene.GWAS.txt",sep="\t",quote=F)
rlt<-data.frame(data1,data1[,1] %in% data3[,1])
```

