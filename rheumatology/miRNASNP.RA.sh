
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
R
data<-read.table("immune_disease.txt",head=T,sep="\t")
data<-read.table("RA.GWAS.SNP.txt",head=F,sep="\t")
grep("chr",data$V4)
snp150<-read.table("~/hpc/db/hg19/snp150.hg19.txt")

rlt<-unique(data.frame(paste("chr",data$V2,sep=""),data$V3,data$V3+1))
write.table(rlt,file="RA.GWAS.SNP.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)

bedtools intersect -wao -a RA.GWAS.SNP.hg19.bed -b ~/hpc/db/hg19/snp150.hg19.txt
```

