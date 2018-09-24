
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
R
data<-read.table("immune_disease.txt",head=T,sep="\t")
data<-read.table("RA.GWAS.SNP.txt",head=F,sep="\t")
grep("chr",data$V4)
snp150<-read.table("~/hpc/db/hg19/snp150.hg19.txt")
```
