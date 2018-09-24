
collect RA-GWAS-SNPS with hg38
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
R
data<-read.table("immune_disease.txt",head=T,sep="\t")
data<-read.table("RA.GWAS.SNP.txt",head=F,sep="\t")
grep("chr",data$V4)
snp150<-read.table("~/hpc/db/hg19/snp150.hg19.txt")

rlt<-unique(data.frame(paste("chr",data$V2,sep=""),data$V3-1,data$V3))
write.table(rlt,file="RA.GWAS.SNP.hg38.bed",sep="\t",col.names=F,row.names=F,quote=F)

bedtools intersect -wao -a RA.GWAS.SNP.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > RA.GWAS.SNP.hg38.commonSNP.bed
```
collect miRNA with hg38-common-SNP150
```
bedtools window -w 500000 -a ~/hpc/db/hg38/miRNA.hg38.bed -b ~/hpc/rheumatology/RA/miRNASNP/RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $4}' | sort -u | wc -l
```
32 miRNA were identified nearby 500K regions of RA-GWAS-Significant-SNPs


