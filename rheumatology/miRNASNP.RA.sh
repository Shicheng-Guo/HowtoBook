
collect RA-GWAS-SNPS with hg38
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
R
data<-read.table("RA.GWAS.SNP.txt",head=F,sep="\t")
grep("chr",data$V4)
snp150<-read.table("~/hpc/db/hg19/snp150.hg19.txt")

rlt<-unique(data.frame(paste("chr",data$V2,sep=""),data$V3-1,data$V3))
write.table(rlt,file="RA.GWAS.SNP.hg38.bed",sep="\t",col.names=F,row.names=F,quote=F)

bedtools intersect -wao -a RA.GWAS.SNP.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > RA.GWAS.SNP.hg38.commonSNP.bed
```
collect miRNA with hg38-common-SNP150
```

bedtools intersect -wo -a ~/hpc/db/hg38/miRNA.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > miRNA.hg38.commonSNP.bed
bedtools window -w 500000 -a miRNA.hg38.commonSNP.bed -b ~/hpc/rheumatology/RA/miRNASNP/RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.hg38.commonSNP.bed -b ~/hpc/rheumatology/RA/miRNASNP/RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $8}' | sort -u
```
32 miRNA were identified nearby 500K regions of RA-GWAS-Significant-SNPs

Obviously, we don't have too many SNPs to be test for RA. in order to make the target to 50, we increase the diease to all immune diease.
```
data<-read.table("immune_disease.txt",head=T,sep="\t")
rlt<-na.omit(unique(data.frame(paste("chr",data$CHR_ID,sep=""),as.numeric(as.character(data$CHR_POS))-1,as.numeric(as.character(data$CHR_POS)))))
write.table(rlt,file="AutoImmue.GWAS.SNP.hg38.bed",sep="\t",col.names=F,row.names=F,quote=F)

bedtools intersect -wo -a ~/hpc/db/hg38/miRNA.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > miRNA.hg38.commonSNP.bed
bedtools intersect -wo -a miRNA.hg38.commonSNP.bed -b ~/hpc/db/hg38/AutoImmue.GWAS.SNP.hg38.bed > AutoImmue.GWAS.SNP.hg38.commonSNP.bed
```
