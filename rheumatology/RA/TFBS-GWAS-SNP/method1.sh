bedtools sort -i  GWAS-RA-792.R2.6.rsSNP.hg19.bed > GWAS-RA-792.R2.6.rsSNP.sort.hg19.bed
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed | sort -u > GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed -b wgEncodeRegDnaseClusteredV3.hg19.bed | sort -u > GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed 
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed  -b ~/hpc/db/hg19/CpGI_Shore.hg19.bed | sort -u >  GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.hg19.bed 
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed  -b ~/hpc/db/hg19/CpGI_Shelf.hg19.bed | sort -u >  GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.hg19.bed 
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed  -b ~/hpc/db/hg19/CpGI.hg19.bed | sort -u >  GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.hg19.bed 

wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.hg19.bed 
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.hg19.bed 
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.hg19.bed 

## R
setwd("./GTEx_Analysis_v7_eQTL")
qval_threshold=0.05
file=list.files(pattern="*.v7.egenes.txt")
qval_threshold=0.05
eqtl<-c()
for(i in file){
data<-subset(read.table(i,head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(eqtl,as.character(data[,19]))
}
length(table(eqtl))
eqtl.snp<-names(table(eqtl))
write.table(eqtl.snp,file="eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)

input<-read.table("../GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.hg19.bed ",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)

input<-read.table("../GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.hg19.bed ",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)

input<-read.table("../GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.hg19.bed",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)

