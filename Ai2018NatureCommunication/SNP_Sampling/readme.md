
1. obtain alternative allele frequency table from 1000 genome dataset
```
perl printMAF.pl ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf > G1000SNP_MAF.bed
```
2. rank by the allele freqency
3. match our pre-set SNPs and get the order
4. sampling with small random numbers nearby the matched orders

```
snpdatabase<-load("/gpfs/home/guosa/hpc/db/hg19/1000Genome/G1000SNP_MAF.RData")
gwasnp<-read.table("~/hpc/db/hg19/1000Genome/GWAS-RA-378.SNPlist.txt",head=F)
snpmatch<-na.omit(match(gwasnp[,1],database[,1]))
for(i in seq(1,100)){
newsnpindex<-snpmatch+(-1)^sample(1000,length(snpmatch))*sample(1000,length(snpmatch))
newsnp<-database[newsnpindex,]
newsnp<-unlist(strsplit(as.character(newsnp[,1]),";"))
if(sum(! is.na(match(gwasnp[,1],newsnp)))>0){
  newsnp<-newsnp[-match(newsnp,gwasnp[,1])]
}
outputfile=paste("/gpfs/home/guosa/nc/gwassnpsampling/gwas-snp-sampling.",i,".txt",sep="")
write.table(newsnp,file=outputfile,sep="\t",col.names=F,row.names=F,quote=F)
}
```
5. calculate the distance and get the distribution of the SNPs to DMERs.
```
```

