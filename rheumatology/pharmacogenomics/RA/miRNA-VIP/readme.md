
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/miRNA-SNP

wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

vip.gene<-read.table("/gpfs/home/guosa/hpc/db/Gnomad/vcf/VIP.hg19.bed")
data1<-read.table("miRNA_targets_loss_by_SNPs_in_seed_regions.txt",head=F,sep="\t")
vip.variant1<-data1[data1[,2]%in%vip.gene[,4],]
data2<-read.table("miRNA_targets_gain_by_SNPs_in_seed_regions.txt",head=T,sep="\t")
vip.variant2<-data2[data2[,2]%in%vip.gene[,4],]


colnames(vip.variant1)[2:5]=colnames(vip.variant2[,2:5])
vip.variant<-rbind(vip.variant1[,2:5],vip.variant2[,2:5])
length(table(vip.variant[,4]))
vip.utr.mirna.snp<-names(table(vip.variant[,4]))
write.table(vip.utr.mirna.snp,file="vip.utr.mirna.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
```

