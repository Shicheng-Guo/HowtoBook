
### PharmacoGenomics of Drug
```
cd /gpfs/home/guosa/hpc/rheumatology/pharmacogenomics
plink --bfile /gpfs/home/guosa/hpc/db/Hapmap/hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode --extract snplist.txt --make-bed --out PGRA

# download all the Variants from Gnomad, hg19: http://gnomad-old.broadinstitute.org/
cd ~/hpc/db/Gnomad/vcf
for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done
```
Method (1): Identify all the SNPs for 66 VIP genes
```
for i in {1..22} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view  -m2 -M2 -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \& INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R VIP.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.vcf.bgz \> VIP.chr$i.vcf >>$i.job
qsub $i.job
done
```
Method (2): Identify all the SNPs for 66 VIP genes
```
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \| INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R VIP.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
```
Step(3) # merge all the VIP SNPs with bcftools concat command
```
ls *rec.vip.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.vip.merge.vcf
```
Step(4) SNPs in the 3-UTR regions of VIP genes (241 SNPs)
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
Step(5) eQTL
```
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL

qval_threshold=0.000000000000000000000005
data1<-subset(read.table("Whole_Blood.v7.egenes.txt",head=T,sep="\t"),qval<0.000000000000000000000001*qval_threshold)
data2<-subset(read.table("Liver.v7.egenes.txt",head=T,sep="\t"),qval<10*qval_threshold)
data3<-subset(read.table("Small_Intestine_Terminal_Ileum.v7.egenes.txt",head=T,sep="\t"),qval<100*qval_threshold)
data4<-subset(read.table("Stomach.v7.egenes.txt",head=T,sep="\t"),qval<0.000000000001*qval_threshold)
eqtl<-c(as.character(data1[,19]),as.character(data2[,19]),as.character(data3[,19]),as.character(data4[,19]))
length(table(eqtl))
vip.eqtl.snp<-names(table(eqtl))
dim(data1)
dim(data2)
dim(data3)
dim(data4)
write.table(vip.eqtl.snp,file="vip.eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
```

Supplementary Materials
(S1) PBS Example
```
#PBS -N 22
#PBS -l nodes=1:ppn=1
#PBS -o /gpfs/home/guosa/hpc/db/Gnomad/vcf/temp/
#PBS -e /gpfs/home/guosa/hpc/db/Gnomad/vcf/temp/
cd /gpfs/home/guosa/hpc/db/Gnomad/vcf
bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz
bcftools view -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.001 | INFO/AF_eas>0.001 & INFO/vep ~ "missense_variant"' -R VIP.hg19.bed gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vip.vcf.bgz
bcftools sort gnomad.exomes.r2.1.sites.chr22.rec.vip.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.vcf.bgz
bcftools norm -d all gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.rmdup.vcf.bgz
bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.rmdup.biallelic.vcf.bgz
```

