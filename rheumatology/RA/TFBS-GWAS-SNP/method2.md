S1: R2>0.6 -> eQTL(PRA) -> TFBS -> DNase -> CpG-Island -> [7 Genomic Regions](S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed), [7 SNPs](S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed)
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/TFBS_GWAS_RA_SNP
for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done

## Step 1. eQTL
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL
###############################
qval_threshold=0.05
eqtl<-c()
for(i in c("Whole_Blood","Whole_Blood","Liver","Small_Intestine_Terminal_Ileum","Stomach","Colon_Sigmoid","Lung","Spleen","Ovary")){
data<-subset(read.table(paste(i,".v7.egenes.txt",sep=""),head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(eqtl,as.character(data[,19]))
}
length(table(eqtl))
eqtl.snp<-names(table(eqtl))
write.table(eqtl.snp,file="eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
#################################
input<-read.table("../GWAS-RA-792.R2.6.rsSNP.hg19.bed",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)
#################################
## Step 2. TFBS-DNase-CpGIsland
cd /gpfs/home/guosa/hpc/rheumatology/RA/TFBS_GWAS_RA_SNP
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip wgEncodeRegTfbsClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3"\t"$4}' wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClusteredV3.hg19.bed
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz
gunzip wgEncodeRegDnaseClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3}' wgEncodeRegDnaseClusteredV3.bed > wgEncodeRegDnaseClusteredV3.hg19.bed

bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.eQTL.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed > GWAS-RA-R2.6.eQTL.tfbs.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.eQTL.tfbs.hg19.bed -b wgEncodeRegDnaseClusteredV3.bed | sort -u > GWAS-RA-R2.6.eQTL.tfbs.DNase.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.eQTL.tfbs.DNase.hg19.bed -b ~/hpc/db/hg19/CpGI.hg19.bed > GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed
bedtools sort -i GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed > GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed
bedtools merge -d 2000 -i GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed > S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed
```
S2: R2>0.6 -> eQTL(Full) -> TFBS -> DNase -> CpG-Island -> [17 Genomic Regions](S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed) and [19 SNPs](S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed)
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/TFBS_GWAS_RA_SNP
for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done
## Step 1. eQTL
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL
qval_threshold=0.05
###############################
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
input<-read.table("../GWAS-RA-792.R2.6.rsSNP.hg19.bed",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)
#################################
## Step 2. TFBS-DNase-CpGIsland
cd /gpfs/home/guosa/hpc/rheumatology/RA/TFBS_GWAS_RA_SNP
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip wgEncodeRegTfbsClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3"\t"$4}' wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClusteredV3.hg19.bed
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz
gunzip wgEncodeRegDnaseClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3}' wgEncodeRegDnaseClusteredV3.bed > wgEncodeRegDnaseClusteredV3.hg19.bed

bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.eQTL.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed > GWAS-RA-R2.6.eQTL.tfbs.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.eQTL.tfbs.hg19.bed -b wgEncodeRegDnaseClusteredV3.bed | sort -u > GWAS-RA-R2.6.eQTL.tfbs.DNase.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.eQTL.tfbs.DNase.hg19.bed -b ~/hpc/db/hg19/CpGI.hg19.bed > GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed
bedtools sort -i GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed > S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed
bedtools merge -d 2000 -i S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed > S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed
```
* S3: R2>0.6 -> mis-sense variants -> [104 SNPs](gnomad.exomes.r2.1.sites.rec.GWAS-RA-792.R2.6.rsSNP.input.hg19.vcf.bed)
```
cp GWAS-RA-792.R2.6.rsSNP.hg19.bed GWAS-RA-792.R2.6.rsSNP.input.hg19.bed
perl -p -i -e "s/chr//" GWAS-RA-792.R2.6.rsSNP.input.hg19.bed

panel="GWAS-RA-792.R2.6.rsSNP.input"
## Function Variants in Exom Regions
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \&INFO/AF_eas\<0.999 \& \(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"frame\"\)\' -R $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
```

