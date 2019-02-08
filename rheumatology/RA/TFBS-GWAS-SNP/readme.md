DNA Methylation in GWAS Significant Rheumatoid Arthritis Associated Regions. 

We download 791 GWAS-Significant SNPs from GWAS Catalog and we collected all the linked SNPs with R2>0.6 in Asian Population. Totally, we identified 21,079 SNPs with above method. After TFBS (18621), DNase (2054) and CpG-island (129) filtering, we obtained the final target: 129 SNPs. We merged the adjust SNPs and found 77 genomic regions (see: GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.merge.hg19.bed)


```
cd /gpfs/home/guosa/hpc/rheumatology/RA/TFBS_GWAS_RA_SNP


for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip wgEncodeRegTfbsClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3"\t"$4}' wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClusteredV3.hg19.bed

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz
gunzip wgEncodeRegDnaseClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3}' wgEncodeRegDnaseClusteredV3.bed > wgEncodeRegDnaseClusteredV3.hg19.bed

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/TFBS-GWAS-SNP/GWAS-RA-792.hg19.bed

bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed > GWAS-RA-R2.6.tfbs.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.tfbs.hg19.bed -b wgEncodeRegDnaseClusteredV3.bed | sort -u > GWAS-RA-R2.6.tfbs.DNase.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.tfbs.DNase.hg19.bed -b ~/hpc/db/hg19/CpGI.hg19.bed > GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.bed
bedtools sort -i GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.bed > GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.bed
bedtools merge -d 2000 -i GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.bed > GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.merge.hg19.bed
```


