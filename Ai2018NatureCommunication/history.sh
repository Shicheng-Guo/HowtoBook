###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo bcftools norm -d both --threads=2 ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O z -o chr$i.vcf.gz >> chr$i.job
qsub chr$i.job
done
###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --exclude mydup.txt --out chr$i.unique >> chr$i.job
qsub chr$i.job
done
###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
awk '{print $4}' /gpfs/home/guosa/hpc/db/hg19/commonsnp150.hg19.bed > commonsnp150.rs.list.tmp
sort -u   commonsnp150.rs.list.tmp >  commonsnp150.rs.list
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --recode --recode-INFO-all --snps commonsnp150.rs.list --out chr$i.common150 >> chr$i.job
qsub chr$i.job
done
###########################
for i in `ls *bed.sort.bed `
do
echo \#PBS -N $i  > $i.job
echo cd $(pwd) >>$i.job
echo bedtools intersect -wao -a $i -b  ~/hpc/db/hg19/commonsnp150.hg19.bed \> $i.SNP150 >> $i.job
qsub $i.job
done
###########################
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in `ls *SNP150`
do
bedtools closest -a $i -b RA_GWAS_475_Catalog_GRCH37.bed | awk '($10!="." && $11!="."){print $1"\t"$2"\t"$3"\t"$10"\t"$18}' > $i.pair.bed  &
echo $i
done
###########################
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
awk '{print $4"\t"$5}' $i.bed.sort.bed.SNP150.pair.bed | sort -u > $i.pairSNP &
done
