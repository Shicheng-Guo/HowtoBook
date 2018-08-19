

1. 1000 Genome phase 3 Data Download (2535 samples, 99 CEU+103CHB+108CHSls ) 
```
./sh download.sh
```
2. uncompress vcf.gz since I need to remove duplications with my own perl script.
```
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo gunzip ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> chr$i.job
qsub chr$i.job
done
```
3. remove duplication
```
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcf2dedup ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \> chr$i.uni.vcf >> chr$i.job
qsub chr$i.job
done
```
4. achieve sub-dataset only conditioning DMER SNPS
```
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
awk '{print $4}' /gpfs/home/guosa/hpc/db/hg19/commonsnp150.hg19.bed > commonsnp150.rs.list.tmp
sort -u   commonsnp150.rs.list.tmp >  commonsnp150.rs.list
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --recode --recode-INFO-all --snps dmer.snp.list --out chr$i.dmer >> chr$i.job
qsub chr$i.job
done
```
5. split the sample by the popultion
```
awk '{print $1,$2 > ""$3".txt"}' 1000GenomeSampleInfo.txt
```
6. calcluate LD between dmer snp and gwas proxy snp
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X Y
do
echo \#PBS -N $i.CEU  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=~/hpc/db/hg19/1000Genome/CEU.txt
vcf=~/hpc/db/hg19/1000Genome/chr$j.uni.vcf
echo plink --vcf $vcf --keep $ceu --r2 --ld-snp-list $i.pairSNP --out ./CEU/$i.chr$j.ld >> $i.$j.job
qsub $i.$j.job
done
done
```
