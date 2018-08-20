

1. 1000 Genome phase 3 Data Download (2535 samples, 99 CEU+103CHB+YRI ) 
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
6. extract the sub-SNPs to decrease the analysis time.
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
7. vcf to vcf.gz and tabix index for vcf.gz files
```
cd /gpfs/home/guosa/nc
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo bgzip -c chr$i.uni.vcf \> chr$i.uni.vcf.gz >> chr$i.job
echo tabix -p vcf chr$i.uni.vcf.gz >> chr$i.job
qsub chr$i.job
done
```
8.  perl script to calculate LD d' pair-by-pair
```
cd /gpfs/home/guosa/nc/
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do
for j in CEU CHB YRI 
do
perl ldcal.pbs.pl $i $j
done
done

```
9. summary all the LD d' and assign it to bed files
```
cd /gpfs/home/guosa/nc/
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in CEU CHB YRI 
do
perl ldsum.pl $i $j
done
done
```

10. genomic shuffling and do the statistic
```
for j in {1..1000}
do
for i in `ls *bed.sort.bed `
do
echo \#PBS -N $i.$j  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo \#PBS -o ./temp/ >>$i.$j.job
echo \#PBS -e ./temp/>> $i.$j.job
echo bedtools shuffle -i $i -excl wgEncodeDukeMapabilityRegionsExcludable.bed -g ~/hpc/db/hg19/hg19.chrom.sizes \> ./Shuffle/$i.$j.shuffle >> $i.$j.job
echo bedtools sort -i ./Shuffle/$i.$j.shuffle  \> ./Shuffle/$i.$j.shuffle.sort  >> $i.$j.job 
echo bedtools closest -a ./Shuffle/$i.$j.shuffle.sort -b RA_GWAS_475_Catalog_GRCH37.bed \> ./Shuffle/$i.$j.GWAS.Cloest >> $i.$j.job
qsub $i.$j.job
done
done
```


######Supplementary Figure and Method
unknown. calcluate LD between dmer snp and gwas proxy snp
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for k in CEU CHB CHS ACB ASW BEB CDX CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI 
do
mkdir $k
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X Y
do
echo \#PBS -N $i.$k  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
pop=~/hpc/db/hg19/1000Genome/$k.txt
vcf=~/hpc/db/hg19/1000Genome/chr$j.uni.vcf
echo plink --vcf $vcf --keep $pop --r2 --ld-snp-list $i.pairSNP --out ./$k/$i.chr$j.ld >> $i.$j.job
qsub $i.$j.job
done
done
done

```
