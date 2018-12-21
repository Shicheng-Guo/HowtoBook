gnomAD (hg19)

1. download gnomAD
```
for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done
```
2. bgz to gz and re-tabix
```
cd /gpfs/home/guosa/hpc/db/Gnomad
for i in `ls *bgz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo cd $(pwd) >> $i.job
echo mv $i \> $i.gz >> $i.job
echo tabix -p vcf $i.gz >> $i.job
qsub $i.job
done
```
3.  extract the VEP annotation to SNV
```
cd /gpfs/home/guosa/hpc/db/Gnomad
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo perl recode.pl gnomad.genomes.r2.1.sites.chr$i.vcf.bgz \> gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.annovar.txt >> $i.job
qsub $i.job
done
```
4. remove non-rs SNPs
```
for i in {1..22}
do
perl -lane '{print $_ if @F[6]=~/rs/}' gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.annovar.txt > gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.annovar.hg19.bed
done
```
5. prepre gnomad annotation to annovar db (4 column bed file is the best choice for annovar-anno-db)
```
cd /gpfs/home/guosa/hpc/db/Gnomad
cat *.txt > gnomad.genomes.r2.1.sites.chr2.vcf.bgz.annovar.db
perl index.pl gnomad.genomes.r2.1.sites.chr2.vcf.bgz.annovar.db 1000
```
6. annotate SNPs with gnomad-annovar-db
```

```



