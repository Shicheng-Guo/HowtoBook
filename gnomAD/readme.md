1. download gnomAD
```
for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done
```
2. extract the VEP annotation to SNV

```
cd /gpfs/home/guosa/hpc/db/Gnomad
for i in `ls *bgz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo cd $(pwd) >> $i.job
echo mv $i > $i.gz >> $i.job
echo tabix -p vcf $i.gz >> $i.job
qsub $i.job
done
```
