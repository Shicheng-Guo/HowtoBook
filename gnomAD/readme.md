gnomAD (hg19)

First, you need download gnomAD data form [GnomAD website](https://gnomad.broadinstitute.org/). In the website, we have two kind of data: 1) genome-wide 2) Exom-wide. They are for different usage. Exom-wide data is for Missense/frameshift/nonsense scanning and genome-wide data is for TFBS, miRNA-binding and so on. In the raw data, multi-allelic SNPs were splited to multiple line which is not the default setting for most of the software, therefore, the first thing after the download is to merge multi-allelic data to one-line in the vcf files. 
```
cd /gpfs/home/guosa/hpc/db/Gnomad/

for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
done
for i in {1..22} X 

do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo tabix -p vcf gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
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
7. prepare eQTL annotation
```
```
8. prepare methylation annotation
```
```


