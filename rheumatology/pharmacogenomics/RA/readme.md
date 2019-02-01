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
