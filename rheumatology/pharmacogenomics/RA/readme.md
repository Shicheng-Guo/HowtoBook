### PharmacoGenomics of Rheumatoid Arthrits
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

# identify all the SNPs for 66 VIP genes

for i in {1..22} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view  -m2 -M2 -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \& INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R VIP.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.vcf.bgz \> VIP.chr$i.vcf >>$i.job
qsub $i.job
done

# merge all the VIP SNPs
perl -lane 'print $_ if $_ !~/#/' VIP*vcf | awk '{print $1,$2,$3,$4,$5}' > VIP.snp.hg19.txt

```
