

In the first step, we collected all the eQTL from whole blood, liver, colon, stomach [etql set1](eqtl.set1.sh). 
* [21098](eQTL.hg19.bed) eQTL SNPs were collected from the GTEx project with FDR<0.05
* [14285](gnomad.genomes.r2.1.sites.rec.eQTL.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1%
* [772](gnomad.exomes.r2.1.sites.rec.eQTL.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1%

However, when I check it again, I found I lost the lung eqtl data. After I add lung eqtl [eqtl set2](eqtl.set2.sh)
* [30517](eQTL.hg19.bed) eQTL SNPs were collected from the GTEx project with FDR<0.05
* [14285](gnomad.genomes.r2.1.sites.rec.eQTL.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1%
* [772](gnomad.exomes.r2.1.sites.rec.eQTL.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1%

```
#####################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP
perl -p -i -e 's/chr//g' eQTL.hg19.bed

panel="eQTL"
bed="innateDbUTR3.hg19.bed"
mkdir temp
## Function Variants in Genome Regions
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \# bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \# tabix -p vcf gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \&INFO/AF_eas\<0.999\' -T $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *genomes*rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
wc -l gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
#####################################################################
```
