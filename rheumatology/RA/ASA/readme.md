
Identify VIP gene bi-allelic SNPs to be genotyped
```
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \& INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R VIP.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
```
