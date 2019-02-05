
```
69   PADI(1-4) Variants in RA Risk
1763 Epigene Functional Variants
3325 GWAS Immnue System Desases
16013 CpGSNP-eQTL-Overlap Variants
```

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

Identify EpiFactors gene bi-allelic SNPs to be genotyped
```
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R Epigene.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.Epi.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.Epi.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.Epi.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf.bed
```

