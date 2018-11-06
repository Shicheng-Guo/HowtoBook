HLA

* ERAP1

```
#HLA-DRB1*0401
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs660895 --snp rs6910071 --snp rs3817964 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out HLA-DRB1-0401 --hap-r2 --recode --geno-r2  --geno-chisq --freq  --counts 
#HLA-B27
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs13202464 --snp rs4349859 --snp rs3819299 --snp rs116488202 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out HLA-B27 --hap-r2  --geno-r2  --geno-chisq --freq  --counts 
#ERAP1
vcftools --gzvcf ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs27044 --snp rs17482078 --snp rs10050860 --snp rs30187 --snp rs2287987 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out ERAP1 --hap-r2  --geno-r2  --geno-chisq --freq  --counts 
```
