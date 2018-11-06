

* ERAP1

rs27044,rs17482078,rs10050860,rs30187,rs2287987

```
#HLA-DRB1*0401
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs660895 --snp rs6910071 --snp rs3817964 --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out HLA-DRB1-0401 --hap-r2  --geno-r2  --geno-chisq --freq  --counts 
#HLA-B27
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs13202464 --snp rs4349859 --snp rs3819299 --snp rs116488202 --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out HLA-B27 --hap-r2  --geno-r2  --geno-chisq --freq  --counts 
```
