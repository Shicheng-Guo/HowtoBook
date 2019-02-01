#### bcftools — utilities for variant calling and manipulating VCFs and BCFs.
#### version: 1.9- 2019-01-31

How to select bi-allelic SNPs from VCF files. (-d and -m cannot be used at the same time)
```
bcftools sort check.temp.17.vcf > check.17.vcf
bcftools norm  -m + check.17.vcf > check.17.recover.vcf
bcftools norm -d both check.17.recover.vcf > check.17.recover.2.vcf
bcftools view -m2 -M2 -v snps check.17.recover.2.vcf > check.17.trim.vcf
grep rs1799966 check.17.trim.vcf | less -S 
```


