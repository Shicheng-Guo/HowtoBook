

```
vcftools --vcf All_samples_Exome_QC.vcf --not-chr 0 --recode --out All_samples_Exome_QC.clean.vcf
bcftools norm -d both --threads=32 All_samples_Exome_QC.clean.vcf.recode.vcf -Ov  -o All_samples_Exome_QC.clean.norm.vcf
java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar impute=false gt=All_samples_Exome_QC.clean.norm.vcf out=All_samples_Exome_QC.clean.norm.vcf.phasing
```
