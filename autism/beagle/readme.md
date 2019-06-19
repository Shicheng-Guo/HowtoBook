

```
vcftools --vcf All_samples_Exome_QC.vcf --not-chr 0 --recode --out All_samples_Exome_QC.clean.vcf
bcftools norm -d both --threads=32 All_samples_Exome_QC.clean.vcf.recode.vcf -Ov  -o All_samples_Exome_QC.clean.norm.vcf
bcftools view -i 'ALT !="-"' All_samples_Exome_QC.clean.norm.vcf.gz -Oz -o All_samples_Exome_QC.norm.vcf.gz
tabix -p vcf All_samples_Exome_QC.norm.vcf.gz
sh remove_VCF_duplicates.sh All_samples_Exome_QC.norm.vcf.gz > All_samples_Exome_QC.clean.norm.undup.vcf
java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar impute=false gt=All_samples_Exome_QC.clean.norm.undup.vcf out=All_samples_Exome_QC.clean.norm.vcf.phasing
```
