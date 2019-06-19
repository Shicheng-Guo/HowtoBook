Step 1. aloft can only be conducted in hpc
```
plink --bfile All_samples_Exome_QC --recode vcf  --out ./aloft/All_samples_Exome_QC
cd ~/hpc/autism/data/aloft
bcftools view -G All_samples_Exome_QC.vcf -Ov -o All_samples_Exome_QC.DG5.vcf
vcftools --vcf All_samples_Exome_QC.DG5.vcf --not-chr 0 --recode --out All_samples_Exome_QC.DG5.clean.vcf
aloft --vcf All_samples_Exome_QC.DG5.clean.vcf --output All_samples_Exome_QC.DG --data ~/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/
```
Step 2. filter out LOF event and variants 
```
```
