Step 1. aloft can only be conducted in hpc
```
cd ~/hpc/autism/data/aloft
awk '$4 !="-" && $5 !="-" {print}' OFS="\t" All_samples_Exome_QC.DG5.vcf > All_samples_Exome_QC.DG5.temp.vcf
vcftools --vcf All_samples_Exome_QC.DG5.temp.vcf --not-chr 0 --recode --out All_samples_Exome_QC.DG5.clean.vcf
aloft --vcf All_samples_Exome_QC.DG5.clean.vcf --output All_samples_Exome_QC.DG --data ~/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/
```
Step 2. filter out LOF event and variants 
```
```
