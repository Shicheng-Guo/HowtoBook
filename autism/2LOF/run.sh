cd /home/guosa/hpc/autism/data/2LOF
mkdir 2LOF
cp All* 2LOF
plink --bfile All_samples_Exome_QC --recode vcf --out All_samples_Exome_QC
 
