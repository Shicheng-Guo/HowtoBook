cd /home/guosa/hpc/autism/data/2LOF
mkdir 2LOF
cp All* 2LOF
plink --bfile All_samples_Exome_QC --recode vcf --out All_samples_Exome_QC
avinput="All_samples_Exome_QC.vcf"
table_annovar.pl $avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out $avinput -remove -protocol refGene,dbnsfp35c,gnomad_genome,gwasCatalog,gtexEqtlCluster -operation g,f,f,r,r -nastring . 



