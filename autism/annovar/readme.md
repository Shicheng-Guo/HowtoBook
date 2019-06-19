Step 1. annotation SNPs with annovar (dbnsfp35c,hg19)
```
mkdir annovar
mkdir temp
for i in {1..22} 
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> chr$i.job
echo \#PBS -m abe  >> chr$i.job
echo \#PBS -o $(pwd)/temp/ >>chr$i.job
echo \#PBS -e $(pwd)/temp/ >>chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq All_samples_Exome_QC.chr$i.vcf.vcf.gz  \> ./annovar/All_samples_Exome_QC.chr$i.avinput >> chr$i.job
echo table_annovar.pl ./annovar/All_samples_Exome_QC.chr$i.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ./annovar/chr$i -remove -protocol ensGene,dbnsfp35c,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,f,r,r,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done
```
Step 2. filter delterious SNPs with dbnsfp35c (6 out of 11 as deleterious)
```

```
