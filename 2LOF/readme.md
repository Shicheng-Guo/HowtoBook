Step 3. Phase the genotyping data with beagle or michigen imputation server
```
cd ~/hpc/project/pmrp/Exom2/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'R2\>0.9\' $i.dose.vcf.gz -Oz -o $i.dose.filter.9.vcf.gz >> $i.job
echo tabix -p vcf $i.dose.trim.9.vcf.gz >> $i.job
echo zcat $i.dose.trim.9.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl ../annovar/$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ../annovar/$i -remove -protocol refGene,dbnsfp30c,gwasCatalog -operation gx,f,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub $i.job
done
```


Step 3. Do the annotation to identify functional variants
```
cd ~/hpc/project/pmrp/Exom2/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'R2\>0.9\' $i.dose.vcf.gz -Oz -o $i.dose.filter.9.vcf.gz >> $i.job
echo tabix -p vcf $i.dose.trim.9.vcf.gz >> $i.job
echo zcat $i.dose.trim.9.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl ../annovar/$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ../annovar/$i -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub $i.job
done
```
Step 4. rebuild vcf files only containning functional loss-of-function variance in the samples
```
# add INFO annotation to LOSS ALLLELS
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla readanno.R $i >> $i.job
qsub $i.job
done
```

step 5. build annotation file for vcf as the input of bcftools
```
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo sort -k1,1 -k2,2n $i.anno.txt \> $i.sort.anno.txt >> $i.job
echo bgzip -c $i.sort.anno.txt \> $i.anno.gz  >> $i.job
echo tabix -s1 -b2 -e2 $i.anno.gz >> $i.job
qsub $i.job
done
```
step 6. add INFO annotation to LOSS ALLLELS
```
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF/
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view ../imputation/$i.dose.trim.9.vcf.gz -R ../annovar/$i.anno.txt -Oz -o ../2LOF/$i.vcf.gz >> $i.job
echo bcftools annotate -x INFO,FILTER ../2LOF/$i.vcf.gz -Oz -o ../2LOF/$i.vcf.tmp.gz >> $i.job
echo bcftools annotate -a ../annovar/$i.anno.gz -h ../annovar/head.hdr -c CHROM,POS,REF,ALT,-,-,TAG3 $i.vcf.tmp.gz -Oz -o $i.update.vcf.gz >> $i.job
echo rm ../2LOF/$i.vcf.tmp.gz >> $i.job
qsub $i.job
done
```
Step 7a. Run the test, each phenotype, each run
```

cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF
for i in {1..22}
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp10_Obesity_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp1_RA_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp7_Iron_C1_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp7_Iron_C2_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp3_PA_rev2_SampleIDs.Michigen.txt >> $i.job 
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp5_Thyroid_C1_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp5_Thyroid_C2_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_SSc_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ANA_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ENA_rev2_SampleIDs.Michigen.txt >> $i.job
qsub $i.job
done
```
Step 7b. Speed up analysis, each chr, each phenotype, each run 
```
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF
for i in 13 19 
do
for j in `ls /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phe*.Michigen.txt`
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R chr$i.update.vcf $j >> $i.job
qsub $i.job
done
done
```
