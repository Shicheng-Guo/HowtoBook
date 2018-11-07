Step 3. Phase the genotyping data with beagle or michigen imputation server


Step 4. Do the annotation to identify functional variants
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



Step 5. Run the test
```
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF
for i in `ls *.update.vcf`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R $i >> $i.job
qsub $i.job
done
```
