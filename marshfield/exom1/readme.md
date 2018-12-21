
```
cd ~/gpfs/home/guosa/hpc/project/pmrp/phase1/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'R2\>0.8\' $i.dose.vcf.gz -Oz -o $i.dose.filter.vcf.gz >> $i.job
echo tabix -p vcf $i.dose.filter.vcf.gz >> $i.job
echo zcat $i.dose.filter.vcf.gz \| awk \'{print \$1,\$2,\$2,\$4,\$5}\' OFS=\"\\t\" \| grep -v \'#\' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl $i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ../annovar/$i -remove -protocol refGene,dbnsfp35c,gwasCatalog -operation gx,f,r -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub $i.job
done
```

