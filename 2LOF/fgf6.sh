
```
bcftools view -i 'R2>0.5' chr12.dose.vcf.gz -Oz -o chr12.dose.filter.vcf.gz
zcat chr12.dose.filter.vcf.gz | awk '{print $1,$2,$3,$4,$5}' OFS="\t" | grep -v '#' > chr12.vcf.avinput
table_annovar.pl chr12.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out chr12 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt
```
