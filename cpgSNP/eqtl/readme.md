### CpG-SNP-Loss-Gain and eQTL in human tissues

* `gene_fullxref.txt` avoid + : [] characteristic in this file.  

```
cd ~/hpc/db/hg19
perl annovar2bed.pl > cpgSNP.hg19.bed.V12.avinput
table_annovar.pl hg19_cpgSNP_eQTL.V6.txt ~/hpc/tools/annovar/humandb/ --thread 24 -buildver hg19 --csvout -out cpgSNP.hg19.bed.V12.avinput -remove -protocol knownGene,gwasCatalog,tfbsConsSites -operation gx,r,r -nastring . -otherinfo  -polish  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 
 ```



