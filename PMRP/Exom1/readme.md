

```
table_annovar.pl FinalRelease_QC_20140311_Team1_Marshfield.avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out FinalRelease_QC_20140311_Team1_Marshfield.avinput -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 
```
