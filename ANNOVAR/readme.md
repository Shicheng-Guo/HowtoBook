How to use ANNOVAR to annotate VCF files:

```
## Step 1:  Download Annotation Files
cd /home/guosa/hpc/tools/annovar
# just for allele frequency
annotate_variation.pl -downdb -webfrom annovar exac03 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2_all humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar gnomad_exome humandb -buildver hg19 &
# whole-exome data
annotate_variation.pl -downdb -webfrom annovar 1000g2015aug humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar kaviar_20150923 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar hrcr1 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cg69 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar gnomad_genome humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar dbnsfp30a humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg19 &
annotate_variation.pl -downdb esp6500siv2 humandb -buildver hg19 &
#  whole-genome data
annotate_variation.pl -downdb -webfrom annovar gerp++ humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cadd humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cadd13 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar fathmm humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar eigen humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar gwava humandb -buildver hg19  &
# for CNV
annotate_variation.pl -downdb -webfrom annovar dbscsnv11 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar spidex humandb -buildver hg19  &
# disease-specific variants
annotate_variation.pl -downdb -webfrom annovar clinvar_20160302 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cosmic70 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar icgc21 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar nci60 humandb -buildver hg19  &
# miSNP
# CpG-SNP
# eQTL

## Step 2:  Creat annovar input from VCF files
echo convert2annovar.pl -format vcf4 chr$i.dose.filter.vcf.gz  \> ../annovar/chr$i.dose.vcf.avinput >> chr$i.job

## Step 3:  Annotation input vcf 
table_annovar.pl ex1.avinput humandb/ -protocol dbnsfp30a -operation f -build hg19 -nastring .

# PBS code
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo bcftools view -i \'R2\>0.9\' chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.vcf.gz>> chr$i.job
echo convert2annovar.pl -format vcf4 chr$i.dose.filter.vcf.gz  \> ../annovar/chr$i.dose.vcf.avinput >> chr$i.job
echo annotate_variation.pl -out ../annovar/chr$i.annovar.txt -build hg19 ../annovar/chr$i.dose.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ >> chr$i.job
qsub chr$i.job
done
```


# Annotation and Compound heterozygotes scanning
'''
cd ~/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo -withzyg chr$i.dose.filter.9.vcf.gz  \> ../annovar/chr$i.dose.vcf.9.avinput >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.dose.vcf.9.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr$i.9 -remove -protocol refGene,ljb23_fathmm,ljb23_metasvm,ljb23_metalr,eigen,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done
'''
