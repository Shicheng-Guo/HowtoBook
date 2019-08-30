#### How to add contig lenght to VCF header downloaded by [Michigan Imputation Server](https://imputationserver.sph.umich.edu/start.html#!pages/home)

Recently, some of my colleagues asked help from me in which they mentioned the contig length information was lost in the imputation files downloaded form [Michigan Imputation Server, 09/01/2019](https://imputationserver.sph.umich.edu/start.html#!pages/home). I do checked it and found MIS do lost contig length information which usually required by lots of downstream analysis packages, such as GATK. 

Here is the solution I recommeneded. 
1. obtain the vcf header and prepare to replance the contig line
```
for i in {1..23}
do
bcftools view -h chr$i.dose.filter.vcf.gz >  chr$i.dose.filter.header
done  
```
2. prepare a perl script to insert contig length
```
```
