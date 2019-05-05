```
cd /home/guosa/hpc/project/pmrp/phase1/imputation
ls chr*dose.filter.vcf.gz > concat.txt
bcftools concat -f concat.txt -O z  -o exom1.imputate.filter.vcf.gz
``` 
