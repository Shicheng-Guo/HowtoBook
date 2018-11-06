### IBD to check ID

```
cd ~/hpc/project/pmrp/Exom2/IBD
perl -lane '{next if /^#/;print "@F[0]\t@F[1]\t@F[3]\t@F[4]"}' chr22.imputation.vcf > chr22.imputation.map
```
