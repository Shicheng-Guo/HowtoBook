Step 1. Single variants based rvtest
```
for i in {1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo rvtest --inVcf All_samples_Exome_QC.chr$i.vcf.vcf.gz --pheno All_samples_Exome_QC.phen --single wald,score --out All_samples_Exome_QC.chr$i.rvtest >>$i.job
qsub $i.job
done
```

Step 2. Gene based rvtest 
```
wget http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat_hg19.txt.gz
gunzip refFlat_hg19.txt.gz
perl -p -i -e 's/chr//' refFlat_hg19.txt
gzip refFlat_hg19.txt

for i in {1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf All_samples_Exome_QC.chr$i.vcf.vcf.gz >>$i.job
echo rvtest --inVcf All_samples_Exome_QC.chr$i.vcf.vcf.gz --pheno All_samples_Exome_QC.phen --out All_samples_Exome_QC.chr$i.rvtest.skat --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac >>$i.job
qsub $i.job
done
````
