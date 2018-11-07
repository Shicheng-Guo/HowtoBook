Step 3. Phase the genotyping data with beagle or michigen imputation server


Step 4. Do the annotation to identify functional variants


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
