```
cd /home/guosa/hpc/project/pmrp/phase1/imputation
ls chr*dose.filter.vcf.gz > concat.txt
bcftools concat -f concat.txt -O z  -o exom1.imputate.filter.vcf.gz

cd ~/hpc/project/pmrp/phase1/imputation
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --vcf chr$i.dose.filter.vcf.gz --make-bed --out Exom1.chr$i.dose.filter.vcf.gz  >>  $i.job
qsub $i.job
done


``` 
