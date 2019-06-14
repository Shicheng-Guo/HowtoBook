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
echo bcftools view -i \'DR2\>0.6\' All_samples_Exome_QC.chr$i.vcf.vcf.gz -Ov -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf >>$i.job
qsub $i.job
done
