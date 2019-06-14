mkdir EUR
mkdir temp
for i in {1..22} X Y
do
echo #PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf chr$i.1kg.phase3.v5a.vcf.gz >> $i.job
echo bcftools view chr$i.1kg.phase3.v5a.vcf.gz -S EUR.SampleList.txt -Oz -o ./EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz >>$i.job
qsub $i.job
done
