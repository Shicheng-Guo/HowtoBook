mkdir temp
for i in {1..22} 
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=8 >> chr$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> chr$i.job
echo \#PBS -m abe  >> chr$i.job
echo \#PBS -o $(pwd)/temp/ >>chr$i.job
echo \#PBS -e $(pwd)/temp/ >>chr$i.job
echo cd $(pwd) >> chr$i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >>chr$i.job
qsub chr$i.job
done
