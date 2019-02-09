cd /gpfs/home/guosa/hpc/db/hg19/1000Genome

mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> $i.job
qsub $i.job
done

