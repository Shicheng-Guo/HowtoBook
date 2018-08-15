for i in chr{1..22} chrX chrY
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedgraph >> $i.job
echo Rscript --vanilla BUR.R $i.tab.matrix.rlt >> $i.job
qsub $i.job
done
