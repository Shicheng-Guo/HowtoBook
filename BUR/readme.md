
```
cd /gpfs/home/guosa/run/bedmethyl
for i in chr{7..22} chrX chrY chr6 chr5 chr4 chr3 chr2 chr1 
do
for j in `ls *bw`
do
echo \#PBS -N $j.$i  > $j.$i.job
echo \#PBS -l nodes=1:ppn=1 >> $j.$i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $j.$i.job
echo bigWigAverageOverBed $j /gpfs/home/guosa/hpc/db/hg38/window200/hg38.$i.win2K.bed $j.$i.tab >> $j.$i.job
qsub $j.$i.job
done
done
```
