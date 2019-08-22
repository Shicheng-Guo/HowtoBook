```
cd /gpfs/home/guosa/run/bedmethyl


for i in `ls *.bed`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $i.job
echo perl bedMethyl2bedgraph.pl $i >> $i.job
echo wigToBigWig $i.bedgraph ~/hpc/db/hg38/hg38.chrom.sizes $i.bedgraph.bw >> $i.job
qsub $i.job
done

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

```
for i in chr{1..22} chrX chrY chrM
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org
echo \#PBS -m abe  > $i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $i.job
echo perl ~/hpc/bin/tab2matrix.pl $i \> $i.tab.matrix.rlt >> $i.job
echo Rscript --vanilla dmr.R $i.tab.matrix.rlt >> $i.job
qsub $i.job
done
```
