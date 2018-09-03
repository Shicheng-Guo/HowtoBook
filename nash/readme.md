

```
cd /home/guosa/hpc/db/hg19
mkdir window2000/
for i in {1..22} X Y M
do
perl ~/hpc/bin/cutchrosome.pl chr$i 500 >  ./window2000/hg19.chr$i.win2K.bed
done

```

```
cp /gpfs/home/guosa/hpc/nash/bam/pool/methyfreq/*bismark.cov.gz /gpfs/home/guosa/hpc/nash/methcov
cp /gpfs/home/guosa/hpc/nash/methyfreq/*bismark.cov.gz  /gpfs/home/guosa/hpc/nash/methcov

```
```
U16_S30_L001_R1_001_00_bismark_bt2_pe.bam
V03_S14_L001_R1_001_00_bismark_bt2_pe.bam
W01_S43_L001_R1_001_00_bismark_bt2_pe.bam
W06_S24_L001_R1_001_00_bismark_bt2_pe.bam
```

```
mkdir ../methyfreq
option1=$(echo --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end)
option2=$(echo --bedGraph --ignore 1 --buffer_size 4G --gzip --comprehensive)
for i in `ls *bam`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/ >> ${i}.job
echo bismark_methylation_extractor ${option1} ${option2} --output ../methyfreq  ./$i >> ${i}.job
echo ${i}.job
done
```
```
qsub U16_S30_L001_R1_001_00_bismark_bt2_pe.bam.job
qsub V03_S14_L001_R1_001_00_bismark_bt2_pe.bam.job
qsub W01_S43_L001_R1_001_00_bismark_bt2_pe.bam.job
qsub W06_S24_L001_R1_001_00_bismark_bt2_pe.bam.job
```

```
mkdir ../methyfreq
option1=$(echo --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end)
option2=$(echo --bedGraph --ignore 1 --buffer_size 4G --gzip --comprehensive)
for i in `ls Pool_*bam`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/pool/bam >> ${i}.job
echo bismark_methylation_extractor ${option1} ${option2} --output ../methyfreq  ./$i >> ${i}.job
echo ${i}.job
done
```
```
cp U16_S30*.cov.gz ../methcov/
cp V03_S14*.cov.gz ../methcov/
cp W01_S43*.cov.gz ../methcov/
cp W06_S24*.cov.gz ../methcov/
```
cov2bedgraph
```
cd /gpfs/home/guosa/hpc/nash/methcov
for i in `ls *.cov`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo awk \'{print \$1,\$2-1,\$3,\$4}\' $i OFS=\"\\t\" \>$i.bedgraph >>$i.job
echo wigToBigWig $i.bedgraph ~/hpc/db/hg19/hg19.chrom.sizes $i.bedgraph.bw >> $i.job
echo $i.job
qsub $i.job
done
```

tab2matrix
```
```

```
cd /gpfs/home/guosa/hpc/nash/methcov
for i in `ls *bw`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo bigWigAverageOverBed $i ~/hpc/db/hg19/LINE1.hg19.bed ./LINE-1/$i.tab >> $i.job
echo $i.job
qsub $i.job
done
```

