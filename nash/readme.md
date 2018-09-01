
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
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/ >> ${i}.job
echo bismark_methylation_extractor ${option1} ${option2} --output ../methyfreq  ./$i >> ${i}.job
echo ${i}.job
done
```