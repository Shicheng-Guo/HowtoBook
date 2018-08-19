

1. 1000 Genome phase 3 Data Download (2535 samples, 99 CEU+103CHB+108CHSls ) 
```
./sh download.sh
```
2. uncompress vcf.gz since I need to remove duplications with my own perl script.
```
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo gunzip ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> chr$i.job
qsub chr$i.job
done
```


