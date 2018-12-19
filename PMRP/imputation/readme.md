
```
for i in {6..22}
do
7za x chr_$i.zip -ppydG2EucK5bIhG &
done

cd /gpfs/home/guosa/hpc/project/pmrp/phase1/imputation
for i in {6..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo 7za x chr_$i.zip -ppydG2EucK5bIhG >>chr$i.job
qsub chr$i.job
done
```

