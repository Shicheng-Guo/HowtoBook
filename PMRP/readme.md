1. CEU Reference from 1000 Genome Project
```
cd /gpfs/home/guosa/hpc/tools/beagle/
for i in {1..24}
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo plink --vcf chr$i.1kg.phase3.v5a.vcf.gz --recode vcf --keep EUR.incl --make-bed --out chr$i.1kg.phase3.v5a >> chr$i.job
qsub chr$i.job
done
```


```
cd /gpfs/home/guosa/hpc/project/pmrp/phase2/phase
for i in {1..24}
do
plink --bfile S_Hebbring_Unr.Guo --recode vcf --chr $i --snps-only just-acgt --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i
done
```

```
pdf("MRCI.PMRP.MAF.distribution.pdf")
hist(log(newdata$MAF,10),xlim=c(-5,0),main="Histogram of Log(MAF, 10)",xlab="Log(MAF, 10)",col="blue")
dev.off()
```
