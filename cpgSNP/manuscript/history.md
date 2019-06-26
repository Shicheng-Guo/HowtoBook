
```
cd /gpfs/home/guosa/hpc/db/hg19/cpgSNP/commonSNP

cpgSNPisland.hg19.bed
cpgSNPisland.hg19.SNP.bed
cpgSNPisland.hg19.binarySNP.bed
```
In order to obtain the summary for different GAP to LOSS or GAIN island. (contious K loss or gain in common CpG-SNPs)
```
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/beagle/cgSNP")
setwd("/home/guosa/hpc/db/hg19/beagle/cgSNP")
head(data)
for(GAP in seq(5,20,1)){
rlt<-c()
for(j in c(1:22,"X","Y")){
data<-read.table(paste("chr",j,".cpgSNP.bin.bed",sep=""))
i=1
print(paste("chr",j,".cpgSNP.bin.bed",sep=""))
while(i<(nrow(data)-100)){
  if(table(data[i:(i+GAP),6])[1]>GAP){
   rlt<-rbind(rlt,c(data[i,1],data[i,2],data[i+GAP,3]))
   print(i)
   i=i+GAP
  }
  i=i+1
}
}
write.table(rlt,file=paste("Long_CpG_Gain.GAP",GAP,"hg19.txt",sep="."),sep="\t",col.names = F,row.names = F)
}
````
extract 1000G VCF by CpG-SNP with `bcftools view -T`

```
cd /gpfs/home/guosa/hpc/db/hg19/beagle/cgSNP/vcf
mkdir temp
for i in {1..23}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mkdir ./temp/chr$i >> $i.job
echo bbcftools view -T ../chr$i.cpgSNP.bin.bed ../../chr$i.1kg.phase3.v5a.vcf.gz -Oz -o chr$i.cpgSNP.vcf.gz >>$i.job
qsub $i.job
done
```
