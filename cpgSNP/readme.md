```
for i in `ls *cpgsnp.bed`
do
awk '{print $4}' $i >> cpgSNP.hg19.list.txt
done
```
```
cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr$i --extract /home/guosa/hpc/db/hg19/plan2/commonCpGSNP/cpgSNP.hg19.list.txt --make-bed --out chr$i.cpgSNP
done
```
