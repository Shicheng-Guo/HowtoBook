```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/All_Target_Locations.hg19.bed
for i in `ls *.bed`; 
do 
perl format.pl $i > $i.txt & 
done
```
