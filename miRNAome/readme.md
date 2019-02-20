364291 uniuqe UTR-3-miRNA-SNPs and [44952](gnomad.genomes.r2.1.sites.rec.UTR3miRNAsNP.merge.vcf.bed) common SNPs (MAF>1%) in East Asian. 


```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/All_Target_Locations.hg19.bed
for i in `ls *.bed`; 
do 
perl format.pl $i > $i.txt & 
done
```
