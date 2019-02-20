364291 uniuqe UTR-3-miRNA-SNPs and [44952](gnomad.genomes.r2.1.sites.rec.UTR3miRNAsNP.merge.vcf.bed) common SNPs (MAF>1%) in East Asian. [4678](4678.UTR3miRNAsNP.EAS.MAF0.01.hg19.bed) SNPs were collected binded by miRNA in UTR3. 

```
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/All_Target_Locations.hg19.bed
for i in `ls *.bed`; 
do 
perl format.pl $i > $i.txt & 
done
```
