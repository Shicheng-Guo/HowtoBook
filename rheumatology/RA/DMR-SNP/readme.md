Identify GWAS LD DMR regions in three population:
```
cat /gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB.txt  /gpfs/home/guosa/hpc/db/hg19/1000Genome/CHS.txt > /gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB_CHS_221.txt

# CHB+CHS
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150
awk '{cmd="plink --bfile ~/hpc/db/hg19/1000Genome/plink/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB_CHS_221.txt' --ld "$4" "$5 " --out './ld/'"$4"."$5".CHB.CHS.r | qsub -N "$4"."$5 " -e ./temp/ -o ./temp/";system(cmd)}'  WGBS.bed.cSNP150.bed.sort.bed.pair.bed

# CEU Sample
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150
awk '{cmd="plink --bfile ~/hpc/db/hg19/1000Genome/plink/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt' --ld "$4" "$5 " --out './ld/'"$4"."$5".CEU.r | qsub -N "$4"."$5 " -e ./temp/ -o ./temp/";system(cmd)}'  WGBS.bed.cSNP150.bed.sort.bed.pair.bed

# YRI Sample
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150
awk '{cmd="plink --bfile ~/hpc/db/hg19/1000Genome/plink/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt' --ld "$4" "$5 " --out './ld/'"$4"."$5".YRI.r | qsub -N "$4"."$5 " -e ./temp/ -o ./temp/";system(cmd)}'  WGBS.bed.cSNP150.bed.sort.bed.pair.bed

grep R-sq *log | awk -F'[=\sD]' '$5>0.1{print}'
 ```
