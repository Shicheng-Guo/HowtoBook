#CEU
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X Y
do
echo \#PBS -N $i.CEU  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=~/hpc/db/hg19/1000Genome/CEU.txt
vcf=~/hpc/db/hg19/1000Genome/chr$j.uni.vcf
echo plink --vcf $vcf --keep $ceu --r2 --ld-snp-list $i.pairSNP --out ./CEU/$i.chr$j.ld >> $i.$j.job
qsub $i.$j.job
done
done
###########################
#CHINA
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X Y
do
echo \#PBS -N $i.CHB  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=~/hpc/db/hg19/1000Genome/CHINA.List.txt
vcf=~/hpc/db/hg19/1000Genome/chr$j.uni.vcf
echo plink --vcf $vcf --keep $ceu --r2 --ld-snp-list $i.pairSNP --out ./CHINA/$i.chr$j.ld >> $i.$j.job
qsub $i.$j.job
done
done
###########################
#JPT
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X Y
do
echo \#PBS -N $i.JPT  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=~/hpc/db/hg19/1000Genome/JPT.List.txt
vcf=~/hpc/db/hg19/1000Genome/chr$j.uni.vcf
echo plink --vcf $vcf --keep $ceu --r2 --ld-snp-list $i.pairSNP --out ./JPT/$i.chr$j.ld >> $i.$j.job
qsub $i.$j.job
done
done
