dmer_snp_db=/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt
for i in chr{1..22} chrX  chrY
do
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink/$i --keep $ceu --extract $dmer_snp_db --make-bed --out $i.dmer
done
