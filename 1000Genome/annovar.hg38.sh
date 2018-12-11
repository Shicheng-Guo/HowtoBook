wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi
gunzip ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz





# split to each chrosome to speed up annotation
cd /gpfs/home/guosa/hpc/db/hg38/1000genome
for i in {1..22} X Y
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo vcftools --gzvcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --chr $i --recode --out ./annovar/chr22 >>chr$i.job
qsub chr$i.job
done

cd ~/hpc/tools/annovar/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/




