
We download 791 GWAS-Significant SNPs from GWAS Catalog and we collected all the linked SNPs with R2>0.6 in Asian Population. Totally, we identified xx SNPs with above method. 
```

for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done

cd /gpfs/home/guosa/hpc/db/hg19/CpGI
cp /gpfs/home/guosa/hpc/db/hg19/CpGI.hg19.bed ./
perl -p -i -e "s/chr//" CpGI.hg19.bed

mkdir temp

for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \# bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \# tabix -p vcf gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.005\' -R CpGI.hg19.bed  gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.vcf.bgz >>$i.job
echo bcftools sort gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.sort.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.sort.rmdup.vcf.bgz -Ov -o gnomad.genomes.r2.1.sites.chr$i.rec.CpGI.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.CpGI.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.bed

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip wgEncodeRegTfbsClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3"\t"$4}' wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClusteredV3.hg19.bed

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz
gunzip wgEncodeRegDnaseClusteredV3.bed.gz

cd /gpfs/home/guosa/hpc/db/hg19/CpGI
awk '{print "chr"$1"\t"$2-1"\t"$2+1"\t"$3"\t"$4"\t"$5}' gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.bed  > gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.hg19.bed
bedtools intersect -a gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed > CpGI.TFBS.SNP.hg19.txt 
bedtools intersect -wa -a CpGI.TFBS.SNP.hg19.sort.uni.txt -b wgEncodeRegDnaseClusteredV3.bed | sort -u > CpGI.TFBS.DNase.SNP.hg19.sort.uni.bed
bedtools intersect -wa -a CpGI.TFBS.DNase.SNP.hg19.sort.uni.bed -b /gpfs/home/guosa/hpc/db/hg19/BUR.GRCH37.hg19.bed | sort -u > CpGI.TFBS.DNase.BUR.hg19.bed
```


