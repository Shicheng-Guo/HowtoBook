wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz
gunzip cytoBand.txt.gz
 
data<-read.table("rheumatoid.gwas_catalog_v1.0-associations_e93_r2019-01-11.tsv",sep="\t")
cytoband<-read.table("~/hpc/db/hg19/cytoband.hg19.bed",sep="\t")



output<-cytoband[cytoband$V4 %in% names(table(data$V11)[table(data$V11)>9]),]
write.table(output,"RA-GWAS-Cytoband.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

cd /gpfs/home/guosa/hpc/rheumatology/RA/GWAS/cytobank
cp /gpfs/home/guosa/hpc/GWAS_Catalog/RA-GWAS-Cytoband.hg19.bed ./

perl -p -i -e 's/chr//g' RA-GWAS-Cytoband.hg19.bed
head RA-GWAS-Cytoband.hg19.bed

panel="RA-GWAS-Cytoband"
mkdir temp
## Function Variants in Exom Regions
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \&INFO/AF_eas\<0.999 \& \(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"frame\"\)\' -R $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
