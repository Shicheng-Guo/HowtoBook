

In the first step, we collected all the eQTL from whole blood, liver, colon, stomach and lung. 

* [14285](gnomad.genomes.r2.1.sites.rec.eQTL.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1%
* [772] (gnomad.exomes.r2.1.sites.rec.eQTL.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1%

```
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL
gunzip *.gz 

qval_threshold=0.05
data1<-subset(read.table("Whole_Blood.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data2<-subset(read.table("Liver.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data3<-subset(read.table("Small_Intestine_Terminal_Ileum.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data4<-subset(read.table("Stomach.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(as.character(data1[,19]),as.character(data2[,19]),as.character(data3[,19]),as.character(data4[,19]))
length(table(eqtl))
eqtl.snp<-names(table(eqtl))
write.table(eqtl.snp,file="vip.eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
```


```
#####################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP
perl -p -i -e 's/chr//g' eQTL.hg19.bed

panel="eQTL"
bed="innateDbUTR3.hg19.bed"
mkdir temp
## Function Variants in Genome Regions
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \# bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \# tabix -p vcf gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \&INFO/AF_eas\<0.999\' -T $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *genomes*rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
wc -l gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
#####################################################################
```
