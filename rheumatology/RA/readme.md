Rheumatoid Arthritis


* Rheumatoid Arthritis Associated Variants


* Rheumatoid Arthritis Associated Genes
* 2018: PTPN22, TLR8, UBASH3A, IRF8, CREM, PAX5, BC017643, KLRC3, MEF2C, FAIM3, CXCR4, ID2, IL2, FOXP3, CD55
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA
for i in PTPN22 TLR8 UBASH3A IRF8 CREM PAX5 BC017643 KLRC3 MEF2C FAIM3 CXCR4 ID2 IL2 FOXP3 CD55
do
grep $i GHRA_ASA.hg19.bed
done
```
* Rheumatoid Arthritis Associated DNA Methylation Loci


* Rheumatoid Arthritis Associated Pathways


```
wget ftp://ftp.broadinstitute.org/pub/rheumatoid_arthritis/Stahl_etal_2010NG/RA_GWASmeta2_20090505-results.txt
wget ftp://ftp.broadinstitute.org/pub/rheumatoid_arthritis/Stahl_etal_2010NG/RA_GWASmeta2_20090505-results-README.txt
```



```
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/PTPN_PADI

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.hg19.bed
vip.gene<-read.table("immuneGWAS.hg19.bed")
data1<-read.table("miRNA_targets_loss_by_SNPs_in_seed_regions.txt",head=F,sep="\t")
vip.variant1<-data1[data1[,2]%in%vip.gene[,4],]
data2<-read.table("miRNA_targets_gain_by_SNPs_in_seed_regions.txt",head=T,sep="\t")
vip.variant2<-data2[data2[,2]%in%vip.gene[,4],]
colnames(vip.variant1)[2:5]=colnames(vip.variant2[,2:5])
vip.variant<-rbind(vip.variant1[,2:5],vip.variant2[,2:5])
length(table(vip.variant[,4]))
vip.utr.mirna.snp<-names(table(vip.variant[,4]))
write.table(vip.utr.mirna.snp,file="vip.utr3.mirna.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)


wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.hg19.bed
grep -v "[;|#|*|SNP|HLA]" immuneGWAS.hg19.bed | grep rs | awk '$2>0{print $1"\t"$2-50000"\t"$3+50000"\t"$6}' > immuneGWAS.hg19.uni.bed
bedtools intersect -wa -a ~/hpc/db/hg19/refGene.hg19.bed -b immuneGWAS.hg19.uni.bed > immuneGWAS.Gene.hg19.bed
awk '{print $5}' immuneGWAS.Gene.hg19.bed | sort -u | wc -l

panel="GHRA_ASA"
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/Epigene.hg19.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/VIP.hg19.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/InnateDB.hg19.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.Gene.hg19.bed

grep PTPN ~/hpc/db/hg19/refGene.hg19.bed > $panel.hg19.bed
grep PADI ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed
grep MUC ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed
grep HLA ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed

cat Epigene.hg19.bed >> $panel.hg19.bed
cat VIP.hg19.bed >> $panel.hg19.bed
cat InnateDB.hg19.bed >> $panel.hg19.bed
cat immuneGWAS.Gene.hg19.bed >> $panel.hg19.bed

sort -u $panel.hg19.bed > $panel.hg19.sort.bed
mv $panel.hg19.sort.bed $panel.hg19.bed
perl -p -i -e "s/chr//" $panel.hg19.bed

## Function Variants in Exom Regions
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
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
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed

### Fine-mapping to GWAS 
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.hg19.bed
```


Reference: 

[1]. Guo, S., etc, Genome-Wide DNA Methylation Patterns in Cd4+ T Cells from Chinese Han Patients with Rheumatoid Arthritis. Mod Rheumatol, 2017. 27(3): p. 441-447.



