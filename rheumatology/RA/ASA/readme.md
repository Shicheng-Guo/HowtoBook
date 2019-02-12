Genetics of Rheumatoid Arthritis: 2019 Status
==========================================================

* [47](48-SNPs-Genetic-Risk-Score-EUR.txt)    SNPs in Genetic Risk Score in European (non-HLA 48 + HLA)

* 69    PADI,PTPN Variants in RA Risk

* 128   [GWAS-Meta Validated SNPs](GWAS-Meta-128-SNPs.20190208.vcf)

* 241   SNPs in the 3-UTR regions of VIP genes (miRNA-UTR3 double pairs)

* 334   [VIP functional variants](gnomad.exomes.r2.1.sites.rec.VIP.hg19.vcf.bed)

* 327   [eQTL for VIP genes](VIP.eQTL.327.SNPs.txt)

* 478   SNPs located in miRNAs seed region

* 907   [HLA functional SNPs](gnomad.exomes.r2.1.sites.rec.HLA.hg19.vcf.bed) (6:29671116-33106890)

* 1375  [Cytoband (>10 records)](../cytoband/1375.gnomad.exomes.r2.1.sites.rec.RA-GWAS-Cytoband.hg19.vcf.bed) from GWAS-Catalog-RA records

* 1763  [Epigene Functional Variants](gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf.bed)

* 2665  GWAS RA signficant; 789 SNPs; LD R2>0.8 SNPs (not include intron,HaploReg)

* 3325  GWAS Immnue System Desases

* 3591  [40K regions of GWAS Immnue System Desases](gnomad.exomes.r2.1.sites.rec.40KGWASAID.merge.vcf.bed)

* 2518  Functional SNPs located in 824 Immune-related genes

* 4215  miRNA-SNPs located in 824 Immune-related genes UTR3 regions

* 4728  Reactome database related immnue system genes (9 pathway and 1571 genes)

* 16013 CpGSNP-eQTL-Overlap Variants

* 6200  TFBS-SNPs in lncRNA promoter region

* 44326 CpGI-TFSB-DNase-BUR SNPs(CpGI.TFBS.DNase.BUR.hg19.bed)

==========================================================

Identify VIP gene biallelic SNPs to be genotyped (N=322) [code](VIP.biallelic.MFS.sh)

Identify all SNPs located in miRNAs (478 SNPs)
```
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

awk '{print $4"\t"$5}' miRNA_targets_gain_by_SNPs_in_seed_regions.txt | sort -u | grep rs > miRNA-SNPs.list.txt
awk '{print $4"\t"$5}' miRNA_targets_loss_by_SNPs_in_seed_regions.txt | sort -u | grep rs >> miRNA-SNPs.list.txt
awk '{print $5}' miRNA_targets_gain_by_SNPs_in_seed_regions.txt | sort -u | grep rs > miRNA-SNPs.snponly.list.txt
awk '{print $5}' miRNA_targets_loss_by_SNPs_in_seed_regions.txt | sort -u | grep rs >> miRNA-SNPs.snponly.list.txt
```

Identify EpiFactors gene bi-allelic SNPs to be genotyped (1763 SNPs)
```
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R Epigene.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.Epi.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.Epi.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.Epi.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf.bed
```

SNPs in the 3-UTR regions of VIP genes (241 SNPs)
```
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/miRNA-SNP
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

vip.gene<-read.table("/gpfs/home/guosa/hpc/db/Gnomad/vcf/VIP.hg19.bed")
data1<-read.table("miRNA_targets_loss_by_SNPs_in_seed_regions.txt",head=F,sep="\t")
vip.variant1<-data1[data1[,2]%in%vip.gene[,4],]
data2<-read.table("miRNA_targets_gain_by_SNPs_in_seed_regions.txt",head=T,sep="\t")
vip.variant2<-data2[data2[,2]%in%vip.gene[,4],]

colnames(vip.variant1)[2:5]=colnames(vip.variant2[,2:5])
vip.variant<-rbind(vip.variant1[,2:5],vip.variant2[,2:5])
length(table(vip.variant[,4]))
vip.utr.mirna.snp<-names(table(vip.variant[,4]))
write.table(vip.utr.mirna.snp,file="vip.utr3.mirna.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
```

eQTL of VIP genes (327 SNPs)
```
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL

qval_threshold=0.000000000000000000000005
data1<-subset(read.table("Whole_Blood.v7.egenes.txt",head=T,sep="\t"),qval<0.000000000000000000000001*qval_threshold)
data2<-subset(read.table("Liver.v7.egenes.txt",head=T,sep="\t"),qval<10*qval_threshold)
data3<-subset(read.table("Small_Intestine_Terminal_Ileum.v7.egenes.txt",head=T,sep="\t"),qval<100*qval_threshold)
data4<-subset(read.table("Stomach.v7.egenes.txt",head=T,sep="\t"),qval<0.000000000001*qval_threshold)
eqtl<-c(as.character(data1[,19]),as.character(data2[,19]),as.character(data3[,19]),as.character(data4[,19]))
length(table(eqtl))
vip.eqtl.snp<-names(table(eqtl))
dim(data1)
dim(data2)
dim(data3)
dim(data4)
write.table(vip.eqtl.snp,file="vip.eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
```

N=2518 Functional SNPs located in 824 Immune-related genes (csv): https://www.innatedb.com/moleculeSearch.do
```
perl -F"\s" -lane  "{print @F[6]}" InnateDB_genes.txt | grep -v name > InnateDB_genes.symbol.txt
### R
genesymbol<-read.table("InnateDB_genes.symbol.txt")
db<-read.table("~/hpc/db/hg19/refGene.hg19.bed",head=F)
rlt<-db[db[,5] %in% genesymbol[,1],]
write.table(rlt,file="InnateDB.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)
### R

perl -p -i -e "s/chr//" InnateDB.hg19.bed
mkdir temp
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R InnateDB.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.InnateDB.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf.bed
```

N=4215 miRNA-SNPs located in 824 Immune-related genes UTR3 regions(csv): https://www.innatedb.com/moleculeSearch.do
```
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

InnateDB<-read.table(file="InnateDB.hg19.bed",sep="\t")
data1<-read.table("miRNA_gain_by_SNPs_in_gene_3utr.txt",sep="\t")
data2<-read.table("miRNA_loss_by_SNPs_in_gene_3utr.txt",sep="\t")
SNP1<-data1[data1[,2] %in% InnateDB[,5],4]
SNP2<-data2[data2[,2] %in% InnateDB[,5],4]
SNP<-c(as.character(SNP1),as.character(SNP2))
write.table(sort(table(SNP),decreasing=T),file="InnateDB.UTR3.snp.txt",sep="\t",quote=F,col.names=F)
```

CpGI-TFBS-DNase-BUR-SNPs
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

N=3591 Functional SNPs within 10K up and down region of GWAS-Catalog autoimmnue SNPs
```
cd /gpfs/home/guosa/hpc/db/Gnomad/vcf/GWAS
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/GWAS-immnue-3325_SNP.hg19.bed
perl -p -i -e "s/chr//" GWAS-immnue-3325_SNP.hg19.bed
awk '{print $1"\t"$2-20000"\t"$3+20000"\t"$4}' GWAS-immnue-3325_SNP.hg19.bed > 5KGWASAID.hg19.bed
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R 5KGWASAID.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.5KGWASAID.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.5KGWASAID.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.5KGWASAID.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.5KGWASAID.merge.vcf.bed

```
N=6200 TFBS-SNPs in lncRNA promoter region
```
wget http://210.46.85.180:8080/SNP_linc_tfbs/download/linc_tfbs_chipseq_snp.txt
wget http://210.46.85.180:8080/SNP_linc_tfbs/download/linc_tfbs_snp.txt
awk '{print $2"\t"$17-1"\t"$17+1"\t"$16}' linc_tfbs_snp.txt | grep rs |sort -u > linc_tfbs_snp.hg19.bed
awk '{print $2"\t"$19-1"\t"$19+1"\t"$18}' linc_tfbs_chipseq_snp.txt |grep rs | sort -u > linc_tfbs_chipseq_snp.hg19.bed
```
N= 4728 Reactome pathway based Functional SNPs
```

PT2Symbol<-function(ReactomePathWayName){
  library("ReactomePA")
  library("graphite")
  library("GOSemSim")
  library("igraph")
org2org <- list(arabidopsis = "athaliana", bovine = "btaurus", 
                canine = "cfamiliaris", chicken = "ggallus", ecolik12 = "ecoli", 
                fly = "dmelanogaster", human = "hsapiens", mouse = "mmusculus", 
                pig = "sscrofa", rat = "rnorvegicus", celegans = "celegans", 
                xenopus = "xlaevis", yeast = "scerevisiae", zebrafish = "drerio")
p <- pathways(org2org[["human"]], "reactome")[[ReactomePathWayName]]
p <- convertIdentifiers(p, "symbol")
g <- pathwayGraph(p)
gg <- igraph.from.graphNEL(g)
gg <- as.undirected(gg)
gg <- set.graph.attribute(gg,"name","RING")
GeneSymbol <- sub("[^:]+:", "", V(gg)$name)
return(GeneSymbol)
}

x1<-PT2Symbol("Interleukin-1 signaling")
x2<-PT2Symbol("Interleukin-2 signaling")
x3<-PT2Symbol("Interleukin-10 signaling")
x4<-PT2Symbol("Interleukin-12 signaling")
x5<-PT2Symbol("Interleukin-6 signaling")
x6<-PT2Symbol("Interleukin-17 signaling")
x7<-PT2Symbol("Regulation of beta-cell development")
x8<-PT2Symbol("Adaptive Immune System")
x9<-PT2Symbol("Immune System")
genelist<-c(x1,x2,x3,x4,x5,x6,x7,x8,x9)

db<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/refGene.hg19.bed",sep="\t")
rlt<-db[db[,5]%in%genelist,]
write.table(rlt,file="ReactomePathWay.immnueGene.hg19.bed",sep="\t",quote=F,row.names = F,col.names = F)
```


