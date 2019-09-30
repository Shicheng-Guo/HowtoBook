

cd ~/hpc/tools/RIblast/extdata
wget http://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
wget http://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz -O Homo_sapiens.GRCh38.mrna.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz

gunzip Homo_sapiens.GRCh38.ncrna.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz 
gunzip Homo_sapiens.GRCh38.mrna.fa.gz 
gunzip gencode.v31.transcripts.fa.gz

../RIblast db -i Homo_sapiens.GRCh38.mrna.fa  -o Homo_sapiens.GRCh38.mrna.fa.db
../RIblast db -i Homo_sapiens.GRCh38.cdna.all.fa  -o Homo_sapiens.GRCh38.cdna.all.fa.db
../RIblast db -i Homo_sapiens.GRCh38.ncrna.fa  -o Homo_sapiens.GRCh38.ncrna.fa.db
../RIblast db -i gencode.v31.transcripts.fa  -o gencode.v31.transcripts.fa.db

grep AC004585.1 Homo_sapiens.GRCh38.ncrna.fa.txt > AC004585.fa
../RIblast ris -i AC004585.fa -o AC004585.txt -d Homo_sapiens.GRCh38.cdna.all.fa.db


plink --file <input_prefix> --recode HV --snps-only just-acgt --out <output_prefix>
plink --bfile ROI.dbsnp --recode HV --out ROI.dbsnp
 
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
bcftools view -i '%QUAL>=20' calls.bcf
cSCC is usually manifested as a spectrum of lesions of progressive malignancy
#########################################################################################################
############################  hampmap phase 2   #####################################
#########################################################################################################
# NCBI36 or hg18
wget http://zzz.bwh.harvard.edu/plink/dist/hapmap_r23a.zip
wget http://zzz.bwh.harvard.edu/plink/dist/hapmap.pop
unzip hapmap_r23a.zip

# download database and script
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/ibdqc.pl

# rebuild plink file to avoid chromsome-miss-order problem
plink --bfile hapmap_r23a --make-bed --out hapmap_r23a.tab

# space to tab to generate bed files for liftOver from hg18 to hg19
plink --bfile hapmap_r23a.sort --recode tab --out hapmap_r23a.tab

# apply liftOverPlink.py to update hg18 to hg19 or hg38
./liftOverPlink.py -m hapmap_r23a.tab.map -p  hapmap_r23a.tab.ped -o hapmap_r23a.hg19 -c hg18ToHg19.over.chain.gz -e ./liftOver
./liftOverPlink.py -m hapmap_r23a.tab.map -p  hapmap_r23a.tab.ped -o hapmap_r23a.hg38 -c hg19ToHg38.over.chain.gz -e ./liftOver

# update plink to binary mode
plink --file hapmap_r23a.hg19 --make-bed --allow-extra-chr --out hapmap2.hg19
plink --file hapmap_r23a.hg38 --make-bed --allow-extra-chr --out hapmap2.hg38

# hapmap3 data cleaning and filtering
plink --bfile hapmap2.hg19 --missing
plink --bfile hapmap2.hg19 --maf 0.01 --make-bed --indep 50 5 2
plink2 --bfile hapmap2.hg19 --king-cutoff 0.125
plink2 --bfile hapmap2.hg19 --remove plink2.king.cutoff.out.id --make-bed -out hapmap2.hg19.deking
plink --bfile  hapmap2.hg19.deking --check-sex


# download database and script
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/ibdqc.pl
# download hapmap3 data in plink format
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/relationships_w_pops_051208.txt
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
# convert from hg18 to hg19 plink file
plink --bfile hapmap3_r1_b36_fwd_consensus.qc.poly.recode --recode --out hapmap3.hg18
./liftOverPlink.py -m hapmap3.hg18.map -p hapmap3.hg18.ped -o hapmap3.hg19 -c hg18ToHg19.over.chain.gz -e ./liftOver
./liftOverPlink.py -m hapmap3.hg18.map -p hapmap3.hg18.ped -o hapmap3.hg38 -c hg19ToHg38.over.chain.gz -e ./liftOver

# update plink to binary mode
plink --file hapmap3.hg19 --make-bed --allow-extra-chr --out hapmap3.hg19
plink --file hapmap3.hg38 --make-bed --allow-extra-chr --out hapmap3.hg38

# hapmap3 data cleaning and filtering
plink --bfile hapmap3.hg19 --missing
plink --file hapmap3.hg19 --maf 0.01 --make-bed --indep 50 5 2 --out hapmap3.hg19.indep
# plink --bfile hapmap3.hg19.indep --extract hapmap3.hg19.indep.in --genome --min 0.185
# perl ibdqc.pl plink
plink2 --bfile hapmap3.hg19 --king-cutoff 0.125
plink2 --bfile hapmap3.hg19 --remove plink2.king.cutoff.out.id --make-bed -out hapmap3.hg19.deking
plink --bfile  hapmap3.hg19.deking --check-sex


### merge personal data with hapmap2 and hapmap3
/home/guosa/hpc/db/hapmap3/hapmap3.hg19.deking
/home/guosa/hpc/db/hapmap2/hapmap2.hg19.deking
cd /home/guosa/hpc/db/hapmap2/
plink2 --bfile ~/hpc/db/hapmap2/hapmap2.hg19.deking --exclude RA2020hapmap2-merge.missnp --max-alleles 2 --make-bed --out hapmap2.hg19
plink2 --bfile ../RA2020-B8.dbsnp --max-alleles 2 --exclude RA2020hapmap2-merge.missnp --make-bed --out RA2020.hg19
plink --bfile hapmap2.hg19 --bmerge RA2020.hg19 --make-bed --out RA2020hapmap2
plink --bfile RA2020hapmap2 --maf 0.01 --geno 0.1 --make-bed --out RA2020hapmap2Update
plink --bfile RA2020hapmap2Update --make-bed --out RA2020hapmap2UpdatePCA
plink2 --bfile RA2020hapmap2UpdatePCA --pca --threads 31
hapmap2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hapmap2.pop",head=F)
head(hapmap2)
head(hapmap3)
eigenvec<-read.table("plink2.eigenvec",head=F)
head(eigenvec)
pop<-as.character(hapmap2[match(eigenvec[,2],hapmap2[,2]),3])
pop[is.na(pop)]<-"Guanghua"
eigenvec$pop=pop
eigenvec$col=as.numeric(as.factor(pop))
head(eigenvec)
pdf("pca.pdf")
plot(eigenvec[,3],eigenvec[,4],pch=16,col=as.factor(eigenvec$pop),xlab="principle component 1",ylab="principle componment 2",cex.axis=1.5,cex.lab=1.5)
legend("topright",pch=16,legend=unique(as.factor(pop)),col=unique(as.factor(eigenvec$pop)),bty="n",cex=1.5)
dev.off()
subset(eigenvec,pop=="Guanghua" & V3<0) # RA478

cd ~/hpc/rheumatology/RA/he2020/hapmap3
plink --bfile ~/hpc/db/hapmap3/hapmap3.hg19.deking --bmerge ../RA2020-B8.dbsnp --make-bed --out RA2020hapmap3
plink2 --bfile ~/hpc/db/hapmap3/hapmap3.hg19.deking --exclude RA2020hapmap3-merge.missnp --max-alleles 2 --make-bed --out hapmap3.hg19
plink2 --bfile ../RA2020-B8.dbsnp --max-alleles 2 --exclude RA2020hapmap3-merge.missnp --make-bed --out RA2020.hg19
plink --bfile hapmap3.hg19 --bmerge RA2020.hg19 --make-bed --out RA2020.hapmap3
plink --bfile RA2020.hapmap3 --maf 0.01 --geno 0.1 --make-bed --out RA2020.hapmap3.update
plink --bfile RA2020.hapmap3.update --make-bed --out RA2020.hapmap3.update.pca
plink2 --bfile RA2020.hapmap3.update.pca --pca --threads 31
### R
hapmap3<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hapmap3.pop",head=T)
head(hapmap3)
eigenvec<-read.table("plink2.eigenvec",head=F)
head(eigenvec)
pop<-as.character(hapmap3[match(eigenvec[,2],hapmap3[,2]),7])
pop[is.na(pop)]<-"SGH"
eigenvec$pop=pop
head(eigenvec)
pdf("Hapmap3.Guanghua.PCA.pdf")
plot(eigenvec[,3],eigenvec[,4],col=as.factor(eigenvec$pop),pch=as.numeric(as.factor(eigenvec$pop)),xlab="principle component 1",ylab="principle componment 2",cex.axis=1.5,cex.lab=1.5)
legend("topright",legend=unique(as.factor(pop)),pch=unique(as.numeric(as.factor(eigenvec$pop))),col=unique(as.factor(eigenvec$pop)),bty="n",cex=1)
dev.off()
subset(eigenvec,pop=="SGH" & V3<0) # RA478
pdf("Hapmap3.Guanghua.PCA.fine.pdf")
plot(eigenvec[,3],eigenvec[,4],col=as.factor(eigenvec$pop),pch=as.numeric(as.factor(eigenvec$pop)),xlab="principle component 1",ylab="principle componment 2",cex.axis=1.5,cex.lab=1.5,xlim=c(0.006,0.01),ylim=c(-0.0075,0.0025))
legend("topright",legend=unique(as.factor(pop)),pch=unique(as.numeric(as.factor(eigenvec$pop))),col=unique(as.factor(eigenvec$pop)),bty="n",cex=1)
dev.off()
subset(eigenvec,pop=="SGH" & V3<0.0075) # RA478
write.table(subset(eigenvec,pop=="SGH" & V3<0.0075),file="GH.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(subset(eigenvec,pop=="SGH" & V3<0.0075)[,1:2],file="pca.iid.remove",sep="\t",col.names=F,row.names=F,quote=F)

#########################################################################################################
############################  GWAS Study for Rheumatoid Arthritis   #####################################
#########################################################################################################
for input in MUC PADI CSF PTPN PLD DNMT TET ZNF MIR MAP GRI FOX IRF ST TTC COL GTF IL7
input="CCL21"
cd ~/hpc/rheumatology/RA/he2020
grep -w "CCL21\|CCR7" ~/hpc/db/hg19/refGene.hg19.V2.bed | awk '{print $1,$2-50000,$3+50000,$5}' OFS="\t" | sort -u | bedtools sort -i > ROI.hg19.bed
mkdir $input
plink --bfile RA2020-B8.dbsnp --extract range ROI.hg19.bed --make-bed --out ./$input/ROI
cd $input
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.R -O manhattan.plot.R
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/make.fancy.locus.plot.unix.R -O make.fancy.locus.plot.unix.R
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/localhit.pl -O localhit.pl
plink --bfile ROI --recode vcf --out ROI
bcftools view ROI.vcf -Oz -o ROI.vcf.gz
tabix -p vcf ROI.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ROI.vcf.gz -Ov -o ROI.refGene.vcf
plink --vcf ROI.vcf.gz --make-bed --out ROI.dbsnp
cp ROI.fam ROI.dbsnp.fam
sleep 30
touch readme.md
echo  "### Novel Genetic susceptibility loci in" $input "associated with rheumatoid arthritis" > readme.md
sleep 30
touch readme.txt
bcftools view -G ROI.refGene.vcf -Ov -o ROI.refGene.vcf.txt
echo  "Title: Novel Genetic susceptibility loci in" $input "associated with seropositive rheumatoid arthritis" > readme.txt
echo "" >> readme.txt
plink --bfile ROI.dbsnp --hardy --out ROI
echo  "ROI.hwe: Hardy-Weinberg test statistics for each SNP for" $input. "SNPs in control group should not signficant which means P should higher than 0.05 or 0.01. this table will be put to supplementary table" >> readme.txt
plink --bfile ROI.dbsnp --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI
echo  "ROI.assoc.logistic: logistic based case-control test for" $input. "default style is to test additive model --ADD-- in logistic regression. this file will be one of most important table in the manuscript" >> readme.txt
echo  "ROI.assoc.logistic.adjusted: this file include all the multiple-test corrected P-value for each SNPs in" $input. " When you prepare the manuscript, this file should be integrate with above file --ROI.assoc.logistic--" >> readme.txt
plink --bfile ROI.dbsnp --assoc --counts --adjust --ci 0.95 --out ROI
echo  "ROI.assoc: Chi-square based case-control test for" $input. "this file will be one of most important table in the manuscript since it showed the number of alleles in case and control" >> readme.txt
echo  "ROI.assoc.adjusted: this file include all the multiple-test corrected P-value --in the file: ROI.assoc-- for each SNPs in" $input. " When you prepare the manuscript, this file should be integrate with above file --ROI.assoc--" >> readme.txt
plink --bfile ROI.dbsnp --fisher --counts --adjust --ci 0.95 --out ROI
echo  'ROI.assoc.fisher: Fisher exact test based case-control association between SNPs and RA. This file will be useful when any cell lt 5. Usually when certain cell have number lt 5, we report fisher P-value not Chi-square P-value' >> readme.txt
echo  "ROI.assoc.fisher.adjusted: this file include all the multiple-test corrected P-value --in the file: ROI.assoc-- for each SNPs in" $input. " When you prepare the manuscript, this file should be integrate with above file --ROI.assoc--" >> readme.txt
plink --bfile ROI.dbsnp --model fisher --ci 0.95 --out ROI
echo  "ROI.model: Fisher's exact test based case-control association with different models for" $input. "this file is one of most important table in the manuscript" >> readme.txt
plink --bfile ROI.dbsnp --logistic --dominant --ci 0.95 --out ROI.dominant
echo  "ROI.dominant.assoc.logistic: logistic regression based association in dominant models for" $input. "this file is one of most important table in the manuscript. In the file of ROI.model, the DOM is based on fisher's exact test" >> readme.txt
plink --bfile ROI.dbsnp --logistic --recessive --ci 0.95 --out ROI.recessive
echo  "ROI.recessive.assoc.logistic: logistic regression based association in recessive models for" $input. "this file is one of most important table in the manuscript" >> readme.txt
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out haplotype --noweb
echo  'haplotype.assoc.hap: chi-square test based haplotype association. This file is important which can be shown with significant haplotype as a table in the manuscript' >> readme.txt
awk '$6=="NMISS"{print}' ROI.assoc.logistic > ROI.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.assoc.logistic >> ROI.assoc.logistic.add

sort -k12n,12 ROI.assoc.logistic | head -n 3
sort -k12n,12 ROI.assoc.logistic | tail -n 3

rs="rs78587665"
chr=17
grep $rs ROI.assoc.logistic
cp ROI.dbsnp.fam $rs.fam
plink --bfile ../RA2020-B8.dbsnp --snp $rs --recode --out $rs
plink --bfile ../RA2020-B8.dbsnp --r2 --ld-snp $rs --ld-window-kb 100 --ld-window 99999 --ld-window-r2 0 --out $rs
perl localhit.pl $rs > $rs.local
Rscript make.fancy.locus.plot.unix.R $rs $rs $chr $rs.local 4 0.05

rm *.R
rm *.pl
echo "" >> readme.txt
echo  'Abstract: The heritability of RA has been shown from twin studies to be 60%. Since 2007, rapid advances in technology underpinning the use of genome-wide association studies have allowed the identification of hundreds of genetic risk factors for many complex diseases. There are now >100 genetic loci that have been associated with RA. In the previous study, the contribution of HLA to heritability has been estimated to be 11–37% while 100 non-HLA loci were shown to explain 4.7% of the heritability of RA in Asians. The majority of the heritability is still missing.' $input ' have xxxxxx function which might be invovled in  pathology of RA. therefore, In this study, we conducted assocation study to investigate the role of xx and its paralog genes and in RA. in the first stage, we colllected 1078 seropositive RA and 1045 matched control. xxx SNPs in xx, xxx, xxx, xx were genotyped. We found SNPs rsxxx in xxx was signifciantly associated with RA, P=xxx, 95%CI.' >> readme.txt
echo "Run title:" $rs "in" $input "and Seropositive Rheumatoid Arthritis" >> readme.txt
echo "Reference: https://github.com/CNAID/Publication/blob/master/2018/1-s2.0-S155608641830090X-main.pdf" >>readme.txt
echo "Reference: https://github.com/CNAID/Publication/blob/master/2019/nejmoa18015622019.pdf" >>readme.txt


wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/UKBiobank_PheWAS_M05.txt
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/UKBiobank_PheWAS_M06.txt
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/UKBiobank_PheWAS_20002_1464.txt
zcat M05.assoc.tsv.gz | awk '$9<4.7e-3{print}' > M05.assoc.diff
zcat M06.assoc.tsv.gz | awk '$9<4.7e-3{print}' > M06.assoc.diff
zcat 20002_1464.assoc.tsv.gz | awk '$9<4.7e-3{print}' > 20002_1464.assoc.diff
awk '{print $2}'  M05.assoc.diff > RA.pheWAS.txt
awk '{print $2}'  M06.assoc.diff >> RA.pheWAS.txt
awk '{print $2}'  20002_1464.assoc.diff >> RA.pheWAS.txt
sort -u RA.pheWAS.txt > RA.pheWAS.uni.txt

wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/RA.pheWAS.uni.txt -O RA.pheWAS.uni.txt
wc -l RA.pheWAS.uni.txt
plink --bfile RA2020-B8.dbsnp --extract RA.pheWAS.uni.txt --assoc --counts --adjust --ci 0.95 --out pheWAS
awk '$1 !=6 && $3<0.05{print}' pheWAS.assoc.adjusted
#########################################################################################################
plink --bfile RA2020-B8 --recode vcf --out RA2020-B8q
bcftools view RA2020-B8.vcf -Oz -o RA2020-B8.vcf.gz
tabix -p vcf RA2020-B8.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/dbSNP/All_20180423.hg19.vcf.gz -c ID RA2020-B8.vcf.gz -Oz -o RA2020-B8.dbsnp.vcf.gz
plink --vcf RA2020-B8.dbsnp.vcf.gz --make-bed --out RA2020-B8.dbsnp
plink --bfile ROI.dbsnp --snp rs9972241 --assoc --counts --adjust --ci 0.95 --out rs9972241
plink --bfile RA2020-B8.dbsnp --snp rs9972241 --assoc counts --adjust --ci 0.95 --out rs9972241
plink --bfile RA2020-B8.dbsnp --snp rs6897932 --assoc counts --adjust --ci 0.95 --out rs6897932
plink --bfile RA2020-B8.dbsnp --snp rs6897932 --model counts --adjust --ci 0.95 --out rs6897932

data<-read.table("ROI.model",head=T,sep="",check.names=F)
sub<-subset(data,P<0.05)
write.table(sub,file="ROI.model.diff",col.names=T,row.names=F,quote=F,sep="\t")

plink --bfile RA2020-B8.dbsnp --from rs145044897 --to rs144054676 --assoc counts --adjust --ci 0.95 --out FOXP1
plink --bfile RA2020-B8.dbsnp --from rs145044897 --to rs144054676 --logistic --adjust --ci 0.95 --out FOXP1

wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.R -O manhattan.plot.R
awk '$6=="NMISS"{print}' ROI.assoc.logistic > ROI.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.assoc.logistic >> ROI.assoc.logistic.add
Rscript manhattan.plot.R ROI.assoc.logistic.add
#########################################################################################################
#########################################################################################################
dmp <- dmpFinder(beta, pheno = phen  , type = "categorical")
designMatrix <- model.matrix(~ phen)
dmr <- bumphunter(GRset.funnorm, design = designMatrix, cutoff = 0.2, B=0, type="Beta")	 

scp nu_guos@submit-1.chtc.wisc.edu:~/*bdg.gz ./

data<-read.table("plink.assoc.logistic",head=T,sep="")
newdata<-data[order(data$P),]
new<-subset(newdata,CHR !=6 & TEST=="ADD")
write.table(new[1:50,],file="plink.assoc.logistic.adjusted.sig.txt",col.names=T,row.names=F,quote=F,sep="\t")


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/pancrease/medip")
data<-read.table("R2.txt",head=T,row.names=1,sep="\t",check.names=F)
sampleList<-unique(unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1])))


use strict;
my @sam=qw/2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102/;
foreach my $sam(@sam){
my @fileC=glob("$sam\_N*narrowPeak");
my @fileT=glob("$sam\_T*narrowPeak");
foreach my $file(@fileT){
my ($out)=split/\_2019/,$file;
print "bedtools intersect -a $fileC[0] -b $file -v > $out.venn\n";
}
}

for i in `ls *venn`
do
awk '{print $1,$2,$3}' OFS="\t" $i > $i.bed
done
conda install -c bioconda intervene

for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
mkdir $i
intervene venn -i /gpfs/home/guosa/hpc/methylation/pancrease/medip/venn/$i*.bed --project $i
done

#########################################################################################################
######################  GWAS Study for Rheumatoid Arthritis   ##################################
#########################################################################################################
for input in MUC PADI CSF PTPN PLD DNMT TET ZNF MIR MAP GRI FOX IRF ST TTC COL GTF
input="GTF"
cd ~/hpc/rheumatology/RA/he2020
grep $input ~/hpc/db/hg19/refGene.hg19.V2.bed | awk '{print $1,$2-20000,$3+20000,$5}' OFS="\t" | sort -u | bedtools sort -i > ROI.hg19.bed
mkdir $input
plink --bfile RA2020-B8 --extract range ROI.hg19.bed --make-bed --out ./$input/ROI
cp ROI.hg19.bed ./$input/
cd $input
plink --bfile ROI --recode vcf --out ROI
bcftools view ROI.vcf -Oz -o ROI.vcf.gz
tabix -p vcf ROI.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/dbSNP/All_20180423.hg19.vcf.gz -c ID ROI.vcf.gz -Oz -o ROI.hg19.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ROI.hg19.vcf.gz -Oz -o ROI.hg19.refGene.vcf.gz
plink --vcf ROI.hg19.vcf.gz --make-bed --out ROI.dbsnp
cp ROI.fam ROI.dbsnp.fam
touch readme.txt
plink --bfile ROI.dbsnp --logistic --adjust --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --assoc --counts --adjust --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --fisher --counts --adjust --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --model fisher --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --logistic --genotypic --ci 0.95 --out ROI.genotype
plink --bfile ROI.dbsnp --logistic --dominant --ci 0.95 --out ROI.dominant
plink --bfile ROI.dbsnp --logistic --recessive --ci 0.95 --out ROI.recessive
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out haplotype --noweb

wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/manhattan.plot.R -O manhattan.plot.R
awk '$5=="ADD"{print}' ROI.assoc.logistic > plink.assoc.logistic.add
Rscript manhattan.plot.R plink.assoc.logistic.add

plink --bfile RA2020-B8 --recode vcf --out RA2020-B8
bcftools view RA2020-B8.vcf -Oz -o RA2020-B8.vcf.gz
tabix -p vcf RA2020-B8.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/dbSNP/All_20180423.hg19.vcf.gz -c ID RA2020-B8.vcf.gz -Oz -o RA2020-B8.dbsnp.vcf.gz
plink --vcf RA2020-B8.dbsnp.vcf.gz --make-bed --out RA2020-B8.dbsnp
plink --bfile ROI.dbsnp --snp rs9972241 --assoc --counts --adjust --ci 0.95 --out rs9972241

plink --bfile RA2020-B8.dbsnp --snp rs9972241 --assoc counts --adjust --ci 0.95 --out rs9972241
plink --bfile RA2020-B8.dbsnp --from rs145044897 --to rs144054676 --assoc counts --adjust --ci 0.95 --out FOXP1
plink --bfile RA2020-B8.dbsnp --from rs145044897 --to rs144054676 --logistic --adjust --ci 0.95 --out FOXP1

wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.R -O manhattan.plot.R
awk '$6=="NMISS"{print}' ROI.assoc.logistic > ROI.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.assoc.logistic >> ROI.assoc.logistic.add
Rscript manhattan.plot.R ROI.assoc.logistic.add
#########################################################################################################
#########################################################################################################
rm ROI.hg19.bed 
for input in CUBN ASTL BMP1 C1R C1S CDCP2 CNTNAP1 CNTNAP2 CNTNAP3 CNTNAP3B CNTNAP3C CNTNAP4 CNTNAP5 CP CSMD1 CSMD2 CSMD3 CUZD1 DCBLD2 DLL3 DMBT1 EDIL3 ENSG00000286088 EYS F5 F8 F9 FBN3 HEPH HEPHL1 LRP1B MASP1 MEP1A MEP1B MFGE8 MFRP NETO1 NETO2 NOTCH2 NOTCH3 NRP1 NRP2 NRXN1 NRXN2 NRXN3 OVCH2 PCOLCE PCOLCE2 PDGFC PDGFD SCUBE2 SEZ6 SEZ6L SEZ6L2 SLIT2 SLIT3 TLL1 TLL2 VWDE 

rm ROI.hg19.bed 
for input in PITX2 ALX1 ALX3 ALX4 ARGFX ARX CRX DLX6 DMBX1 DPRX DRGX DUX4 DUX4L2 DUXA DUXB ESX1 EVX1 EVX2 GSC GSC2 HESX1 HOPX HOXA4 HOXB1 ISX MIXL1 MNX1 MSX1 NANOGP8 NKX2-5 NOBOX OTP OTX1 OTX2 PAX1 PAX2 PAX3 PAX4 PAX5 PAX6 PAX7 PAX8 PAX9 PHOX2A PHOX2B PITX1 PITX3 POU6F1 PROP1 PRRX1 PRRX2 RAX RAX2 RHOXF1 RHOXF2 RHOXF2B SEBOX SHOX SHOX2 TPRX2P UNCX VSX2 

do
grep $input ~/hpc/db/hg19/refGene.hg19.V2.bed | awk '{print $1,$2-20000,$3+20000,$5}' OFS="\t" | sort -u | bedtools sort -i >> ROI.hg19.bed
done
input="PITX2"
mkdir $input
plink --bfile RA2020-B8 --extract range ROI.hg19.bed --make-bed --out ./$input/ROI
cp ROI.hg19.bed ./$input/
cd $input
plink --bfile ROI --recode vcf --out ROI
bcftools view ROI.vcf -Oz -o ROI.vcf.gz
tabix -p vcf ROI.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/dbSNP/All_20180423.hg19.vcf.gz -c ID ROI.vcf.gz -Oz -o ROI.hg19.vcf.gz
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') ROI.hg19.vcf.gz -Oz -o ROI.hg19.refGene.vcf.gz
plink --vcf ROI.hg19.vcf.gz --make-bed --out ROI.dbsnp
cp ROI.fam ROI.dbsnp.fam
touch readme.txt
plink --bfile ROI.dbsnp --logistic --adjust --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --assoc --counts --adjust --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --fisher --counts --adjust --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --model fisher --ci 0.95 --out ROI
plink --bfile ROI.dbsnp --logistic --genotypic --ci 0.95 --out ROI.genotype
plink --bfile ROI.dbsnp --logistic --dominant --ci 0.95 --out ROI.dominant
plink --bfile ROI.dbsnp --logistic --recessive --ci 0.95 --out ROI.recessive
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out haplotype --noweb
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/manhattan.plot.R -O manhattan.plot.R
awk '$5=="ADD"{print}' ROI.assoc.logistic > plink.assoc.logistic.add
Rscript manhattan.plot.R plink.assoc.logistic.add
#########################################################################################################
#########################################################################################################	 
plink --bfile ROI.dbsnp.fam --logistic --covar plink.eigenvec --covar-number 1-20 --adjust
cd ~/hpc/rheumatology/RA/he2020
plink --bfile RA2020-B8.dbsnp --snp rs35469986 --window 50 --make-bed --out ./TAB1/TAB1
plink --bfile TAB1 --logistic --covar ../plink.eigenvec --covar-number 1-20 --adjust
plink --bfile TAB1 --r2 --ld-snp rs35469986 --ld-window-kb 100 --ld-window 99999  --ld-window-r2 0
perl local.pl  > local.new
wc -l plink.ld
grep ADD plink.assoc.logistic | wc -l 

mkdir temp
wget https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/hg19.chrom.sizes -O hg19.chrom.sizes

for i in `ls *.bdg`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # gunzip $i >> $i.job
echo bedtools sort -i $i \> $i.sort >> $i.job
echo bedGraphToBigWig $i.sort hg19.chrom.sizes $i.bw >> $i.job
qsub $i.job
done

for i in `ls *.bw`
do
for j in `ls *.bw`
do
echo \#PBS -N $i.$j  > $i.$j.job
echo \#PBS -l nodes=1:ppn=1 >> $i.$j.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.$j.job
echo \#PBS -m abe  >> $i.$j.job
echo \#PBS -o $(pwd)/temp/ >>$i.$j.job
echo \#PBS -e $(pwd)/temp/ >>$i.$j.job
echo cd $(pwd) >> $i.$j.job
echo bigWigCorrelate $i $j >> $i.$j.job
qsub $i.$j.job
done
done

use strict;
my @file=glob("*.e*");
foreach my $file(@file){
if($file ~= /_(\w\d*)_2019\.+_(\w\d*)_2019/){
print "$1\t$2\n";
}
}

wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.minfi.R -O manhattan.qqplot.minfi.R
Rscript manhattan.qqplot.minfi.R dmr.bumphunter.minfi.hg19.bed.txt
cp cutadapt ~/miniconda3/bin/cutadapt /gpfs/home/guosa/tools/TrimGalore-0.5.0/
grep "Case_to_Control.logFC\|PAX2" AtrialFibrillation.CaseControl.myDMP.adjp.0.5.txt
#########################################################################################################
######################  GWAS Study for Rheumatoid Arthritis   ##################################
#########################################################################################################
cd ~/hpc/rheumatology/RA/he2020
plink --bfile RA2020 --mind 01 --make-bed --out RA2020-B1
plink --bfile RA2020-B1 --geno 0.1 --make-bed --out RA2020-B2
plink --bfile RA2020-B2 --maf 0.01 --make-bed --out RA2020-B3
plink --bfile RA2020-B3 --hwe 0.00001 --make-bed --out RA2020-B4
plink2 --bfile RA2020-B4 --king-cutoff 0.125
plink2 --bfile RA2020-B4 --remove plink2.king.cutoff.out.id --make-bed -out RA2020-B5
plink --bfile RA2020-B5 --check-sex
plink --bfile RA2020-B5 --impute-sex --make-bed --out RA2020-B6
plink --bfile RA2020-B6 --check-sex
grep PROBLEM plink.sexcheck | awk '{print $1,$2}' > sexcheck.remove
plink --bfile RA2020-B6 --remove sexcheck.remove --make-bed --out RA2020-B7
plink --bfile RA2020-B7 --test-missing midp 
awk '$5<0.000001{print}' plink.missing | awk '{print $2}' > missing.imblance.remove
plink --bfile RA2020-B7 --exclude missing.imblance.remove --make-bed --out RA2020-B8
plink --bfile RA2020-B8 --pca --threads 31
perl phen.pl RA2020-B8.fam > RA2020-B8.fam.new
mv RA2020-B8.fam.new RA2020-B8.fam
plink --bfile RA2020-B8 --logistic --covar plink.eigenvec --covar-number 1-5 --adjust
plink --bfile RA2020-B8 --assoc mperm=1000000 --adjust gc --threads 31
grep "ADD\|NMISS" plink.assoc.logistic > plink.assoc.logistic.add
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/manhattan.plot.R -O manhattan.plot.R
Rscript manhattan.plot.R plink.assoc.logistic.add
## local
head -n 50 plink.assoc.logistic.adjusted

rs="rs17499655"
chr="6"
plink --bfile RA2020-B8.dbsnp --snp $rs --window 500 --make-bed --out $rs
plink --bfile $rs --logistic --covar plink.eigenvec --covar-number 1-5 --adjust --out $rs
plink --bfile $rs --r2 --ld-snp $rs --ld-window-kb 1000 --ld-window 99999  --ld-window-r2 0  --out $rs
perl local.pl  > $rs
wc -l plink.ld
grep "ADD\|NMISS" plink.assoc.logistic > plink.assoc.logistic.add
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/localhit.pl -O localhit.pl
perl localhit.pl $rs > $rs.local
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/make.fancy.locus.plot.unix.R -O make.fancy.locus.plot.unix.R
Rscript make.fancy.locus.plot.unix.R rs9972241 rs9972241 $chr $rs.local 20 0.00005


grep "PL" ~/hpc/db/hg19/refGene.hg19.V2.bed | awk '{print $1,$2-10000,$3+10000,$5}' OFS="\t" | sort -u | bedtools sort -i > PLD.hg19.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/RA2020/PLD/Phospholipase.txt -O Phospholipase.txt
grep -w -f Phospholipase.txt PLD.hg19.bed > PLD.hg19.vcf.bed
input="RA2020-B8.vcf"
bcftools view $input -Oz -o $input.gz
tabix -p vcf $input.gz
bcftools annotate -a ~/hpc/db/hg19/dbSNP/All_20180423.hg19.vcf.gz -c ID $input.gz -Oz -o $input.hg19.gz
plink --bfile RA2020-B8.dbsnp --pca --threads 31
plink --bfile RA2020-B8.dbsnp --logistic --covar plink.eigenvec --covar-number 1-20 --adjust
cd ~/hpc/rheumatology/RA/he2020
plink --bfile RA2020-B8.dbsnp --snp rs35469986 --window 50 --make-bed --out ./TAB1/TAB1
plink --bfile TAB1 --logistic --covar ../plink.eigenvec --covar-number 1-20 --adjust
plink --bfile TAB1 --r2 --ld-snp rs35469986 --ld-window-kb 100 --ld-window 99999  --ld-window-r2 0
perl local.pl  > local.new
wc -l plink.ld
grep ADD plink.assoc.logistic | wc -l 

#########################################################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3
plink --vcf chr22.dose.clean.hg19.vcf.gz --keep ../../RA2020-B8.dbsnp.fam --make-bed --out chr22.dose.dbsnp
perl phen.pl chr22.dose.dbsnp.fam > chr22.dose.dbsnp.fam.2
mv  chr22.dose.dbsnp.fam.2  chr22.dose.dbsnp.fam
plink --bfile chr22.dose.dbsnp --pca --threads 31
plink --bfile chr22.dose.dbsnp --logistic --covar plink.eigenvec --covar-number 1-20 --adjust

# PBS START
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # plink --vcf chr$i.dose.clean.hg19.vcf.gz --keep ../../RA2020-B8.dbsnp.fam --make-bed --out ./plink/chr$i.dose.dbsnp >>$i.job
echo # cd ./plink >>$i.job
echo # perl phen.pl chr$i.dose.dbsnp.fam \> chr$i.dose.dbsnp.fam.2 >>$i.job
echo # mv chr$i.dose.dbsnp.fam.2  chr$i.dose.dbsnp.fam >>$i.job
echo # plink --bfile chr$i.dose.dbsnp --pca --threads 31 --out chr$i >>$i.job
echo # plink --bfile chr$i.dose.dbsnp --logistic --covar chr$i.eigenvec --covar-number 1-20 --adjust --out chr$i >>$i.job
echo grep \"ADD\\\|NMISS\" chr$i.assoc.logistic \> chr$i.assoc.logistic.add >>$i.job
qsub $i.job
done

PLD4, chr14:105391,187-105399,573

plink --bfile chr14.dose.dbsnp --chr 14 --from-kb 105389 --to-kb 105400 --logistic --covar chr14.eigenvec --covar-number 1-20
grep "ADD\|NMISS" plink.assoc.logistic > plink.assoc.logistic.add
plink --bfile plink --snp $rs --window 50 --make-bed
rs=""
chr=""
plink --bfile RA2020-B8.dbsnp --snp $rs --window 50 --make-bed
plink --bfile plink --logistic --covar $chr.eigenvec --covar-number 1-20 --adjust
plink --bfile plink --r2 --ld-snp $rs --ld-window-kb 100 --ld-window 99999  --ld-window-r2 0
perl local.pl  > $rs
wc -l plink.ld
grep "ADD\|NMISS" plink.assoc.logistic > plink.assoc.logistic.add

grep rs191997515 chr7.assoc.logistic.add
grep rs191997515 plink.assoc.logistic.add

plink --bfile RA2020-B8.dbsnp --snp rs191997515 --hardy
plink --bfile chr7.dose.dbsnp --snp rs191997515 --hardy
#########################################################################################################
#########################################################################################################

# PBS END
for i in {15..61}
do
qdel 5602$i.bright
done

#########################################################################################################

#####################################
###  local.pl
#####################################
my %db;
my $rs=shift @ARGV;
open DB1,"$rs.assoc.logistic";
while(<DB1>){
chomp;
next if $_ !~ /ADD/i;
my @line=split/\s+/;
$db{$line[2]}="$line[2]\t$line[3]\t$line[9]\ttyped";
}
close DB1;
my %snp;
open DB2, "$rs.ld";
print "SNP\tPOS\tPVAL\tTYPE\tRSQR\n";
while(<DB2>){
chomp;
next if /CHR_A/;
my @line=split /\s+/,$_;
next if defined $snp{$line[6]};
print "$db{$line[6]}\t$line[7]\n";
$snp{$line[6]}=$line[6];
}
#####################################


data<-read.table("plink.assoc.logistic.add",head=T,sep="",check.names=F)
ss<-subset(data,P<10^-7)
write.table(ss,file="plink.assoc.logistic.add.significant",sep="\t",quote=F,col.names=T,row.names=F)

wget https://www.broadinstitute.org/files/shared/diabetes/scandinavs/known_genes_build35_050307.tar.gz
wget http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/HapmapII_GRCh37_RecombinationHotspots.tar.gz

my %db;
open DB1,"plink.assoc.logistic";
while(<DB1>){
chomp;
next if $_=! ~/ADD/i;
my @line=split/\s+/;
$db{$line[1]}="$line[1]\t$line[2]\t$line[8]\ttyped";
}
close DB;
open F, "plink.ld";
while(<F>){
chomp;
my @line=split /\s+/,$_;
print "$db{$line[1]}\t$line[5]\t$line[6]";
}


		  
plink --vcf chr22.dose.clean.hg19.vcf.gz --keep ../../RA2020-B8.dbsnp.fam --make-bed --out chr22.dose.dbsnp
perl phen.pl chr22.dose.dbsnp.fam > chr22.dose.dbsnp.fam.2
mv  chr22.dose.dbsnp.fam.2  chr22.dose.dbsnp.fam
plink --bfile chr22.dose.dbsnp --pca --threads 31
plink --bfile chr22.dose.dbsnp --logistic --covar plink.eigenvec --covar-number 1-20 --adjust

plink --file data --keep mylist.txt


/gpfs/home/guosa/hpc/methylation/clep/wgbs
perl ~/bin/smartbismark.pl --input SraRunTable.txt --genome hg19 --phred=33 --server MCRI --queue shortq --submit no
for i in `ls *.pbs`
do
qsub $i
done

for i in SRR98883{01..04}.pbs
do
qsub $i
done



mkdir temp
for i in `ls *.fastq`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo gzip $i >> $i.job
qsub $i.job
done


for i in SRR98883{05..41}
do
fastq-dump --skip-technical --split-files --gzip $i &  
echo $i
done

mkdir temp
for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # fastq-dump -I --split-files $i >> $i.job
echo # bismark --bowtie2 --phred33-quals -p 6 --fastq -L 32 -N 0 -D 5 -R 1 ~/hpc/db/hg19/bismark -1 $i\_1.fastq -2 $i\_2.fastq -o ./  >>$i.job
echo gzip $i
qsub $i.job
done

grep 'C methylated in CpG context' *.txt
grep 'C methylated in CHG context' *.txt
grep 'C methylated in CHH context' *.txt


cp R05_S18_L001_R1_001_00_bismark_bt2_pe.bam
cp B01_S21_L001_R1_001_00_bismark_bt2_pe.bam

~/src/sambamba/sambamba_v0.3.3 view -h -t $numThreads -s $fractionOfReads -f bam --subsampling-seed=$seed $testBam -o $subsampledTestBam

for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126

for i in SRR9888301 SRR9888302 SRR9888303 SRR9888304
do
fastq-dump -I --split-files $i &
done 

for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126
do
seqtk sample -s100 $i\_1.fastq 5000000 > ./subsampling/RRS$i.fq
seqtk sample -s100 $i\_1.fastq 5000000 > ./subsampling/RRS$i.fq
echo $i
done

for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126
do
head -n 20000000  $i\_1.fastq > ./subsampling/RRS$i\_1.fq
head -n 20000000  $i\_2.fastq > ./subsampling/RRS$i\_2.fq
echo $i
done

for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126
do
bismark --bowtie2 --phred33-quals -p 6 --fastq -L 32 -N 0 -D 5 -R 1 ~/hpc/db/hg19/bismark -1 RRS$i\_1.fq -2 RRS$i\_2.fq -o ./
done

for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126
do
echo RRS$i
done


ftp -i ftp-private.ncbi.nlm.nih.gov
Username: subftp
Password: w4pYB9VQ
cd uploads/shicheng.guo_hotmail.com_45OXdp2r
mkdir new_folder
cd new_folder
put file_name




Circulating cfDNA based Low-pass WGBS in hepatocellular carcinoma surveillance 

for i in SRR10070129 SRR10070128 SRR10070127 SRR10070126
do
fastq-dump -I --split-files $i &
done


for i in `<SRR_Acc_List.txt`
do
fastq-dump -I --split-files $i &
done


#!/usr/bin/perl

git clone https://github.com/taoliu/MACS.git
cd ~/hpc/tools/MACS
python setup.py install --user

pip install MACS2
pip install -U MACS2


scp -o ProxyCommand='ssh -A  ssh -A B -W %h:%p' /tmp/a C:/tmp/a


scp id_rsa.pub guosa@localhost:./
ssh -p 8809 guosa@localhost

ssh-keygen -t rsa
scp -P 8809 /home/biostaff/.ssh/id_rsa.pub  guosa@localhost:~/.ssh/authorized_keys
ssh -p 8809 guosa@localhost
chmod 600 ~/.ssh/authorized_keys
chmod 700 ~/.ssh/

ssh-keygen -t rsa
scp -P 37122 /home/nu_guos/.ssh/id_rsa.pub biostaff@ada.zettadom.com:~/.ssh/authorized_keys
ssh -p 37122 biostaff@ada.zettadom.com
chmod 600 ~/.ssh/authorized_keys
chmod 700 ~/.ssh/

# BIRC10-LC
ssh-keygen -t rsa
scp ~/.ssh/id_rsa.pub nu_guos@submit-1.chtc.wisc.edu:~/.ssh/authorized_keys
ssh nu_guos@submit-1.chtc.wisc.edu
chmod 600 ~/.ssh/authorized_keys
chmod 700 ~/.ssh/

for i in `ls *.r`
do
Rscript $i
done

scp nu_guos@submit-1.chtc.wisc.edu:~/*bdg.gz ./

scp -o "ProxyJump -P 37122 biostaff@ada.zettadom.com" -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/*.sorted ./


scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/narrowpeak.tar.gz ./
scp -P 37122 biostaff@ada.zettadom.com:~/bin/* ./
scp nu_guos@submit-1.chtc.wisc.edu:~/narrowpeak.tar.gz ./

scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/bdg.tar.gz ./

scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/*_N_*.bdg.gz ./
scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/*_T1_*.bdg.gz ./
scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/*_T2_*.bdg.gz ./
scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/*_T3_*.bdg.gz ./
scp -P 8809 guosa@localhost:/home/guosa/hpc/pancreatic/*_T4_*.bdg.gz ./

scp -P 37122 biostaff@ada.zettadom.com:~/bin/* ./
 
 
 ChAMPdata’, ‘IlluminaHumanMethylationEPICanno.ilm10b2.hg19
 
 
qsub SRR9888301.pbs
qsub SRR9888302.pbs
qsub SRR9888303.pbs
qsub SRR9888304.pbs
qsub SRR9888305.pbs
qsub SRR9888306.pbs
qsub SRR9888307.pbs
qsub SRR9888308.pbs
qsub SRR9888309.pbs
qsub SRR9888310.pbs
qsub SRR9888311.pbs
qsub SRR9888312.pbs
qsub SRR9888313.pbs
qsub SRR9888314.pbs
qsub SRR9888315.pbs
qsub SRR9888316.pbs
qsub SRR9888317.pbs
qsub SRR9888318.pbs

fastq_screen 
fastq_screen --bisulfite --get_genomes

qsub SRR9888319.pbs
qsub SRR9888320.pbs
qsub SRR9888321.pbs
qsub SRR9888322.pbs
qsub SRR9888323.pbs
qsub SRR9888324.pbs
qsub SRR9888325.pbs
qsub SRR9888326.pbs
qsub SRR9888327.pbs
qsub SRR9888328.pbs
qsub SRR9888329.pbs
qsub SRR9888330.pbs
qsub SRR9888331.pbs
qsub SRR9888332.pbs
qsub SRR9888333.pbs
qsub SRR9888334.pbs
qsub SRR9888335.pbs
qsub SRR9888336.pbs
qsub SRR9888337.pbs
qsub SRR9888338.pbs
qsub SRR9888339.pbs
qsub SRR9888340.pbs
qsub SRR9888341.pbs

 

Library preparation method (Swift, TruSeq, QIAseq)
DNA input (200ng, 50ng)
DNA fragmentation before bisulfite treatment (Yes, No)
Bisulfite conversion kit (EZ DNA Methylation-Gold Kit, EpiTect Fast Bisulfite Kit)
Sequencing Platform (HiSeq X)
Spike-in (%) PhiX(5%) or WGS(10%-90%)
PCR enzyme (Enzyme R3, FailSafe, VeraSeq Ultra)
PCR cycles (6, 7 9 10 12)
Average data output per sample ( 66-200 million)


Host jump
    HostName ada.zettadom.com
    User biostaff
    Port 37122
Host localhost
    HostName guosa
    Port 8809
    ProxyJump jump

scp guosa@localhost:/home/guosa/hpc/pancreatic/* ./

echo "ls -l" | at 14:27:22

for i in `ls *sorted`
do
macs2 callpeak -t $i -f BAM -g hs -n $i.macs -B -q 0.05  &
done




wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://code.google.com/p/bedtools/


BiocManager::install("ChIPseeker")
BiocManager::install("SeqPlots")
BiocManager::install("Genomation")

cd /home/guosa/hpc/project/pmrp/phase1/plink
cp /home/guosa/hpc/project/pmrp/phase2/S_Hebbring_Unr.Guo.bim ./
cp /home/guosa/hpc/project/pmrp/phase2/S_Hebbring_Unr.Guo.fam ./
cp /home/guosa/hpc/project/pmrp/phase2/S_Hebbring_Unr.Guo.bed ./
cp /home/guosa/hpc/project/pmrp/phase1/plink/FinalRelease_QC_20140311_Team1_Marshfield.bim ./
cp /home/guosa/hpc/project/pmrp/phase1/plink/FinalRelease_QC_20140311_Team1_Marshfield.fam ./
cp /home/guosa/hpc/project/pmrp/phase1/plink/FinalRelease_QC_20140311_Team1_Marshfield.bed ./


  cp ~/hpc/db/hg19/hg19.fa ./
  perl -p -i -e 's/chr//g' hg19.fa
  gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict 
  cat hg19.dict 
  samtools faidx hg19.fa
  cat hg19.fa.fai 
 
  bgzip -c file.vcf > file.vcf.gz
  tabix -p vcf file.vcf.gz

scp nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/All_20180423* ./

scp nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/GCF_0* ./

gawk -v RS="(\r)?\n" 'BEGIN { FS="\t" } !/^#/ { if ($10 != "na") print $7,$10; else print $7,$5 }' GCF_000001405.38_GRCh38.p12_assembly_report.txt > dbSNP-to-UCSC-GRCh38.p12.map
bcftools annotate --rename-chrs dbSNP-to-UCSC-GRCh38.p12.map GCF_000001405.38.gz | gawk '/^#/ && !/^##contig=/ { print } !/^#/ { if( $1!="na" ) print }' | bgzip -@4 -l9 -c > GCF_000001405.38.dbSNP152.GRCh38p12b.GATK.vcf.gz


wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.gz
wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.gz.tbi
wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_assembly_report.txt 
awk -v RS="(\r)?\n" 'BEGIN { FS="\t" } !/^#/ { if ($10 != "na") print $7,$10; else print $7,$5 }' GCF_000001405.38_GRCh38.p12_assembly_report.txt > dbSNP-to-UCSC-GRCh38.p12.map
perl -p -i -e '{s/chr//}' dbSNP-to-UCSC-GRCh38.p12.map
bcftools annotate --rename-chrs dbSNP-to-UCSC-GRCh38.p12.map GCF_000001405.38.gz | gawk '/^#/ && !/^##contig=/ { print } !/^#/ { if( $1!="na" ) print }' | bgzip -c > GCF_000001405.38.dbSNP153.GRCh38p12b.GATK.vcf.gz
python CrossMap.py vcf hg38Tohg19.over.chain.gz GCF_000001405.38.dbSNP153.GRCh38p12b.GATK.vcf.gz ~/hpc/db/hg19/hg19.fa  GCF_000001405.38.dbSNP153.hg19.gz


/home/guosa/hpc/tools/jre1.8.0_221/bin/java -jar 

wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.gz
wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.gz.tbi
wget https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/GCF_000001405.25_GRCh37.p13_assembly_report.txt -O GCF_000001405.25_GRCh37.p13_assembly_report.txt
awk -v RS="(\r)?\n" 'BEGIN { FS="\t" } !/^#/ { if ($10 != "na") print $7,$10; else print $7,$5 }' GCF_000001405.25_GRCh37.p13_assembly_report.txt > dbSNP-to-UCSC-GRCh37.p13.map
perl -p -i -e '{s/chr//}' dbSNP-to-UCSC-GRCh37.p13.map

bcftools annotate --rename-chrs dbSNP-to-UCSC-GRCh37.p13.map GCF_000001405.25.gz | gawk '/^#/ && !/^##contig=/ { print } !/^#/ { if( $1!="na" ) print }' | bgzip -c > dbSNP153.hg19.vcf.gz

time bcftools view GCF_000001405.25.gz -r NC_000001.10 -Oz -o dbSNP153.chr1.hg19.vcf.gz
time tabix GCF_000001405.25.gz NC_000001.10 | bgzip -c > dbSNP153.chr1.hg19.vcf.gz



mkdir temp
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/contigReplace.pl
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view GCF_000001405.25.gz -r chr$i -Oz -o dbSNP152.chr$i.hg19.vcf.gz >>$i.job
echo tabix -p vcf chr$i.dose.contig.vcf.gz >>$i.job
qsub $i.job
done


input="RA2020-B8.vcf"
bcftools view $input -Oz -o $input.gz
tabix -p vcf $input.gz
bcftools annotate -a ~/hpc/db/hg19/dbSNP/All_20180423.hg19.vcf.gz -c ID $input.gz -Oz -o $input.hg19.gz


source("")
grep ADD 
d<-read.table("ADD.p.txt")
colnames(mylimma)=c("CHR","probeID","MAPINFO","","","","","","P.Value")
ManhattanPlot(mylimma)



exmAims<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/exmAIMs.txt")
db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/Illumina_CoreExome_Beadchip.hg19.exm2rs.bed.txt")
rsAims<-db[na.omit(match(exmAims[,1],db$V3)),]
write.table(rsAims,file="rsAims.bed",quote=F,row.names = F,col.names = F,sep="\t")
write.table(rsAims[,6],file="rsAims.txt",quote=F,row.names = F,col.names = F,sep="\t")

cd /home/guosa/hpc/rheumatology/RA/he2020
input="RA2020"
plink –-noweb --file $input --make-bed --out $input
plink --bfile $input --maf 0.1 --check-sex --make-bed --out $input.1
wget https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/Gender.R -O Gender.R
wget https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/rsAims.txt -O rsAims.txt
Rscript Gender.R $input.1.sexcheck sex.checking.jpeg
plink --bfile $input.1 --extract rsAims.txt --recode --out AIMs

perl -p -i -e '{/-9/1/g}' AIMs.map

cd ~/hpd/tools/
git clone https://github.com/chrchang/eigensoft.git

vim par.PED.EIGENSTRAT

genotypename:    AIMs.ped
snpname:         AIMs.map
indivname:       AIMs.ped
outputformat:    EIGENSTRAT
genotypeoutname: AIMs.geno
snpoutname:      AIMs.snp
indivoutname:    AIMs.ind
familynames:     NO

smartpca.perl -i AIMs.geno -a AIMs.snp -b AIMs.ind -o AIMs.pca -e AIMs.eval -l AIMs.log -p AIMs -m 0

smartpca -p AIMs.pca.par >AIMs.log
ploteig -i AIMs.pca.evec -c 1:2  -p Control  -x  -y  -o AIMs.xtxt
evec2pca.perl 10 AIMs.pca.evec AIMs.ind AIMs.pca

wget https://raw.githubusercontent.com/Shicheng-Guo/ExonChipProcessing/master/PCAPlot.R
Rscript PCAPlot.R AIMs.pca

awk '$1==23{print $2}' RA2020.1.bim > chr23_26.txt
awk '$1==24{print $2}' RA2020.1.bim >> chr23_26.txt
awk '$1==25{print $2}' RA2020.1.bim >> chr23_26.txt
awk '$1==26{print $2}' RA2020.1.bim >> chr23_26.txt

plink --bfile RA2020 --maf 0.1 --exclude chr23_26.txt --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile RA2020 --extract indepSNP.prune.in --genome --out RAindepSNP



ind<-read.table("RAindepSNP.genome",head=T,sep="")
ind1<-subset(ind,PI_HAT>0.25)
write.table(ind1,file="Relatives.RA2020.txt",col.names=NA,row.names=T,sep="\t",quote=F)
pdf("PI_HAT.pdf")
hist(ind1$PI_HAT,xlim=c(0,1),col="blue",breaks = seq(0,1,by=0.001),border="blue",ylim=c(0,250),xlab="PI_HAT, width=0.001",main="")
dev.off()


plink --bfile RA2020 --mind 0.05 --make-bed --out RA2020-B1
plink --bfile RA2020-B1 --geno 0.95 --make-bed --out RA2020-B2
plink --bfile RA2020-B2 --maf 0.01 --make-bed --out RA2020-B3
plink --bfile RA2020-B3 --hwe 0.00001 --make-bed --out RA2020-B4
plink2 --bfile RA2020-B4 --king-cutoff 0.125
plink2 --bfile RA2020-B4 --remove plink2.king.cutoff.out.id --make-bed -out RA2020-B5
plink --bfile RA2020-B5 --check-sex
plink --bfile RA2020-B5 --impute-sex --make-bed --out RA2020-B6
plink --bfile RA2020-B6 --check-sex
grep PROBLEM plink.sexcheck | awk '{print $1,$2}' > sexcheck.remove
plink --bfile RA2020-B6 --remove sexcheck.remove --make-bed --out RA2020-B7
plink --bfile RA2020-B7 --test-missing midp 
awk '$5<0.000001{print}' plink.missing | awk '{print $2}' > missing.imblance.remove
plink --bfile RA2020-B7 --exclude missing.imblance.remove --make-bed --out RA2020-B8
plink --bfile RA2020-B8 --assoc mperm=5000 --adjust gc --threads 31
plink --bfile RA2020-B8 --pca --threads 31
plink --bfile RA2020-B8 --logistic --covar plink.eigenvec --covar-number 1-20 --adjust

plink --bfile hapmap1 --assoc --adjust --out as2


fastq-dump -X 1000  -I --split-files SRR390728
fastq-dump -X 1000 SRR390728



p<-read.table("plink.assoc",head=T,sep="\t")
P<-p$P
library("Haplin")
png("qqplot.ra.png")
pQQ(na.omit(P), nlabs =length(na.omit(P)), conf = 0.95)
dev.off()

p<-read.table("plink.assoc.mperm",head=T)
P<-p$EMP1
png("qqplot.ra.emp1.png")
pQQ(na.omit(P), nlabs =length(na.omit(P)), conf = 0.95)
dev.off()

p<-read.table("plink.assoc.logistic",head=T)
p<-subset(p,TEST=="ADD")
P<-p$P
png("qqplot.ra.glm.png")
pQQ(na.omit(P), nlabs =length(na.omit(P)), conf = 0.95)
dev.off()


gatk LiftoverVcf -I GCF_000001405.38.dbSNP152.GRCh38p12b.GATK.vcf.gz -O dbSNP153.hg19.vcf -C hg38ToHg19.over.chain.gz --REJECT rejected.vcf -R ~/hpc/db/hg19/hg19.fa
gatk LiftoverVcf -I RA2020.dbSNP.hg19.vcf.gz -O RA2020.dbSNP.hg38.vcf.gz -C hg19ToHg38.over.chain.gz --REJECT rejected.vcf -R ~/hpc/db/hg38/hg38.fa

cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R
mkdir temp
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/contigReplace.pl
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo perl contigReplace.pl $i \| bgzip -c \> chr$i.dose.contig.vcf.gz >>$i.job
echo tabix -p vcf chr$i.dose.contig.vcf.gz >>$i.job
qsub $i.job
done


cd /gpfs/home/guosa/hpc/db/hg19/dbSNP152
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view dbSNP152.hg19.chr$i.vcf.recode.vcf -Oz -o dbSNP152.chr$i.hg19.vcf.gz >>$i.job
echo tabix -p vcf dbSNP152.chr$i.hg19.vcf.gz >>$i.job
qsub $i.job
done


cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R
mkdir temp
for i in 1
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools annotate -a ~/hpc/db/hg19/dbSNP152/dbSNP152.chr$i.hg19.vcf.gz -c ID  chr$i.dose.contig.vcf.gz -Oz -o chr$i.dose.dbSNP.hg19.vcf.gz >>$i.job
echo bcftools view -i \'\(IMP=1 \& R2\>0.6\)\|IMPUTED=0\' chr$i.dose.dbSNP.hg19.vcf.gz \|  bcftools annotate -x \^FORMAT/GT -Oz -o chr$i.dose.dbSNP.clean.hg19.vcf.gz  >>$i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -f -p vcf chr$i.dose.dbSNP.hg19.vcf.gz >> $i.job
echo bcftools view -i \'R2\>0.6\|TYPED=1\|TYPED_ONLY=1\' -R MUC.hg19.sort.bed chr$i.dose.dbSNP.hg19.vcf.gz \|  bcftools annotate -x \^FORMAT/GT -Oz -o chr$i.dose.MUC.clean.hg19.vcf.gz  >>$i.job
qsub $i.job
done


cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -f -p vcf chr$i.dose.dbSNP.hg19.vcf.gz >> $i.job
echo bcftools view -i \'R2\>0.6\|TYPED=1\|TYPED_ONLY=1\' chr$i.dose.dbSNP.hg19.vcf.gz \|  bcftools annotate -x \^FORMAT/GT -Oz -o chr$i.dose.clean.hg19.vcf.gz  >>$i.job
qsub $i.job
done



bcftools annotate -a ~/hpc/db/hg19/ -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') 

bcftools view -i '%iD=="rs35705950"' MUC.anno.hg19.vcf.gz | less -S 
bcftools view -i '%iD=="rs79920422"' MUC.anno.hg19.vcf.gz | less -S 

 
 
 
perl -p -i -e '{s/chr//g}' MUC

/gpfs/home/guosa/hpc/db/hg19/dbSNP152

mkdir ./temp/chr1/
mkdir ./temp/chr7/
mkdir ./temp/chr8/
mkdir ./temp/chr9/
zcat dbSNP152.chr1.hg19.vcf.gz | vcf-sort -p 16 -t ./temp/ | bgzip -c > dbSNP152.chr1.hg19.sort.vcf.gz &
zcat dbSNP152.chr7.hg19.vcf.gz | vcf-sort -p 16 -t ./temp/ | bgzip -c > dbSNP152.chr7.hg19.sort.vcf.gz &
zcat dbSNP152.chr8.hg19.vcf.gz | vcf-sort -p 16 -t ./temp/chr8 | bgzip -c > dbSNP152.chr8.hg19.sort.vcf.gz &
zcat dbSNP152.chr9.hg19.vcf.gz | vcf-sort -p 16 -t ./temp/chr9 | bgzip -c > dbSNP152.chr9.hg19.sort.vcf.gz &
mv dbSNP152.chr8.hg19.sort.vcf.gz dbSNP152.chr8.hg19.vcf.gz
mv dbSNP152.chr9.hg19.sort.vcf.gz dbSNP152.chr9.hg19.vcf.gz
tabix -p vcf dbSNP152.chr8.hg19.vcf.gz &
tabix -p vcf dbSNP152.chr9.hg19.vcf.gz &


cd /gpfs/home/guosa/hpc/db/hg19/dbSNP152

mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -f -p vcf dbSNP152.chr$i.hg19.vcf.gz >>$i.job
qsub $i.job
done




gatk VariantAnnotator -R ~/hpc/db/hg19/hg19.fa -V chr22.dose.contig.vcf.gz -O chr22.dose.contig.dbSNP.vcf.gz --dbsnp ~/hpc/rheumatology/RA/he2019/imphase/dbSNP152.GSC.hg19.vcf

java -jar GenomeAnalysisTK.jar -T VariantAnnotator 

--help

java -jar GenomeAnalysisTK.jar -R ~/hpc/db/hg19/hg19.fa -T VariantAnnotator -L chr22.dose.contig.vcf.gz -o chr22.dose.contig.dbSNP.vcf.gz
   
 java -jar ~/hpc/tools/SnpSift.jar annotate ~/hpc/rheumatology/RA/he2019/imphase/dbSNP152.GSC.hg19.vcf  chr22.dose.contig.vcf.gz > chr22.dose.contig.dbsnp.vcf.gz

snpsift annotate ~/hpc/rheumatology/RA/he2019/imphase/dbSNP152.GSC.hg19.vcf chr22.dose.contig.vcf.gz

zcat chr22.dose.contig.vcf.gz | head -n 500 | bgzip -c > x.vcf.gz

snpsift annotate -id ~/hpc/rheumatology/RA/he2019/imphase/dbSNP152.GSC.hg19.vcf x.vcf.gz

 alias snpsift="java -jar ~/hpc/tools/snpEff/SnpSift.jar"

 
/gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3/test

  cp /gpfs/home/guosa/hpc/rheumatology/RA/he2019/imphase/dbSNP152.GSC.hg19.vcf ~/hpc/db/hg19 &
  

To change the SNP name, you need a list with all SNPs available.
1) Go to UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables) and download a list of all available SNPs with chr, pos and rs number.
2) Create a file with chr:pos in column 1 and rs number in column 2.
3) In PLINK, use options --update-map and --update-name to change your SNPs name. WARNING: Some position have more than one SNP.




cd /gpfs/home/guosa/hpc/project/pmrp/merge
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile PMRP-Phase1-phase2-Full --chr $i --recode vcf --out PMRP2019.chr$i >>$i.job
echo bcftools view PMRP2019.chr$i.vcf -Oz -o PMRP2019.chr$i.vcf.gz >>$i.job
qsub $i.job
done

awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' hg19.fa.fai 



cd /gpfs/home/guosa/hpc/rheumatology/RA/he2019/imphase
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools reheader -h chr$i.dose.filter.newheader chr$i.dose.filter.vcf.gz -o chr$i.dose.filter.newheader.vcf >>$i.job
echo bcftools view -G chr$i.dose.filter.newheader.vcf -Oz -o chr$i.dose.filter.newheader.vcf>>$i.job
qsub $i.job
done



cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo unzip -P aaViTz8M6tkZAU chr_$i.zip >>$i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo unzip -P \'dOBdaLMv\>H49s\' chr_$i.zip >>$i.job
qsub $i.job
done
  
  
 
  
plink --vcf chr2.dose.vcf.gz 'dosage='HDS 

https://www.dropbox.com/s/qv61mgtx6pz54fz/chr1_phase3.pgen.zst?dl=1

bgzip -c 

  
open F,"hg19.contig.txt" || die print "hg19.contig.txt";
my %len;
while(<F>){
chomp;
my ($chr,$len)=split/\s+/;
$len{$chr}=$len;
}
close F;

foreach my $chr(1..22){
open F1,"chr$i.dose.filter.header" || die "Could not open file 'chr$i.dose.filter.header'. $!";
while(<F1>){
print "$_";
}

for i in {1..22}
do 
bcftools reheader -h chr$i.dose.filter.newheader chr$i.dose.filter.vcf.gz 
done

  
for i in {1..22}
do
bcftools view -h chr$i.dose.filter.vcf.gz >  chr$i.dose.filter.header
done  
  
open F,"hg19.contig.txt";
my %len;
while(<F>){
chomp;
my ($chr,$len)=split/\s+/;
$len{$chr}="##contig=<ID=$chr,length=$len>";
}

foreach my $chr(1..22){
print $chr;
}


cd /gpfs/home/guosa/hpc/rheumatology/RA/he2019/impute
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # bcftools view chr$i.phased.vcf.gz -R MUC.hg19.bed -Ov -o MUC.chr$i.vcf
echo bcftools view -i '(IMP=1 & R2>0.6)|IMPUTED=0' All_samples_Exome_QC.chr22.vcf.gz -Oz -o test.vcf.gz
tabix -p vcf test.vcf.gz



for i in {1..22};
  do
	echo \#PBS -N $i  > $i.job
	echo \#PBS -l nodes=1:ppn=1 >> $i.job
	echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
	echo \#PBS -m abe  >> $i.job
	echo \#PBS -o $(pwd)/temp/ >>$i.job
	echo \#PBS -e $(pwd)/temp/ >>$i.job
	echo cd $(pwd) >> $i.job
	echo gunzip chr$i.info.gz >>$i.job
	echo bcftools view -i \'R2\>0.6\|TYPED=1\|TYPED_ONLY=1\' -Oz chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.vcf.gz >>$i.job
	echo plink --vcf chr$i.dose.filter.vcf.gz --make-bed --out chr$i.s1 >>$i.job
	echo plink --bfile chr$i.s1 --list-duplicate-vars --out chr$i >>$i.job
	echo plink --bfile chr$i.s1 --exclude plink.dupvar --make-bed --out ../plinkout/chr$i >> $i.job
	qsub $i.job
done




wget https://imputationserver.sph.umich.edu/share/results/bef4440e3b1075b1e36132516faa3ef3/chr_1.zip
wget https://imputationserver.sph.umich.edu/share/results/1336b8967d31c223a569adf29412caa4/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/91c6e45418f942d69073e281cf3c5904/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/35920d77b5210b944ab61e218e73972c/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/64faf60dec19e796a254b8e5e4376e3d/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/7b6243c571ed0296769241a47ae74ed/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/99cd67c29c8b3fe19ff5f18806909d71/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/95d6d2dc584325fa929b999a6354a0e1/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/d8fa74efe670d6c065df3e80f78b282f/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/6f8fd584e0811bcc63986414a4c41328/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/21fe31b051d2db6af92033b118d5cbe/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/c30b986014f813fbdaab8be6264803e0/chr_2.zip
wget https://imputationserver.sph.umich.edu/share/results/e213c4bd06a8f7e8d7fc2afaf18b56a9/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/95908a009e3e59bf7ea642f65fb920c7/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/f0eae19cf181a38ecf9938733b2c08d7/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/d6645b933656c261046390e6bb5f9696/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/a9221a9b2ec5b5d6f1b791e6e9450cc5/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/bd4f613061782d49115de51570019fe1/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/e5681641a2ed8f2010dd1c2a767d8d2f/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/468e27dea2d50fec15d06440a6ade64d/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/f6b602737cd162c836b94bc52f51503e/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/d6231c89c5a760226aab316cc3b14ddd/chr_9.zip

for i in {1..22}
do
unzip -P 0S8RedytR%:iDH chr_$i.zip 
done


for i in {1..22}
do
tabix -p vcf chr$i.dose.vcf.gz &
done

for i in {1..22}
do
echo $i chr$i >>chr_name_conv.txt
done

for i in {1..22};
  do
	echo \#PBS -N $i  > $i.job
	echo \#PBS -l nodes=1:ppn=1 >> $i.job
	echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
	echo \#PBS -m abe  >> $i.job
	echo \#PBS -o $(pwd)/temp/ >>$i.job
	echo \#PBS -e $(pwd)/temp/ >>$i.job
	echo cd $(pwd) >> $i.job
	echo bcftools annotate --rename-chrs chr_name_conv.txt chr$i.dose.vcf.gz \| bgzip \> chr$i.dose.chr.vcf.gz  >> $i.job
	qsub $i.job
done

  cd ~/hpc/rheumatology/RA/he2019/impute
  cp ~/hpc/db/hg19/hg19.fa ./
  perl -p -i -e 's/chr//g' hg19.fa
  gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict 
  cat hg19.dict 
  samtools faidx hg19.fa
  cat hg19.fa.fai 
  
  dbsnp152="/home/guosa/hpc/db/hg19/dbSNP152.hg19.vcf"
  gatk VariantAnnotator -R hg19.fa -V chr20.dose.vcf.gz -O chr20.dose.rs.vcf.gz --dbsnp dbSNP152.GSC.hg19.vcf

  awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' hg19.fa.fai > contig.len.hg19.txt
  
  awk '/^#CHROM/ { printf("##contig=<ID=1,length=195471971>\n##contig=<ID=2,length=182113224>\n");} {print;}' in.vcf > out.vcf

  gatk SelectVariants -V chr22.dose.vcf.gz -R hg19.fa -O chr22.dose.contig.vcf.gz

  
gatk FixVcfHeader -I chr22.dose.vcf.gz -O chr22.dose.contig.vcf.gz
	 
plink --file RA56 --merge RA57.ped RA57.map --recode --make-bed --out RA113
plink --bfile RA113 --bmerge he2019.bed he2019.bim he2019.fam --recode --make-bed --out RA2020


plink --bfile he2019 --bmerge RA113.bed RA113.bim RA113.fam --recode --make-bed --out RA2020


plink --bfile ../RA2020 --bmerge hapmap3_r1_b37_fwd_consensus.qc.poly.recode.r3 --make-bed --geno 0.1 --maf 0.05 --hwe 0.00001 --out RA2204-Hapmap


plink --bfile ~/hpc/db/Hapmap/hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode --exclude RA2204-Hapmap-merge.missnp --make-bed --out hapmap3_r1_b37_fwd_consensus.qc.poly.recode.r3





awk '/^#CHROM/ { printf("##contig=<ID=1,length=195471971>\n##contig=<ID=2,length=182113224>\n");} {print;}' in.vcf > out.vcf

awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai
awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' hg19.fa.fai

bcftools view RA113.chr1.vcf -Oz -o RA113.chr1.vcf.gz

cd ~/hpc/rheumatology/RA/he2019/plink/vcf
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view RA113.chr$i.vcf -Oz -o RA113.chr$i.vcf.gz >>$i.job
echo tabix -f -p vcf RA113.chr$i.vcf.gz >> $i.job
qsub $i.job
done

cd ~/hpc/rheumatology/RA/he2019/imphase
mkdir temp
for i in {1..24}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo gunzip chr$i.dose.vcf.gz >>$i.job
qsub $i.job
done


  bcftools reheader -f hg19.fa.fai chr22.dose.vcf.gz -Oz -o chr22.dose.contig.vcf.gz
  fetchChromSizes hg38 > hg38.chrom.sizes
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz -O refGene.hg38.txt.gz
  
  
  perl ref2tss.pl -refGene refGene.hg38.txt -output refGene.hg38 -chromsize hg38.chrom.sizes
  
  
  
  for i in {1..22}
  do
  echo NC_000001.$i $i >> chr_name_conv.txt
  done
  
  vcftools --vcf dbSNP152.GSC.hg19.vcf --plink --out dbSNP152.hg19

  perl addfakegenotype.pl > dbSNP152.GSC.hg19.vcf
  vcftools --vcf dbSNP152.GSC.hg19.vcf --plink --out dbSNP152.hg19

  use strict;
  open F,"/home/guosa/hpc/db/hg19/dbSNP152.hg19.vcf";
  while(<F>){
  chomp;
  print "$_\n" if /^#/;
  my $line=$_;
  my $line=~s/NC_000001//;
  print "$line\n";
  }
  
    
  cd /home/guosa/hpc/db/hg19/
  gatk CreateSequenceDictionary -R hg19.fa -O /home/guosa/hpc/db/hg19/hg19.dict
  cat hg19.dict 
  samtools faidx hg19.fa 


#!/usr/bin/perl
use Cwd;
my $FF1="/thinker/aid/udata/bing/beijingYixian/fastqs/Methy";
my $FF2="/thinker/storage/udata/muyl/genomes/human/ucsc.hg19.nb/bowtie2Index";
my $FF3="/thinker/storage/udata/muyl/genomes/human/ucsc.hg19.nb/ucsc.hg19.fasta";
open F,"saminfo.txt";
while(<F>){
chomp;
my ($file1,$file2)=split/\s+/;
my ($sam)=split/_R1/,$file1;
open F2,">$sam.sh";
print F2 "bowtie2 -q -x $FF2/hg19 -1 $FF1/$file1 -2 $FF1/$file2 -S $sam.sam\n"; 
print F2 "samtools view -bS $sam.sam > $sam.bam\n";
print F2 "samtools sort $sam.bam $sam.sorted\n";
print F2 "samtools mpileup -uf $FF3 $sam.sorted.bam -v -o $sam.vcf\n";
}



scp /home/nu_guos/.ssh/id_rsa.pub 
/home/biostaff/.ssh
scp -P 37122 biostaff@ada.zettadom.com:/home/biostaff/COSCE129LIN64.bin ./
ssh -p 37122 biostaff@ada.zettadom.com
cat /home/nu_guos/.ssh/id_rsa.pub | ssh -p 37122 biostaff@ada.zettadom.com 'cat >> .ssh/authorized_keys'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh -p 37122 biostaff@ada.zettadom.com 'chmod 700 ~/.ssh'
alias ada="ssh -p 37122 biostaff@ada.zettadom.com"
alias host="ssh -p 8809 guosa@localhost"
python get-pip.py --user
ln -s /thinker/storage/udata/guosa hpc

condor_q -hold
ssh-keygen -t rsa
ssh shg047@23.99.137.107 'mkdir -p .ssh'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'cat >> .ssh/authorized_keys'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'chmod 700 ~/.ssh'

cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'cat >> .ssh/authorized_keys2'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'chmod 700 ~/.ssh/authorized_keys2'
ssh shg047@23.99.137.107
alias chtc="ssh nu_guos@submit-3.chtc.wisc.edu"

plink --bfile he2019 --assoc --count --maf 0.05 --allow-no-sex --1

plink<-read.table("plink.assoc",head=T)
ManhattanPlot<-function(mylimma){
library(qqman)
res <- mylimma
SNP=res$SNP
CHR=res$CHR
if(length(grep("X",CHR))>0){
  CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
  CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
}
CHR<-as.numeric(CHR)
BP=res$BP
P=res$P
manhattaninput=na.omit(data.frame(SNP,CHR,BP,P))
genomewideline=0.05/nrow(manhattaninput)
pdf("assoc.manhattan.pdf")
manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,30),genomewideline=F,lwd=2, suggestiveline=F)
dev.off()
}
ManhattanPlot(plink)
plink$logP<--log(plink$P,10)

sigplink<-subset(plink, logP>10)

plink --bfile he2019 --logistic --covar plink.eigenvec --covar-number 1-6 --allow-no-sex --1
ManhattanPlot<-function(mylimma){
library(qqman)
res <- mylimma
SNP=res$SNP
CHR=res$CHR
if(length(grep("X",CHR))>0){
  CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
  CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
}
CHR<-as.numeric(CHR)
BP=res$BP
P=res$P
manhattaninput=na.omit(data.frame(SNP,CHR,BP,P))
genomewideline=0.05/nrow(manhattaninput)
head(manhattaninput)
dim(manhattaninput)
pdf("manhattan.pdf")
manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,20),genomewideline=-log10(genomewideline),lwd=2, suggestiveline=F)
dev.off()
}
plink<-read.table("plink.assoc.logistic",head=T)
mylimma<-subset(plink,TEST=="ADD")
ManhattanPlot(plink)
/home/guosa/hpc/rheumatology/RA/he2019/sigplink.bed

for i in {1..26}
do
plink --bfile he2019 --chr $i --recode vcf --out ./vcf/he2009.chr$i
done

cd /gpfs/home/guosa/hpc/rheumatology/RA/he2019/impute
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -f -p vcf chr$i.phased.vcf.gz >> $i.job
echo bcftools view chr19.phased.vcf.gz -R MUC.hg19.bed -Ov -o MUC.chr19.vcf
qsub $i.job
done




cd /gpfs/home/guosa/hpc/rheumatology/RA/he2019/impute
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo wget https://imputationserver.sph.umich.edu/share/results/f3bd4c751822ad5b595f84ddce8633dd/chr_$i.zip  >>$i.job
qsub $i.job
done

wget https://imputationserver.sph.umich.edu/share/results/f3bd4c751822ad5b595f84ddce8633dd/chr_1.zip
wget https://imputationserver.sph.umich.edu/share/results/4a40083b89985685aa7497410b76de2e/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/bed7f5623a79df56011c415de914bea8/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/578ff3fa346da5c04804abdb0eb93900/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/1ad44aa0169dacd18d2328a4d08cd9f5/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/1182c8c1af253aa56b9b5b9002a2575a/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/86825c0065032220e8d74c811396debe/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/fb068bed6ec28f6c19c08c1f62ac508c/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/f1d226ea4115c6c5c05f84bb6a7fedb8/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/54bdee5c147c31f77b2ce393eb1be295/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/220100a69422665935a46cd3733112b0/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/2c20372088ebc27233341ab964641baa/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/bfca52b60d40dbe2965464e9b93dcb1c/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/beabe65757681d5f4e5f5cb068a4e58b/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/4593e713bc69c4072b380bc7a4c037a9/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/e42511c91c695d1153ac5a48c32231e2/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/26565945d9d95c950f0db8c9afb11f3c/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/c5c56e404c061669b4e26d0bdefdbd49/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/8f6da5b2dda3e67f3933170ae23b10bd/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/39678907b07f5517f6af388ff1d962a1/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/3c6ed7d7031cf5e28ef9b0cf9b6d3992/chr_9.zip

for i in {1..22}
do
unzip -P gK?9sQr5bTtJZR chr_$i.zip 
done

for i in {1..22}
do
unzip -P ci7PvDMsT8zqsA chr_$i.zip 
done


cd /gpfs/home/guosa/hpc/rheumatology/SLE/BCR/vdj
mkdir temp
for i in $(cat fq.txt)
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mixcr align -f -s hsa -p rna-seq $i\_R1.fq $i\_R2.fq $i.vdjca -r $i  >>$i.job
echo mixcr assemble -f --write-alignments $i.vdjca $i.clna >>$i.job
echo mixcr assembleContigs -f $i.clna $i.clns >>$i.job
echo mixcr exportClones -f $i.clns $i.txt >>$i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/rheumatology/SLE/BCR/vdj
mkdir temp
for i in $(cat fq.txt)
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bbmerge.sh in1=$i\_R1.fq in2=$i\_R2.fq out=../imgt/$i.fq outu1=../imgt/$i.u1map outu2=../imgt/$i.u2map >>$i.job
qsub $i.job
done

for i in `ls *.fq`
do
sed -n '1~4s/^@/>/p;2~4p' $i > $i.fasta
echo $i
done

for i in `ls *.fasta`
do
seqtk sample -s seed=110 $i 150000 > $i.imgt
echo $i
done

for i in $(cat fq.txt)
do
sed -e 's/\s,\+/\t/g' -e 's/,\+\s/\t/g' $i.txt | sed -e 's/;,\+/;/g' > $i.sed
done

vdjtools Convert -S mixcr -m vdjtools.m.txt  metadata.txt
vdjtools CalcSegmentUsage -p -m metadata.txt vdjtools

vdjtools CalcSegmentUsage -f CellType -p -m metadata.txt vdjtools.celltype
vdjtools CalcSegmentUsage -f SampeName -p -m metadata.txt vdjtools.samplename
vdjtools CalcSegmentUsage -f SampeType -p -m metadata.txt vdjtools.sampletype

vdjtools CalcSpectratype -a -m metadata.txt vdjtools
vdjtools CalcPairwiseDistances -p -m metadata.txt vdjtools

for i in `ls metadata.*.txt`
do
vdjtools CalcPairwiseDistances -p -m $i vdj.$i
done

vdjtools CalcPairwiseDistances -p -m metadata.S.txt vdj.S
vdjtools ClusterSamples -p vdj.S  vdj.S
vdjtools CalcPairwiseDistances -p -m metadata.H.txt vdj.H
vdjtools ClusterSamples -p vdj.H  vdj.H

for i in `ls metadata.HC*.txt`
do
vdjtools CalcPairwiseDistances -p -m $i  $i.circos 
done

for i in `ls metadata.S*.txt`
do
vdjtools CalcPairwiseDistances -p -m $i  $i.circos 
done



sed -e 's/\s,+/\t/g' -e 's/,+\s/\t/g' 16S1710LF01.txt | sed -e 's/;,+/;/g' >16S1710LF01.txt.sed.txt
sed -e 's/\s,\+/\t/g' -e 's/,\+\s/\t/g' 16S1710LF01.txt | sed -e 's/;,\+/;/g' > 16S1710LF01.txt.sed.txt

vdjtools PlotFancySpectratype [options] sample.txt output_prefix
vdjtools Convert -S mixcr -m vdjtools.txt  metadata.txt
qqplot(P)
vdjtools PlotFancyVJUsage vdjtools.segments.wt.V.txt 16S1710LF01
vdjtools CalcPairwiseDistances -p -m metadata.txt vdjtools

vdjtools Convert -S mixcr 16S1710LF01.txt.sed.txt  16S1710LF01.vdjtools













bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>
for i in `ls *extendedFrags.fastq`
do
sed -n '1~4s/^@/>/p;2~4p'  $i > $i.fa
done


cp 1st/* vdj/
cp 2nd/* vdj/
gunzip *.gz 
http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip
unzip hisat2-2.1.0-source.zip
cd /gpfs/home/guosa/hpc/tools/hisat2-2.1.0
make


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-2")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
data[1:12,1:12]
methdata<-data.matrix(data[,c(12:180)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rownames(methdata))
head(colnames(methdata))

cat 180205LSJ32.fq.extendedFrags.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $i.fa


ouggnehcihs

https://shicheng-guo.github.io/assets/images/Bewerbungsfoto.jpg
https://shicheng-guo.github.io/assets/images/my_story.png
https://shicheng-guo.github.io/machine_learning/2017/dashboard_screenshot.png
https://shicheng-guo.github.io/
https://raw.githubusercontent.com/ShirinG/ShirinG.github.io/master/assets/images/logo.png
avatar: https://raw.githubusercontent.com/ShirinG/ShirinG.github.io/master/assets/images/logo.png
https://shiring.github.io/assets/images/doublehelix.png
https://shicheng-guo.github.io/


# 2019-07-03
cd ~/hpc/tools
git clone https://github.com/carjed/helmsman.git
cd helmsman
pip install virtualenv --user
 
mv S001B_Tumor2-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf   001B_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
mv S001D_Tumor2_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf  001D_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
mv S001L_Tumor2-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001L_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
mv S001N_Tumor2_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001N_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf

mv S001B_Tumor1-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf    001B_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf
mv S001D_Tumor1_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf  001D_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf
mv S001L_Tumor1-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001L_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf
mv S001N_Tumor1_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001N_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf

ls -l *Tumor2*.vcf | wc -l
ls -l *Tumor1*.vcf | wc -l

cd /home/guosa/hpc/project/LungBrainMetastasis
awk '{print $1,$2,$4,$5}' OFS="\t" *Normal-Tumor1_mem*.bed  | sort -u > Normal-Tumor1.uni.hg19.bed
awk '{print $1,$2,$4,$5}' OFS="\t" *Normal-Tumor2_mem*.bed  | sort -u > Normal-Tumor2.uni.hg19.bed


## step 1.1 check sample names and remove or keeps samples of interest. 
rm Tumor1_VarD.txt
cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
bcftools query -l 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf | grep Tumor1 | grep VarD >> Tumor1_VarD.txt
done

rm Tumor2_VarD.txt
cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
bcftools query -l 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf | grep Tumor2 | grep VarD >> Tumor2_VarD.txt
done

## Step 2. copy raw vcf to new folder to keep raw data safe
cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
cp 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf $i.T1.vcf
done

cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
cp 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf $i.T2.vcf
done

## step 3. 
for i in A B C E F G H I J K L M O
do
echo $i
bcftools view --force-samples -S Tumor1_VarD.txt $i.T1.vcf | bcftools annotate -x ID,^INFO/AN,INFO/DP,FORMAT -I +'%CHROM:%POS' | bcftools sort -Ov -o $i.T1.sort.vcf
perl -p -i -e 's/\.\/\./0\/1/g' $i.T1.sort.vcf
bgzip -f $i.T1.sort.vcf
tabix -f -p vcf $i.T1.sort.vcf.gz
done
ls *.T1.sort.vcf.gz > T1.txt
bcftools merge -l T1.txt -Ov | bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') > T1.vcf


cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
echo $i
bcftools view --force-samples -S Tumor2_VarD.txt $i.T2.vcf | bcftools annotate -x ID,^INFO/AN,INFO/DP,FORMAT -I +'%CHROM:%POS' | bcftools sort -Ov -o $i.T2.sort.vcf
perl -p -i -e 's/\.\/\./0\/1/g' $i.T2.sort.vcf
bgzip -f $i.T2.sort.vcf
tabix -f -p vcf $i.T2.sort.vcf.gz
done
ls *.T2.sort.vcf.gz > T2.txt
bcftools merge -l T2.txt -Ov | bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') > T2.vcf

touch -m T*.vcf

## step 4. Rename chrosome and re-annoate VCF with ANNOVAR

rm rename-chrs.txt
for i in {1..24} X Y
do
echo -e chr$i'\t'$i >> rename-chrs.txt 
done 

bcftools annotate --rename-chrs rename-chrs.txt T1.vcf -Ov -o T1.chr.vcf
bcftools annotate --rename-chrs rename-chrs.txt T2.vcf -Ov -o T2.chr.vcf
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') T1.chr.vcf 
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') T2.chr.vcf 

table_annovar.pl -vcfinput T1.chr.vcf ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 -out T1 -remove -protocol refGene,dbnsfp33a -operation gx,f -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt
table_annovar.pl -vcfinput T2.chr.vcf ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 -out T2 -remove -protocol refGene,dbnsfp33a -operation gx,f -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt
 
## step 4. apply maftools do analysis. 



# 2019-06-26
cd  /gpfs/home/guosa/hpc/db/hg19/beagle
/gpfs/home/guosa/hpc/db/hg19/cpgSNP.hg19.bed

cp /gpfs/home/guosa/hpc/db/hg19/cpgSNP/commonSNP/*cpgSNP.bed ./
for i in {1..22}
do
perl -p -i -e 's/chr//' chr$i.cpgSNP.bed
done
perl -p -i -e 's/chrX/23/' chrX.cpgSNP.bed
perl -p -i -e 's/chrY/24/' chrY.cpgSNP.bed

for i in `ls *.bed`
do
bedtools sort -i $i > $i.sort.bed
bgzip $i.sort.bed
tabix -p bed  $i.sort.bed.gz
done

for i in {1..22} X Y
do
perl -lane '{print $_ if $_=~/\t[A-Z]\/[A-Z]\t/}' chr$i.cpgSNP.bed.sort.bed > chr$i.cpgSNP.bin.bed 
done

cd /gpfs/home/guosa/hpc/db/hg19/beagle/cgSNP/vcf
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -T ../chr$i.cpgSNP.bin.bed ../../chr$i.1kg.phase3.v5a.vcf.gz -Oz -o chr$i.cpgSNP.vcf.gz >>$i.job
qsub $i.job
done



#### CpG-SNPs in human genomes
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/beagle/cgSNP")
head(data)
for(GAP in seq(5,20,1)){
rlt<-c()
for(j in c(1:22,"X","Y")){
data<-read.table(paste("chr",j,".cpgSNP.bin.bed",sep=""))
i=1
print(paste("chr",j,".cpgSNP.bin.bed",sep=""))
while(i<(nrow(data)-100)){
  if(table(data[i:(i+GAP),6])[1]>GAP){
   rlt<-rbind(rlt,c(data[i,1],data[i,2],data[i+GAP,3]))
   print(i)
   i=i+GAP
  }
  i=i+1
}
}
write.table(rlt,file=paste("Long_CpG_Gain.GAP",GAP,"hg19.txt",sep="."),sep="\t",col.names = F,row.names = F)
}


# 2019-06-25
mkdir ~/hpc/tools/bcftools-1.9/db
cp ~/hpc/db/hg19/hg19.fa ~/hpc/tools/bcftools-1.9/db
perl -p -i -e 's/>chr/>/' ~/hpc/tools/bcftools-1.9/db/hg19.fa
perl -p -i -e 's/>X/>23/' ~/hpc/tools/bcftools-1.9/db/hg19.fa
perl -p -i -e 's/>Y/>24/' ~/hpc/tools/bcftools-1.9/db/hg19.fa
samtools faidx ~/hpc/tools/bcftools-1.9/db/hg19.fa
# REF/ALT total/modified/added:   2601741/397967/114
bcftools norm -t "^24,25,26" -m-any --check-ref s -f ~/hpc/tools/bcftools-1.9/db/hg19.fa All_samples_Exome_QC.clean.norm.vcf.gz -Ov | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS' -Oz -o All_samples_Exome_QC.clean.vcf.gz

mkdir chr
mkdir temp
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar -O beagle.16May19.351.jar
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar

inputvcf="All_samples_Exome_QC.clean.vcf.gz"
tabix -f -p vcf All_samples_Exome_QC.clean.vcf.gz
for i in {1..23}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=24 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mkdir ./temp/chr$i >> $i.job
echo bcftools view -e \'ALT =\"-\" \| REF =\"-\"\' -t $i $inputvcf -Oz -o ./chr/All_samples_Exome_QC.clean.chr$i.vcf.gz >>$i.job
echo java -Djava.io.tmpdir=./temp/chr$i -Xmx64g -jar beagle.16May19.351.jar gt=./chr/All_samples_Exome_QC.clean.chr$i.vcf.gz ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz map=~/hpc/db/hg19/beagle/plink.chr$i.GRCh37.map out=./chr/All_samples_Exome_QC.chr$i.vcf >>$i.job
qsub $i.job
done




HRC (Version r1.1 2016)
This HRC panel consists of 64,940 haplotypes of predominantly European ancestry.

12:163819-191525

Window 2 (16:17024548-57696735)
Reference markers:     291,891
Study markers:          17,271
##contig=<ID=16,length=90274696>
q

ERROR: missing REF or ALT allele at 22:17426007
bcftools view -t "22:17619023" ./chr/All_samples_Exome_QC.clean.chr22.vcf.gz | less -S

bcftools view -t "22:17619023" All_samples_Exome_QC.clean.norm.vcf.gz | less -S


./chr/All_samples_Exome_QC.temp.vcf.recode.clean.chr13.vcf.recode.vcf

java -jar ./conform-gt.24May16.cee.jar gt=All_samples_Exome_QC.clean.vcf.gz ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz  out=All_samples_Exome_QC.clean.gtc.vcf.gz

./chr/All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf

All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf

mkdir ./temp/chr22
bcftools view -t 22 All_samples_Exome_QC.clean.vcf.gz -Oz -o ./chr/All_samples_Exome_QC.clean.chr22.vcf.gz
java -Djava.io.tmpdir=./temp/chr22 -Xmx32g -jar beagle.16May19.351.jar gt=./chr/All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf ref=/gpfs/home/guosa/hpc/db/hg19/beagle/EUR/chr22.1kg.phase3.v5a.EUR.vcf.gz map=/gpfs/home/guosa/hpc/db/hg19/beagle/plink.chr22.GRCh37.map out=All_samples_Exome_QC.chr22.vcf



[PolyPhen-2]()
[LS-SNP/PDB]()
[SNPeffect 4.0]()
[Missense3D](http://www.sbg.bio.ic.ac.uk/~missense3d/index.html)
 




# 2019-06-24
001O_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
bcftools annotate -x ID, INFO, FORMAT LungBrain.vcf.gz -Oz -o LungBrain.trim.vcf.gz
bcftools annotate -x ID,FORMAT,INFO,FILTER 001O_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf | grep 29486505


 
bcftools annotate -Ov -I +'%ID' #leaves it as the existing ID
bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' #sets it to chr:pos:ref:alt



for i in `ls *snpEff.vcf`
do
bcftools norm -m-any $i | bcftools norm -Ov --check-ref w -f ~/hpc/db/hg19/hg19.fa | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS:%REF:%ALT' -Oz -o $i.trim.vcf.gz
done


for i in `ls *snpEff.vcf`
do
bcftools view -T ^remove.txt $i -Ov | bcftools norm -m-any -Ov --check-ref w -f ~/hpc/db/hg19/hg19.fa | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS' -Oz -o $i.trim.vcf.gz
done

for i in `ls *.gz`
do
bcftools query -l $i | grep VarD | grep Tumor
done

for i in `ls *snpEff.vcf`
do
bcftools view --force-samples -T ^remove.txt -S sample.txt $i -Ov | bcftools norm -m-any -Ov --check-ref w -f ~/hpc/db/hg19/hg19.fa | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS' -Oz -o $i.trim.vcf.gz
done

for i in `ls *.gz`
do
done

for i in `ls *.merged.HighConf.snpEff.vcf.trim.vcf.gz`
do
bcftools view --force-samples -S sample.txt $i -Ov -o $i.sin.vcf
done

for i in `ls *.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf`
do
perl updata.pl $i > $i.chan.vcf
done

for i in `ls *.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf`
do
bcftools sort $i -Oz -o $i.sort.vcf.gz
done

for i in `ls *mem.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf.sort.vcf.gz`
do
tabix -f -p vcf $i
done


bcftools query -l LungBrainAnnovar.hg19_multianno.trim.vcf | grep Tumor1 > Lung.txt
bcftools query -l LungBrainAnnovar.hg19_multianno.trim.vcf | grep Tumor2 > Brain.txt

bcftools view -S Lung.txt LungBrainAnnovar.hg19_multianno.vcf -Ov -o LungBrainAnnovar.Lung.hg19.vcf
bcftools view -S Brain.txt LungBrainAnnovar.hg19_multianno.vcf -Ov -o LungBrainAnnovar.Brain.hg19.vcf

bcftools annotate -x "^INFO/AN,INFO/AC,INFO/GENE,INFO/ExonicFunc.refGene" LungBrainAnnovar.Lung.hg19.vcf | less -S 

for i in `ls *.vcf`
do 
grep -v '#' $i | wc -l 
done

for i in `ls *mem.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf`
do 
echo $i
# grep -v '#' $i | wc -l 
done





ls *npEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf.sort.vcf.gz > input.new.txt
bcftools merge -l input.new.txt -Ov -o test.vcf
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') test.vcf -o test.anno.vcf

awk '{gsub(/^chr/,""); print}' test.vcf > test.temp.vcf
awk '{gsub(/=chr/,"="); print}' test.temp.vcf > test.update.vcf

bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') test.update.vcf -o test.anno.vcf


table_annovar.pl -vcfinput LungBrain.vcf ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 -out LungBrainAnnovar -remove -protocol refGene,dbnsfp33a -operation gx,f -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 

bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') test.update.vcf -o test.anno.vcf

bcftools annotate -x "^INFO/GENE,INFO/ExonicFunc.refGene" LungBrainAnnovar.hg19_multianno.vcf| less -S 


001E_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf

bcftools view -t ^'chr17:80396771' 001E_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
bcftools view -T ^remove.txt 001E_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf



bcftools norm -m-any 001C_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf -Ov -o test.vcf
| grep 120714401

##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=249250621>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER  INFO	FORMAT	Sample1
chr1	977330	.	T	C,G 225 PASS	.	GT	1/2


Error: wrong number of fields in FMT/F1R2 at chr2:120714401, expected 24, found 18
[E::vcf_parse_format] Number of columns at chr1:242383170 does not match the number of samples (5 vs 6)
Error: wrong number of fields in FMT/F1R2 at chr14:105181012, expected 18, found 12
Error: wrong number of fields in FMT/F1R2 at chr19:7977928, expected 18, found 12
[E::bcf_write] Broken VCF record, the number of columns at chr17:80396771 does not match the number of samples (0 vs 6)
Error: wrong number of fields in FMT/F1R2 at chr12:77216330, expected 18, found 12
[E::bcf_write] Broken VCF record, the number of columns at chr11:70317174 does not match the number of samples (0 vs 6)
Error: wrong number of fields in FMT/F1R2 at chr17:45664573, expected 18, found 12
[E::bcf_write] Broken VCF record, the number of columns at chr17:8093805 does not match the number of samples (0 vs 6)
Lines   total/split/realigned/skipped:  251/0/17/0
Error: wrong number of fields in FMT/F1R2 at chr17:45664573, expected 24, found 18
[E::vcf_parse_format] Number of columns at chr17:36357124 does not match the number of samples (4 vs 6)
Error: wrong number of fields in FMT/F1R2 at chr12:29486505, expected 18, found 12
[E::bcf_write] Broken VCF record, the number of columns at chr10:26534837 does not match the number of samples (0 vs 6)
Lines   total/split/realigned/skipped:  165/0/2/0






wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz
aloft --vcf=ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz --data ~/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/
cd aloft_output


for i in `ls *.lof`
do
awk '{print $1,$2,$3,$4,$5}' $i | sort -u | wc -l 
done

for i in `ls *.splice`
do
awk '{print $1,$2,$3,$4,$5}' $i | sort -u | wc -l 
done



for i in `ls *Stop*rmdup*.gz`
do
zcat $i | grep -v '#' | wc -l 
done

for i in `ls *Frame*rmdup*biallelic*.gz`
do
zcat $i | grep -v '#' | awk '{print $1,$2,$3}' | sort -u | wc -l 
done

for i in `ls *ExomeSplice*rmdup*biallelic*.gz`
do
zcat $i | grep -v '#' | awk '{print $1,$2,$3}' | sort -u | wc -l 
done

for i in `ls *Regulatory*rmdup*.gz`
do
zcat $i | grep -v '#' | wc -l 
done

for i in `ls *ExomeMissense.*.rmdup.biallelic.*gz`
do
zcat $i | grep -v '#' | awk '{print $1,$2,$3}' | sort -u | wc -l 
done


for i in `ls *ExomeMissense.*.rmdup.biallelic.*gz`
do
zcat $i | grep -v '#' | awk '{print $1,$2,$3}' | sort -u | wc -l 
done

for i in `ls *ExomeMissense.*.rmdup.biallelic.*gz`
do
echo $i
done

*.ExomeFrame.sort.rmdup.biallelic.vcf.gz
perl fixRefVCF.pl All_samples_Exome_QC.temp.vcf.recode.clean.chr21.vcf.recode.vcf  /gpfs/home/guosa/hpc/db/hg19/fa/chr21.fa  All_samples_Exome_QC.temp.vcf.recode.clean.chr21.vcf.fix.recode.vcf
perl fixRefVCF.pl All_samples_Exome_QC.temp.vcf.recode.clean.chr14.vcf.recode.vcf  ~/hpc/db/hg19/fa/chr14.fa  All_samples_Exome_QC.temp.vcf.recode.clean.chr14.vcf.fix.recode.vcf
14      rs1191543       30009227        A       G

cd ~/hpc/db/Gnomad/exome/aloft-exome-rec
mkdir aloft
mkdir temp
for i in 1 22
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mkdir aloft/chr$i >>$i.job
# echo bcftools view gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.vcf.gz >>$i.job
# echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.vcf.gz >>$i.job
# echo bcftools view -G gnomad.exomes.r2.1.sites.chr$i.rec.vcf.gz -Oz -o gnomad.exomes.r2.1.sites.chr$i.dq.rec.vcf.gz >> $i.job
echo ~/hpc/tools/aloft/aloft-annotate/aloft --vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.gz --output aloft/chr$i --data /gpfs/home/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/ >>$i.job
qsub $i.job
done


for i in {1..22} X Y
do
rm -rf chr$i
done

~/hpc/tools/aloft/aloft-annotate/aloft --vcf gnomad.exomes.r2.1.sites.chr1.rec.vcf.gz --output aloft/chr1 --data ~/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/

cd /gpfs/home/guosa/hpc/db/Gnomad/VCF
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo mkdir aloft/chr$i >>$i.job
# echo bcftools view gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.vcf.gz >>$i.job
# echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.vcf.gz >>$i.job
# echo bcftools view -G gnomad.exomes.r2.1.sites.chr$i.rec.vcf.gz -Oz -o gnomad.exomes.r2.1.sites.chr$i.dq.rec.vcf.gz >> $i.job
# echo zcat gnomad.exomes.r2.1.sites.chr$i.vcf.gz \| grep -v \'\#\' \| wc -l \> gnomad.exomes.r2.1.sites.chr$i.vcf.gz.number.txt >>$i.job
echo zcat gnomad.exomes.r2.1.1.sites.$i.vcf.bgz \| grep -v \'\#\' \| wc -l \> gnomad.exomes.r2.1.1.sites.$i.vcf.bgz.number.txt >>$i.job
qsub $i.job
done

for i in {557000..557023}.bright
do
qdel $i
done

54506 LOF (30685 LOF + 5993 Splice)

for i in `ls *.lof`
do
awk '{print $1,$2,$3,$4,$5}' $i  |sort -u | wc -l
done

for i in `ls *.splice`
do
awk '{print $1,$2,$3,$4,$5}' $i  |sort -u | wc -l
done


awk '{print $1,$2,$3,$4,$5}' gnomad.exomes.r2.1.sites.chr11.dq.rec.vcf.gz.vat.aloft.splice  | sort -u  > splice.txt
awk '{print $1,$2,$3,$4,$5}' gnomad.exomes.r2.1.sites.chr11.dq.rec.vcf.gz.vat.aloft.vcf | grep -v 'pos' | grep -v '#' > lof.vcf

data1<-read.table("lof.txt")
data2<-read.table("splice.txt")
data3<-read.table("lof.vcf")

data3=data3[-which(data3[,3] %in% data2[,3]),]
data3=data3[-which(data3[,3] %in% data1[,3]),]




      V1     V2           V3 V4    V5
1  chr11 180210 rs1279064406  A     G
2  chr11 180211 rs1323177122  A     G
8  chr11 193152  rs762487881  G   A,C
9  chr11 193154  rs201079312  C     T
12 chr11 193713  rs776547973  G A,C,T
20 chr11 193910 rs1189752583  T     C
21 chr11 193911  rs767802769  G     A
25 chr11 194418  rs371592375  G     A
26 chr11 194420 rs1437241589  G   A,T
29 chr11 197413  rs768048702  G   A,C
36 chr11 198464  rs766746003  G     A
37 chr11 198583  rs778946478  A     G
43 chr11 199507 rs1371341421  A     G
45 chr11 199946  rs772124586  G     A
48 chr11 205337  rs766186224  T     A
49 chr11 205344 rs1238039215  G   A,C
50 chr11 205345 rs1475051584  G     A
51 chr11 205432  rs775471050  C     T
53 chr11 205439 rs1435984912  G     A
54 chr11 205440 rs1300694307  C     T



mkdir annovar
mkdir temp
for i in X Y
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=8 >> chr$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> chr$i.job
echo \#PBS -m abe  >> chr$i.job
echo \#PBS -o $(pwd)/temp/ >>chr$i.job
echo \#PBS -e $(pwd)/temp/ >>chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq gnomad.exomes.r2.1.sites.chr$i.rec.vcf.gz  \> ./annovar/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.avinput >> chr$i.job
echo table_annovar.pl ./annovar/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ./annovar/chr$i -remove -protocol refGene,dbnsfp33a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,f,r,r,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done





# 2019-06-18

annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp30a   ~/hpc/tools/annovar/humandb/
annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp40a   ~/hpc/tools/annovar/humandb/
annotate_variation.pl --buildver hg19 --downdb seq /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_seq
retrieve_seq_from_fasta.pl /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_ensGene.txt -seqdir /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_seq -format ensGene -outfile /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_ensGeneMrna.fa
	
mkdir annovar
mkdir temp
for i in {1..22} 
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=8 >> chr$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> chr$i.job
echo \#PBS -m abe  >> chr$i.job
echo \#PBS -o $(pwd)/temp/ >>chr$i.job
echo \#PBS -e $(pwd)/temp/ >>chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq All_samples_Exome_QC.chr$i.vcf.vcf.gz  \> ./annovar/All_samples_Exome_QC.chr$i.avinput >> chr$i.job
echo table_annovar.pl ./annovar/All_samples_Exome_QC.chr$i.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ./annovar/chr$i -remove -protocol refGene,dbnsfp33a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,f,r,r,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done


mkdir annovar
mkdir temp
for i in {1..22} 
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=8 >> chr$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> chr$i.job
echo \#PBS -m abe  >> chr$i.job
echo \#PBS -o $(pwd)/temp/ >>chr$i.job
echo \#PBS -e $(pwd)/temp/ >>chr$i.job
echo cd $(pwd) >> chr$i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >>chr$i.job
qsub chr$i.job
done


echo convert2annovar.pl -format vcf4 -allsample -withfreq All_samples_Exome_QC.chr$i.vcf.vcf.gz  \> ./annovar/All_samples_Exome_QC.chr$i.avinput >> chr$i.job
echo table_annovar.pl ./annovar/All_samples_Exome_QC.chr$i.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ./annovar/chr$i -remove -protocol refGene,dbnsfp33a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,f,r,r,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job





plink2 --threads 32 --vcf All_samples_Exome_QC.temp.vcf.recode.clean.chr14.vcf.recode.vcf --a2-allele dbSNP152.hg19.fix.vcf 4 3 '#' --recode vcf --out All_samples_Exome_QC.temp.vcf.recode.clean.chr14.vcf.fix.recode.vcf 

wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.bgz -O ~/hpc/db/hg19/dbSNP152.hg19.vcf.bgz
tabix -p vcf dbSNP152.hg19.vcf.bgz
zcat dbSNP152.hg19.vcf.bgz > dbSNP152.hg19.vcf

wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz -O ~/hpc/db/hg38/dbSNP152.hg38.vcf.bgz
tabix -p vcf dbSNP152.hg38.vcf.bgz
zcat dbSNP152.hg19.vcf.bgz > dbSNP152.hg38.vcf.bgz

1) how many novel SNPs existed in GNMODA (EXOME and GENOME) compared with dbSNP152

bcftools 


wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz  -O ~/hpc/db/hg19/dbSNP151.hg19.vcf.bgz
awk '{print $1,$2,$3,$4,$5}' OFS="\t" ~/hpc/temp/LOF/prematureLOF.hg19.txt > prematureLOF.1K.hg19.vcf
bcftools view prematureLOF.1K.hg19.vcf -Oz -o prematureLOF.1K.hg19.vcf.gz
tabix -p vcf prematureLOF.1K.hg19.vcf.gz

bcftools view ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz.vat.aloft.vcf -Oz -o ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz.vat.aloft.vcf.gz
tabix -p vcf ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz.vat.aloft.vcf.gz


All_samples_Exome_QC

cd $HOME/.vep
curl -O ftp://ftp.ensembl.org/pub/release-96/variation/indexed_vep_cache/homo_sapiens_vep_96_GRCh37.tar.gz
tar xzf homo_sapiens_vep_96_GRCh38.tar.gz

cd ~/hpc/autism/data/aloft
awk '$4 !="-" && $5 !="-" {print}' OFS="\t" All_samples_Exome_QC.DG5.vcf > All_samples_Exome_QC.DG5.temp.vcf
vcftools --vcf All_samples_Exome_QC.DG5.temp.vcf --not-chr 0 --recode --out All_samples_Exome_QC.DG5.clean.vcf
aloft --vcf All_samples_Exome_QC.DG5.clean.vcf --output All_samples_Exome_QC.DG --data ~/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/


vcftools --vcf All_samples_Exome_QC.vcf --not-chr 0 --recode --out All_samples_Exome_QC.clean.vcf
bcftools norm -d both --threads=32 All_samples_Exome_QC.clean.vcf.recode.vcf -Ov  -o All_samples_Exome_QC.clean.norm.vcf
bcftools view -i 'ALT !="-" & REF !="-" ' All_samples_Exome_QC.clean.norm.vcf.gz -Oz -o All_samples_Exome_QC.norm.vcf.gz
tabix -p vcf All_samples_Exome_QC.norm.vcf.gz
sh remove_VCF_duplicates.sh All_samples_Exome_QC.norm.vcf.gz > All_samples_Exome_QC.clean.norm.undup.vcf
java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar impute=false gt=All_samples_Exome_QC.clean.norm.undup.vcf out=All_samples_Exome_QC.clean.norm.vcf.phasing


mkdir chr
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo vcftools --vcf All_samples_Exome_QC.clean.norm.undup.vcf --chr $i --recode --out ./chr/All_samples_Exome_QC.temp.vcf.recode.clean.chr$i.vcf >>$i.job
echo java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar gt=All_samples_Exome_QC.temp.vcf.recode.clean.chr$i.vcf.recode.vcf ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz map=~/hpc/db/hg19/beagle/plink.chr$i.GRCh37.map out=All_samples_Exome_QC.chr$i.vcf >>$i.job
qsub $i.job
done

for i in {1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -f -p vcf All_samples_Exome_QC.chr$i.vcf.vcf.gz >> $i.job
echo bcftools view -i \'AF \> 0\' All_samples_Exome_QC.chr$i.vcf.vcf.gz -Oz -o All_samples_Exome_QC.chr$i.vcf.gz >>$i.job
qsub $i.job
done

bcftools view -i '(IMP=1 & DR2>0.3)|IMP=0' All_samples_Exome_QC.chr22.vcf.gz -Oz -o test.vcf.gz
tabix -p vcf test.vcf.gz
plink --vcf [name of plink-exported VCF with incorrect reference alleles] --a2-allele [name of RefSeq file] [1-based column index of ref alleles] [1-based column index of variant IDs] --recode vcf --real-ref-alleles

mkdir chr
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo vcftools --vcf dbSNP152.hg19.fix.vcf --chr $i --recode --out ./dbSNP152/dbSNP152.hg19.chr$i.vcf >>$i.job
qsub $i.job
done

mkdir chr
mkdir temp
for i in {23..24} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo vcftools --vcf dbSNP152.hg19.fix.vcf --chr $i --recode --out ./dbSNP152/dbSNP152.hg19.chr$i.vcf >>$i.job
qsub $i.job
done



bcftools annotate -a /gpfs/home/guosa/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -x DS -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.vcf.gz -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.anno.vcf


grep -v '#' All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf | awk '{print $1"\t"$2}'

bgzip All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf
tabix -p vcf All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf.gz
bcftools view -R All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf.gz All_samples_Exome_QC.chr22.vcf.gz




# Stop_gain and Stop_loss 
cd ~/hpc/db/Gnomad/exome/aloft-exome-rec
mkdir vep
panel="ExomeStop"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -f PASS -i \'\(INFO\/vep \~ \"stop_gained\"\|\INFO\/vep \~ \"stop_lost\"\)\' gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz >>$i.job
echo bcftools sort ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -T ./temp/ >> $i.job
echo bcftools norm -d all ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz >> $i.job
echo bcftools view ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.gz >>$i.job
qsub $i.job
done

cd vep
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed


#frameshift_variant (not include  inframe_insertion, inframe_deletion )
cd /gpfs/home/guosa/hpc/db/Gnomad/exome/aloft-exome-rec
panel="ExomeFrame"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -f PASS -i \'\(INFO\/vep \~ \"frameshift\"\)\' gnomad.exomes.r2.1.sites.chr$i.rec.vcf.gz -Oz -o  ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz >>$i.job
echo bcftools sort ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -T ./temp/ >> $i.job
echo bcftools norm -d all ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz >> $i.job
echo bcftools view ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.gz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed


#splice_acceptor and splice_donor (not include   splice_region_variant)
cd /gpfs/home/guosa/hpc/db/Gnomad/exome/aloft-exome-rec
panel="ExomeSplice"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -f PASS -i \'\(INFO\/vep \~ \"splice_acceptor\"\|INFO\/vep \~ \" splice_donor\"\)\' gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz >>$i.job
echo bcftools sort ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -T ./temp/ >> $i.job
echo bcftools norm -d all ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz >> $i.job
echo bcftools view ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.gz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed

#missense_variant (not include) , 
cd ~/hpc/db/Gnomad/vcf
panel="ExomeMissense"
mkdir temp
for i in 3 5 7 8 12 15 16 19
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -f PASS -i \'\(INFO\/vep \~ \"missense_variant\"\)\' gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz >>$i.job
echo bcftools sort ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -T ./temp/ >> $i.job
echo bcftools norm -d all ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz >> $i.job
echo bcftools view ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.gz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed




#regulatory_region_variant  (not include) , 
cd /gpfs/home/guosa/hpc/db/Gnomad/exome/aloft-exome-rec
panel="ExomeRegulatory"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -f PASS -i \'\(INFO\/vep \~ \"regulatory_region_variant\"\)\' gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz >>$i.job
echo bcftools sort ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -T ./temp/ >> $i.job
echo bcftools norm -d all ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz >> $i.job
echo bcftools view ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.gz -Oz -o ./vep/gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed






awk '$4 !="-" && $5 !="-" {print}' OFS="\t" All_samples_Exome_QC.DG.vcf > All_samples_Exome_QC.DG5.temp.vcf
vcftools --vcf All_samples_Exome_QC.DG5.temp.vcf --not-chr 0 --recode --out All_samples_Exome_QC.DG5.clean.vcf
java -Xmx4g -jar ~/hpc/tools/snpEff/snpEff.jar GRCh37.75 All_samples_Exome_QC.DG5.clean.vcf.recode.vcf > All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.vcf

java -jar  ~/hpc/tools/snpEff/SnpSift.jar filter "ANN[0].EFFECT has 'missense_variant'" All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.vcf > All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.missense.vcf
java -jar  ~/hpc/tools/snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.vcf > All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.HIGH.vcf

bcftools sort All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.HIGH.vcf -o All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.HIGH.sort.vcf

perl -p -i -e 's/chr//' All_samples_Exome_QC.DG.vcf.vat.aloft.vcf
bcftools sort All_samples_Exome_QC.DG.vcf.vat.aloft.vcf -o All_samples_Exome_QC.DG.vcf.vat.aloft.sort.vcf

aloft --vcf All_samples_Exome_QC.DG.vcf --output All_samples_Exome_QC.DG --data /home/local/MFLDCLIN/guosa/hpc/tools/aloft/aloft-annotate/data.txt

cd /gpfs/home/guosa/hpc/autism/data/aloft
cd /gpfs/home/guosa/hpc/autism/data/2lof

grep -v '#' All_samples_Exome_QC.DG.vcf.vat.aloft.vcf | awk '{print $3}' > All_samples_Exome_QC.DG.vcf.vat.aloft.vcf.snpid
grep -v '#' All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.HIGH.sort.vcf | awk '{print $3}' > All_samples_Exome_QC.DG5.clean.vcf.recode.ANN.HIGH.sort.vcf.snpid
vcftools --vcf All_samples_Exome_QC.vcf --snps finalist.txt --recode --out All_samples_Exome_QC.LOF.vcf
bcftools annotate -a /gpfs/home/guosa/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') All_samples_Exome_QC.LOF.vcf.recode.vcf -o All_samples_Exome_QC.LOF.vcf.recode.anno.vcf
bcftools sort All_samples_Exome_QC.LOF.vcf.recode.anno.vcf -o All_samples_Exome_QC.LOF.vcf.recode.anno.sort.vcf


vcftools --vcf All_samples_Exome_QC.vcf --not-chr 0 --recode --out All_samples_Exome_QC.temp.vcf
awk '$4 !="-" && $5 !="-"' OFS="\t" All_samples_Exome_QC.temp.vcf.recode.vcf > All_samples_Exome_QC.temp.vcf.recode.clean.vcf

mkdir chr
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo vcftools --vcf All_samples_Exome_QC.temp.vcf.recode.clean.vcf --chr $i --recode --out ./chr/All_samples_Exome_QC.temp.vcf.recode.clean.chr$i.vcf >>$i.job
qsub $i.job
done


cd ~/hpc/db/hg19/beagle
for i in {1..22} X Y
do
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$i.1kg.phase3.v5a.vcf.gz
done
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/20140625_related_individuals.txt
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_male_samples_v3.20130502.ALL.panel
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel

mkdir EUR
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf chr$i.1kg.phase3.v5a.vcf.gz >> $i.job
echo bcftools view chr$i.1kg.phase3.v5a.vcf.gz -S EUR.List.txt -Oz -o ./EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz >>$i.job
qsub $i.job
done

mkdir temp
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar gt=All_samples_Exome_QC.temp.vcf.recode.clean.chr$i.vcf.recode.vcf ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz map=~/hpc/db/hg19/beagle/plink.chr$i.GRCh37.map out=All_samples_Exome_QC.chr$i.vcf >>$i.job
qsub $i.job
done


for i in X Y
do
java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar gt=All_samples_Exome_QC.temp.vcf.recode.clean.chr$i.vcf.recode.vcf ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz map=~/hpc/db/hg19/beagle/plink.chr$i.GRCh37.map out=All_samples_Exome_QC.chr$i.vcf  &
done


rvtest --inVcf input.vcf --pheno phenotype.ped --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac




for i in {1..45} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo rvtest --inVcf All_samples_Exome_QC.gz --pheno All_samples_Exome_QC.phen --mpheno $i --single wald,score --out All_samples_Exome_QC.mphen$i.rvtest >>$i.job
qsub $i.job
done


wget http://qbrc.swmed.edu/zhanxw/seqminer/data/refFlat_hg19.txt.gz
gunzip refFlat_hg19.txt.gz
perl -p -i -e 's/chr//' refFlat_hg19.txt
gzip refFlat_hg19.txt

for i in {1..45} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf All_samples_Exome_QC.chr$i.vcf.vcf.gz >>$i.job
echo rvtest --inVcf All_samples_Exome_QC.gz --pheno All_samples_Exome_QC.phen --mpheno $i --covar All_samples_Exome_QC.cov --covar-name AvgAge --out All_samples_Exome_QC.rvtest.burden.mphen$i --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac >>$i.job
qsub $i.job
done

rvtest --inVcf All_samples_Exome_QC.gz --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov --covar-name AvgAge --gene NEGR1 --out NEGR1 --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac


cat All_samples_Exome_QC.chr*.rvtest.skat.CMC.assoc >> All_samples_Exome_QC.rvtest.skat.CMC.assoc
cat All_samples_Exome_QC.chr*.rvtest.skat.VariableThresholdPrice.assoc >> All_samples_Exome_QC.rvtest.skat.VariableThresholdPrice.assoc
cat All_samples_Exome_QC.chr*.rvtest.skat.Skat.assoc >> All_samples_Exome_QC.rvtest.skat.Skat.assoc
cat All_samples_Exome_QC.chr*.rvtest.skat.Kbac.assoc >> All_samples_Exome_QC.rvtest.skat.Kbac.assoc

awk '$7<0.0005' All_samples_Exome_QC.rvtest.skat.CMC.assoc  | wc -l 
awk '$13<0.0005' All_samples_Exome_QC.rvtest.skat.VariableThresholdPrice.assoc | wc -l 
awk '$13<0.0005' All_samples_Exome_QC.rvtest.skat.Skat.assoc | wc -l 
awk '$6<0.0005' All_samples_Exome_QC.rvtest.skat.Kbac.assoc | wc -l 

awk '$7<0.0005' All_samples_Exome_QC.rvtest.skat.CMC.assoc  
awk '$13<0.0005' All_samples_Exome_QC.rvtest.skat.VariableThresholdPrice.assoc 
awk '$13<0.00000000005' *.Skat.assoc
awk '$7<0.00000000005' *.Skat.assoc

awk '$6<0.0005' All_samples_Exome_QC.rvtest.skat.Kbac.assoc 






for i in {1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf All_samples_Exome_QC.chr$i.vcf.vcf.gz >>$i.job
qsub $i.job
done

bcftools concat -f concate.list.txt -Oz -o All_samples_Exome_QC.gz











mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'DR2\>0.6\' All_samples_Exome_QC.chr$i.vcf.vcf.gz -Ov -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf >>$i.job
qsub $i.job
done


mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view  All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf -Oz -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf.gz >>$i.job
echo tabix -p vcf All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf.gz >>$i.job
echo bcftools view -G All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf.gz -Ov -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.vcf >>$i.job
echo ~/hpc/tools/aloft/aloft-annotate/aloft --vcf All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.vcf --output All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.vcf.aloft --data /gpfs/home/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/ >>$i.job
qsub $i.job
done

mkdir temp
for i in 7
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo ~/hpc/tools/aloft/aloft-annotate/aloft --vcf All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.vcf --output All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.vcf.aloft --data /gpfs/home/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/ >>$i.job
qsub $i.job
done

mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx4g -jar ~/hpc/tools/snpEff/snpEff.jar GRCh37.75 All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.vcf \> All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.snpEff.vcf >>$i.job
qsub $i.job
done

 
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -jar  ~/hpc/tools/snpEff/SnpSift.jar filter \"\(ANN[0].IMPACT has \'HIGH\'\)\&\& \(AF\< 0.2\)\" All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.snpEff.vcf \> All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.snpEff.High.vcf >>$i.job
qsub $i.job
done


cat *.vcf.DR2L0.8.DG.snpEff.High.vcf | grep -v '#' |awk '{print $1,$2}' OFS="\t" > LOF.hg19.txt
cat All_samples_Exome_QC.chr*.vcf.DR2L0.8.DG.vcf.vat.aloft.vcf | grep -v '#' |awk '{print $1,$2}' OFS="\t" | sort -u >>  LOF.hg19.txt
perl -p -i -e 's/chr//' LOF.hg19.txt
sort -u LOF.hg19.txt > LOF.hg19.sort.txt

cat All_samples_Exome_QC.chr*.vcf.DR2L0.8.DG.snpEff.High.vcf | grep -v '#' |awk '{print $3}' OFS="\t" > LOF.hg19.txt
cat All_samples_Exome_QC.chr*.vcf.DR2L0.8.DG.vcf.vat.aloft.vcf | grep -v '#' |awk '{print $3}' OFS="\t" | sort -u >>  LOF.hg19.txt
perl -p -i -e 's/chr//' LOF.hg19.txt
sort -u LOF.hg19.txt | grep -v '[.]' > LOF.hg19.sort.txt


for i in {1..22}
do
vcftools --vcf All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf --snps LOF.hg19.sort.txt --recode --stdout | bgzip -c > All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.vcf.gz
tabix -p vcf All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.vcf.gz
bcftools annotate -a /gpfs/home/guosa/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -x DS -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.vcf.gz -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.anno.vcf
done


for i in {1..22}
do
bcftools annotate -x FORMAT/DS -a /gpfs/home/guosa/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.vcf.gz -o All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.anno.vcf
done

for i in {1..22}
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen >> $i.job
qsub $i.job
done


for i in {1..22}
do
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr$i.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
done

Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr1.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr2.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr3.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr4.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr5.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr6.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr7.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr8.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr9.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr10.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr11.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr12.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr13.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr14.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr15.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr16.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr17.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr18.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr19.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr20.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr21.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen
Rscript --vanilla 2LOF.R All_samples_Exome_QC.chr22.vcf.DR2L0.8.lof.anno.vcf All_samples_Exome_QC.phen

cp ~/hpc/project/pmrp/Exom2/2LOF/chr22.update.vcf  ./
cp ~/hpc/project/pmrp/Exom2/2LOF/2LOF.R  ./


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF

for i in {1..22}
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp10_Obesity_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp1_RA_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp7_Iron_C1_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp7_Iron_C2_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp3_PA_rev2_SampleIDs.Michigen.txt >> $i.job 
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp5_Thyroid_C1_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp5_Thyroid_C2_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_SSc_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ANA_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ENA_rev2_SampleIDs.Michigen.txt >> $i.job
qsub $i.job
done







 
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -jar  ~/hpc/tools/snpEff/SnpSift.jar filter \"\(ANN[0].IMPACT has \'HIGH\'\)\&\& \(AF\< 0.2\)\" All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.snpEff.vcf \> All_samples_Exome_QC.chr$i.vcf.DR2L0.8.DG.snpEff.High.vcf >>$i.job
qsub $i.job
done


zcat gnomad.genomes.r2.1.sites.chr20.vcf.bgz | grep 'frame' | awk '{print $3}' > 
#frameshift_variant (not include  inframe_insertion, inframe_deletion )
cd ~/hpc/db/Gnomad/vcf
panel="ExomeSplice"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"frameshift\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -T ./temp/ >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed


zcat gnomad.genomes.r2.1.sites.chr20.vcf.bgz | grep 'splice' | awk '{print $3}' > splice.rs.txt
#frameshift_variant (not include  inframe_insertion, inframe_deletion )
cd ~/hpc/db/Gnomad/vcf
panel="ExomeFrame"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"frameshift\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -T ./temp/ >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed



file=list.files(pattern="*.predict")
output<-c()
for(i in 1:22){
input=paste("chr",i,".vcf.vat.aloft.lof.predict",sep="")
data<-read.table(input,head=F,sep="\t")
newdata<-subset(data,data[,13]<0.05 & data[,15] !="Tolerant")[,c(1,2,9,15)]
output<-rbind(output,newdata)
}
colnames(output)<-c("CHR","POS","GENE","MODEL")
write.table(output,file="aloft.hg19.txt",col.names=T,row.names=F,quote=F,sep="\t")
perl -p -i -e 's/chr//i' aloft.hg19.txt
perl -p -i -e 's/Dominant/Dom/i' aloft.hg19.txt
perl -p -i -e 's/Recessive/Rec/i' aloft.hg19.txt
gzip aloft.hg19.txt




# Tutorial: Comparison between VEP and ALoFT
i="chr21"
panel="ExomeStop"
bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"stop\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ov -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf
aloft --vcf gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf --data /home/local/MFLDCLIN/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/  --output gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop
perl -p -i -e 's/chr//i' gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf.vat
perl -p -i -e 's/chr//i' gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf.vat.aloft.vcf
perl -p -i -e 's/chr//i' gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf.vat.aloft.lof
perl -p -i -e 's/chr//i' gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf.vat.aloft.splice
vcftools --vcf gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf --diff gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf.vat --diff-site-discordance
awk '$3==1' out.diff.sites | wc -l
awk '$3==1{print $1"\t"$2}' out.diff.sites > Target.txt
bgzip -c Target.txt > Target.txt.gz && tabix -s1 -b2 -e2 Target.txt.gz
bcftools view -T Target.txt.gz gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf | grep -v '#' | awk '{print $3}'






cd /home/guosa/hpc/autism/data/2LOF
avinput="All_samples_Exome_QC.vcf"
table_annovar.pl $avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out $avinput -remove -protocol refGene,dbnsfp35c,gnomad_genome,gwasCatalog,gtexEqtlCluster -operation g,f,f,r,r -nastring . 

bcftools_view view gnomad.genomes.r2.1.sites.chr1.vcf.bgz -t rs28362286 
bcftools view gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz -t rs764400971 

https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.1.vcf.bgz

for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done

cd /home/guosa/hpc/db/Gnomad/vcf
panel="refGene"
mkdir temp
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed



cd /home/guosa/hpc/db/Gnomad/vcf
panel="LOF"
mkdir temp
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice_acce\" \| INFO\/vep \~ \"splice_d\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo tabix -p vcf  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >> $i.job
echo bcftools sort gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.sort.rec.$panel.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.genomes.r2.1.sites.chr$i.sort.rec.$panel.vcf.bgz  -Oz -o gnomad.genomes.r2.1.sites.chr$i.sort.rmdup.rec.$panel.vcf.bgz  >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.genomes.r2.1.sites.chr$i.sort.rmdup.rec.$panel.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.sort.rmdup.biallelic.rec.$panel.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.final.vcf > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.genomes.r2.1.sites.rec.$panel.hg19.vcf.bed


panel="LOF"
mkdir temp
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo bcftools norm -m \+ ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
# echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice_acce\" \| INFO\/vep \~ \"splice_d\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  gnomad.exomes.r2.1.sites.chr$i.$panel.rec.vcf.bgz >>$i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.$panel.rec.vcf.bgz >> $i.job 
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.$panel.rec.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.$panel.sort.rec.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.$panel.sort.rec.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.$panel.sort.rmdup.rec.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.$panel.sort.rmdup.rec.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.$panel.sort.rmdup.biallelic.rec.vcf >>$i.job
qsub $i.job
done

ls *$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.$panel.hg19.vcf.bed


panel="LOF" 
mkdir temp
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed

#PBS -N 22
#PBS -l nodes=1:ppn=1
#PBS -M Guo.shicheng@marshfieldresearch.org
#PBS -m abe
#PBS -o /gpfs/home/guosa/hpc/db/Gnomad/vcf/temp/
#PBS -e /gpfs/home/guosa/hpc/db/Gnomad/vcf/temp/
cd /gpfs/home/guosa/hpc/db/Gnomad/vcf
bcftools view -v snps -f PASS -i '(INFO/vep ~ "stop" | INFO/vep ~ "lost" | INFO/vep ~ "splice_acce" | INFO/vep ~ "splice_d" | INFO/vep ~ "gain" | INFO/vep ~ "frame")' /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr22.LOF.rec.vcf.bgz
tabix -p vcf gnomad.exomes.r2.1.sites.chr22.LOF.rec.vcf.bgz
bcftools sort gnomad.exomes.r2.1.sites.chr22.LOF.rec.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr22.LOF.sort.rec.vcf
bcftools norm -d all gnomad.exomes.r2.1.sites.chr22.LOF.sort.rec.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr22.LOF.sort.rmdup.rec.vcf.bgz
bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr22.LOF.sort.rmdup.rec.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr22.LOF.sort.rmdup.biallelic.rec.vcf


# 2019-06-04

for i in `ls *H3k27ac*.bigWig`
do
bigWigToBedGraph $i $i.hg19.bed
done

for i in `ls *H3k4me1*.bigWig`
do
bigWigToBedGraph $i $i.hg19.bed
done

for i in `ls *H3k4me3*.bigWig`
do
bigWigToBedGraph $i $i.hg19.bed
done

for i in `ls *H3k27ac*.bigWig.hg19.bed`
do
cat $i >> wgEncodeBroadHistone.H3k27ac.hg19.bed
done

for i in `ls *H3k4me1*.bigWig.hg19.bed`
do
cat $i >> wgEncodeBroadHistone.H3k4me1.hg19.bed
done

for i in `ls *H3k4me3*.bigWig.hg19.bed`
do
cat $i >> wgEncodeBroadHistone.H3k4me3.hg19.bed
done

for i in `ls *wgEncodeBroadHistone.*.hg19.bed`
do
echo $i
bedtools sort -i $i > $i.sort.bed &
done

bedtools merge -i wgEncodeBroadHistone.H3k27ac.hg19.bed > wgEncodeBroadHistone.H3k27ac.merge.hg19.bed &
bedtools merge -i wgEncodeBroadHistone.H3k4me1.hg19.bed > wgEncodeBroadHistone.H3k4me1.merge.hg19.bed &
bedtools merge -i wgEncodeBroadHistone.H3k4me3.hg19.bed > wgEncodeBroadHistone.H3k4me3.merge.hg19.bed &

cp wgEncodeBroadHistone.H3k27ac.merge.hg19.bed ../
cp wgEncodeBroadHistone.H3k4me1.merge.hg19.bed ../
cp wgEncodeBroadHistone.H3k4me3.merge.hg19.bed ../

GPL13534_450K_hg19_PBMC_BUR.bed
wgEncodeRegDnaseClusteredV3.hg19.bed
wgEncodeBroadHistone.H3k27ac.merge.hg19.bed
wgEncodeBroadHistone.H3k4me1.merge.hg19.bed
wgEncodeBroadHistone.H3k4me3.merge.hg19.bed




wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_expected_count.txt.gz
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_exon_reads.parquet

Anaplastic astrocytoma (WHO grade III)
Anaplastic Astrocytoma (WHO grade III)
Diffuse astrocytoma (WHO grade II)
Diffuse Astrocytoma (WHO grade II)
Ganglioneuroblastom
Hirnstammgliom
Maligner glialer Tumor
Pilocytic Astrocytoma (WHO grade I)



python configureStrelkaSomaticWorkflow.py --normalBam baihaixing_control.recal.bam --tumorBam baihaixing_T1.recal.bam --referenceFasta /thinker/storage/udata/bing/GATK/hg19/ucsc.hg19.fasta --exome --callRegions /thinker/storage/udata/bing/bed/xGenExomeResearchPaneltargets.sorted.ip100.merged.bed.gz --runDir baihaixing_T1_strelka_somatic


mkdir temp
perl -p -i -e "{s/chr//}" $panel.hg19.bed
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ ~/hpc/project/pmrp/phase1/imputation/chr$i.dose.filter.vcf.gz -Oz -o chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -f PASS -i \'INFO\/STATUS \' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
bcftools view -i 'INFO/STATUS ~ "Somatic" | INFO/STATUS ~ "LOH"'  
bcftools view -f ., PASS P2_T3.strelka.somatic.vcf.gz
for i in `ls *.strelka.somatic.vcf.gz`; do bcftools view -f PASS $i | wc -l ; done
for i in `ls *.strelka.somatic.vcf.gz`
do
bcftools view -H -f PASS $i | awk '{print $1,$2,$3,$4,$5}' OFS="\t" >> prostate.vep.input
echo $i
done


for i in `ls *.strelka.somatic.vcf.gz.pass.vcf.gz.vardit.overlap.vcf.gz`
do
bcftools view -H $i | awk '{print $1,$2,$3,$4,$5}' OFS="\t" >> prostate.vep.input
echo $i
done






### pick up all the somatic mutations from strelka
for i in `ls *.strelka.somatic.vcf.gz`
do
bcftools view -f PASS $i -o $i.pass.vcf 
bgzip $i.pass.vcf
tabix -p vcf $i.pass.vcf.gz
bcftools view -H -f PASS $i.pass.vcf.gz | awk '{print $1,$2-1,$2,$1":"$2-1"-",$3,$4,$5}' OFS="\t" | sort -u > $i.bed
echo $i
done 

### pick up all the somatic mutations from vardict
for i in `ls *.vardict.vcf.gz`
do
bcftools view -i 'INFO/STATUS !="Germline"' $i -o $i.pass.vcf 
bgzip $i.pass.vcf
tabix -p vcf $i.pass.vcf.gz
bcftools view -H -f PASS $i.pass.vcf.gz |awk '{print $1,$2-1,$2,$1":"$2-1"-"$2,$3,$4,$5}' OFS="\t" |sort -u > $i.bed
echo $i
done 

###### merge vardict and strelka
for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
bedtools intersect -wa -a $i.vardict.vcf.gz.bed -b $i.strelka.somatic.vcf.gz.bed | sort -u > $i.vardict.strelka.overlap.bed
done

wc -l *vardict.strelka.overlap.bed > mutation.counts.txt

for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
echo $i
done

###### merge vardict and strelka
for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
bedtools intersect -wa -a $i.vardict.vcf.gz.bed -b $i.strelka.somatic.vcf.gz.bed | sort -u | awk '{print $1,$3}' OFS="\t" > $i.vardict.strelka.overlap.pos
echo $i
done


ls *.vcf.gz.bed | awk -F"." '{print $1}' | sort -u > FileUnit.txt
cat *vardict.strelka.overlap.bed | awk '{print $1,$3}' OFS="\t" | sort -u > valid.pos.txt

for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
bcftools view -T $i.vardict.strelka.overlap.pos $i.strelka.somatic.vcf.gz.pass.vcf.gz -Oz -o $i.strelka.vardit.overlap.vcf.gz
tabix -p vcf $i.strelka.vardit.overlap.vcf.gz
echo $i
done

for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
bcftools view -T $i.vardict.strelka.overlap.pos $i.vardict.vcf.gz.pass.vcf.gz -Oz -o $i.vardict.strelka.overlap.vcf.gz
tabix -p vcf $i.vardict.strelka.overlap.vcf.gz
echo $i
done

mkdir temp
for i in `ls *strelka.vardit.overlap.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx16g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> $i.anno.vcf >> $i.job 
qsub $i.job
done

mkdir temp
for i in `ls *.vardict.strelka.overlap.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx16g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> $i.anno.vcf >> $i.job 
qsub $i.job
done

#################################################################################################################################################
#################################################################################################################################################
for i in `ls *.strelka.vardit.overlap.vcf.gz.anno.vcf`
do
java -jar ~/hpc/tools/snpEff/SnpSift.jar extractFields $i CHROM POS REF ALT "ANN[*].GENE:" | awk -F'[\t|]' '{print $1,$2,$3,$4,$8}' OFS="\t" > $i.hg19.anno.bed
java -jar ~/hpc/tools/snpEff/SnpSift.jar tstv $i > $i.hg19.tstv.bed
done

for i in `ls *vardict.strelka.overlap.vcf.gz.anno.vcf`
do
grep -v '#' $i | perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1\t$2\t$3\t$4\t$5" if $_=~/TYPE=(\w+);.*\|(\w+)\|(\w+)\|(\w+)\|(ENSG\w+)/}' > $i.hg19.bed
echo $i
done

for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
grep -v '#' $i.vardict.strelka.overlap.vcf.gz.anno.vcf| perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1" if $_=~/TYPE=(\w+);/}' > $i.vardit.strelka.overlap.type.hg19.bed
echo $i
done


mkdir temp
for i in `ls *vardict.strelka.overlap.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx16g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> $i.anno.vcf >> $i.job 
qsub $i.job
done

for i in `ls *vardict.strelka.overlap.vcf.gz.anno.vcf` 
do
grep -v '#' $i | perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1\t$2\t$3\t$4\t$5" if $_=~/TYPE=(\w+);.*\|(\w+)\|(\w+)\|(\w+)\|(ENSG\w+)/}' > $i.hg19.bed
echo $i
done

grep -v '#' P3_T2.vardict.strelka.overlap.vcf.gz.anno.vcf | perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1\t$2" if $_=~/TYPE=(\w+);.*(ENSG\w+)/}'






###### mapping genomic to gene symbol-14185


for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14ming_T1 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
bedtools intersect -wa -a $i.vardict.vcf.gz.bed -b $i.strelka.somatic.vcf.gz.bed > $i.vardict.strelka.overlap.bed
wc -l $i.vardict.strelka.overlap.bed
done






for i in P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P14ming P15 P16
do
for j in T1 T2 T3 T4
do
for k in T1 T2 T3 T4
do
bedtools intersect -wa -a $i\_$j.vardict.strelka.overlap.bed -b $i\_$k.vardict.strelka.overlap.bed | wc -l 
done
done 
done

bcftools merge --force-samples P1_T1*pass*.gz -Oz -o P1_T1.vcf.gz
tabix -p vcf P1_T1.vcf.gz
bcftools view -i 'TUMOR !="./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:."'  P1_T1.vcf.gz

perl vep_merge2.pl > prostate.somaticMutation.gene.txt
perl -p -i -e 's/;-//g' prostate.somaticMutation.gene.txt

mkdir temp
for i in `ls *strelka.somatic.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mkdir $i.rlt >> $i.job
echo java -Xmx16g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> ./$i.rlt/$i.anno.vcf >> $i.job 
qsub $i.job
done

java -jar ~/hpc/tools/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar FlagStat
java -jar ~/hpc/tools/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar FlagStat
	 

for i in `ls *strelka.somatic.vcf.gz`
do
bcftools stats -f PASS $i | grep TSTV | grep  '#' >> TSTV.header.txt
done


for i in `ls *strelka.somatic.vcf.gz`
do
bcftools stats -f PASS $i | grep TSTV | grep -v '#' >> TSTV.txt
done

for i in `ls *strelka.somatic.vcf.gz`
do
echo $i >> TSTV.Filename.txt
done

###############################################################################################################################
!/usr/bin

/mnt/bigdata/Genetic/Projects/Jixia_Liu/BreastCancer/Harold/condor_compute

#!/bin/bash

## change attribute rvtest to executable
chmod a+x rvtest

## define computation parameters
PHEN_NAME="disease"
COV_NAME="sex,erh,genopl"  ## "sex,age,erh"
STAT_TEST="SKATO"

VCFIN="BC_gene_all_variant.vcf.gz"
PEDIN="phe$1.in"
SETIN="BC_gene_Harold_set.list"

ASSCOUT="BC_gene_all_variant_EX_ALL_SKATO_DX$1"

## with dosage
./rvtest --pheno $PEDIN --pheno-name $PHEN_NAME --covar $PEDIN --covar-name $COV_NAME --inVcf $VCFIN --dosage DS --setFile $SETIN --kernel $STAT_TEST --out $ASSCOUT #$SKATOUT

## add error indicator

rvtestexit=$?

if [ $rvtestexit -ne 0 ]; then
  exit $rvtestexit;
fi

#gzip $ASSCOUT'.SingleFirth.assoc'
#rm $ASSOOUT'.SingleScore.assoc'

gzip $ASSCOUT'.SkatO.assoc'

rm $VCFIN
rm $VCFIN'.tbi'


#################################################################################################################################
perl -lane '{print $1 if $_=~/(GPL\d+)/}' gds_result.txt | sort -u
GPL10332
#################################################################################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/IBD

grep FSTL1 ~/hpc/db/hg19/refGene.hg19.bed > RAcandidate.hg19.bed
grep CCR7 ~/hpc/db/hg19/refGene.hg19.bed >> RAcandidate.hg19.bed
# perl -lane "{print if /\s+IL\d+/}"  ~/hpc/db/hg19/refGene.hg19.bed >>  RAcandidate.hg19.bed
# perl -lane "{print if /\s+MUC\d+/}"  ~/hpc/db/hg19/refGene.hg19.bed >>  RAcandidate.hg19.bed
perl -p -i -e "s/chr//" RAcandidate.hg19.bed
sort -u RAcandidate.hg19.bed > RAcandidate.uni.hg19.bed
mv RAcandidate.uni.hg19.bed  RAcandidate.hg19.bed
panel="refGene"
mkdir temp
perl -p -i -e "{s/chr//}" $panel.hg19.bed
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ ~/hpc/project/pmrp/phase1/imputation/chr$i.dose.filter.vcf.gz -Oz -o chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF\>0.001 \& INFO/AF\<0.999 \& \(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' -R $panel.hg19.bed ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
#################################################################################################################################
cd ~/hpc/rheumatology/RA/IBD
i=3
bcftools norm -m + ~/hpc/project/pmrp/phase1/imputation/chr$i.dose.filter.vcf.gz -Oz -o chr$i.rec.vcf.bgz
tabix -p vcf chr$i.rec.vcf.bgz
tabix -p vcf ~/hpc/project/pmrp/phase1/imputation/chr$i.dose.filter.vcf.gz 3:120080201-120285582 -h > chr$i.rec.fstl1.vcf.bgz
bcftools view -i 'R2>0.8' chr$i.rec.fstl1.vcf.bgz -Ov -o chr$i.rec.fstl1.slim.vcf

# rvtest
rvtest --noweb --inVcf ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --pheno ~/hpc/project/pmrp/phase1/pmrp.exom1.2ALOF.phen --pheno-name PheTyp1_RA_C1 --out PheTyp1_RA_C1 --geneFile ~/hpc/project/pmrp/phase1/refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac
rvtest --noweb --inVcf ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --pheno ~/hpc/project/pmrp/phase1/pmrp.exom1.2ALOF.phen --gene FSTL1 --pheno-name PheTyp1_RA_C2 --out PheTyp1_RA_C2 --geneFile ~/hpc/project/pmrp/phase1/refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac

# plink + Chi-Square
cd /home/guosa/hpc/rheumatology/RA/IBD/plink
plink --bfile ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --pheno ~/hpc/project/pmrp/phase1/pmrp.exom1.2ALOF.phen --pheno-name PheTyp1_RA_C1 --assoc counts --out PheTyp1_RA_C1
plink --bfile ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --pheno ~/hpc/project/pmrp/phase1/pmrp.exom1.2ALOF.phen --pheno-name PheTyp1_RA_C2 --assoc counts --out PheTyp1_RA_C2
awk '{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" PheTyp1_RA_C1.assoc | grep -v "NA" | grep -v CHISQ > PheTyp1_RA_C1.hg19.bed
awk '{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" PheTyp1_RA_C2.assoc | grep -v "NA" | grep -v CHISQ > PheTyp1_RA_C2.hg19.bed
# FSTL1_hg19.bed chr3	120080201	120285582	FSTL1
bedtools intersect -header -a PheTyp1_RA_C1.hg19.bed -b FSTL1.hg19.bed > PheTyp1_RA_C1.FSTL1.hg19.bed
bedtools intersect -header -a PheTyp1_RA_C1.hg19.bed -b FSTL1.hg19.bed > PheTyp1_RA_C2.FSTL1.hg19.bed
head -n 1 PheTyp1_RA_C1.assoc | awk '{{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}}' OFS="\t" > PheTyp1_RA_C1.FSTL1.assoc.hg19.bed
head -n 1 PheTyp1_RA_C2.assoc | awk '{{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}}' OFS="\t" > PheTyp1_RA_C2.FSTL1.assoc.hg19.bed
cat PheTyp1_RA_C1.FSTL1.hg19.bed >> PheTyp1_RA_C1.FSTL1.assoc.hg19.bed
cat PheTyp1_RA_C2.FSTL1.hg19.bed >> PheTyp1_RA_C2.FSTL1.assoc.hg19.bed

# plink + Fisher
cd /home/guosa/hpc/rheumatology/RA/IBD/plink
plink --bfile ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --pheno ~/hpc/project/pmrp/phase1/pmrp.exom1.2ALOF.phen --pheno-name PheTyp1_RA_C1 --assoc fisher  --out PheTyp1_RA_C1
plink --bfile ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --pheno ~/hpc/project/pmrp/phase1/pmrp.exom1.2ALOF.phen --pheno-name PheTyp1_RA_C2 --assoc fisher  --out PheTyp1_RA_C2
awk '{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" PheTyp1_RA_C1.assoc.fisher | grep -v "NA" | grep -v CHISQ | grep -v "SNP" > PheTyp1_RA_C1.hg19.bed
awk '{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" PheTyp1_RA_C2.assoc.fisher | grep -v "NA" | grep -v CHISQ | grep -v "SNP" > PheTyp1_RA_C2.hg19.bed
# FSTL1_hg19.bed chr3	120080201	120285582	FSTL1
bedtools intersect -a PheTyp1_RA_C1.hg19.bed -b FSTL1.hg19.bed > PheTyp1_RA_C1.FSTL1.hg19.bed
bedtools intersect -a PheTyp1_RA_C2.hg19.bed -b FSTL1.hg19.bed > PheTyp1_RA_C2.FSTL1.hg19.bed
head -n 1 PheTyp1_RA_C1.assoc.fisher | awk '{{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}}' OFS="\t" > PheTyp1_RA_C1.FSTL1.assoc.hg19.bed
head -n 1 PheTyp1_RA_C2.assoc.fisher | awk '{{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11}}' OFS="\t" > PheTyp1_RA_C2.FSTL1.assoc.hg19.bed
cat PheTyp1_RA_C1.FSTL1.hg19.bed >> PheTyp1_RA_C1.FSTL1.assoc.hg19.bed
cat PheTyp1_RA_C2.FSTL1.hg19.bed >> PheTyp1_RA_C2.FSTL1.assoc.hg19.bed
head PheTyp1_RA_C1.FSTL1.assoc.hg19.bed
head PheTyp1_RA_C2.FSTL1.assoc.hg19.bed

awk '{print $4}'  PheTyp1_RA_C1.FSTL1.assoc.hg19.bed | grep -v SNP > FSTL1.SNPs.txt
plink --bfile ~/hpc/project/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --extract FSTL1.SNPs.txt --make-bed --out FSTL1
plink --bfile FSTL1 --allow-no-sex --hap FSTL1.hlist 


#################################################################################################################################

awk '{print "chr"$2,$3-1,$3,$2":"$3}' OFS="\t" ASA1.txt | grep -v "NAME" | grep -v 'POS' > ASA1.bed
awk '{print "chr"$2,$3-1,$3,$2":"$3}' OFS="\t" ASA2.txt | grep -v "INFO" | grep -v 'POS' > ASA2.bed
grep -v '?' ASA1.bed | sort -u > ASA1.hg19.bed 
grep -v '?' ASA2.bed | sort -u > ASA2.hg19.bed 
wc -l ASA1.hg19.bed 
wc -l ASA2.hg19.bed 
bedtools intersect -a ASA1.hg19.bed -b ASA2.hg19.bed > ASA.bed
wc -l ASA.bed

track type=bedGraph name=H1299_NSCLC description=Human_NSCLC_H1299 
 autoScale=off maxHeightPixels=128:64:32 visibility=full 

browser position chr19:49302001-49304701
browser pack refGene encodeRegions
track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=0,100,200 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full 
chr19 49302000 49302300 -1.0
chr19 49302300 49302600 -0.75

# bw2wig
mkdir temp
for i in `ls *bw`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job 
echo bigWigToWig $i $i.wig >> ./temp/$i.job
qsub ./temp/$i.job
done

awk '{print FILENAME,$7,$1":"$11,$13}' OFS="\t" *multianno.txt.bed.vep.bed > cBioPortal.input.txt

ADRA1A chr8:26,723,345-26,723,990
# PBMC 
wget http://smithdata.usc.edu/methbase/data/Li-PBMC-2010/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw -O Human_PBMC_Li2010.bw &
wget http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw -O Human_PBMC_Heyn2012.meth.bw &
wget http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-100yr/tracks_hg19/Human_CD4T-100yr.meth.bw -O Human_CD4T-100yr.meth.bw
wget http://smithdata.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-Newborn/tracks_hg19/Human_CD4T-Newborn.meth.bw -O 
wget http://smithdata.usc.edu/methbase/data/Pacis-Human-2015/Human_DendriticCell/tracks_hg19/Human_DendriticCell.meth.bw -O Human_DendriticCell_hg19.meth.bw

cat Human_PBMC_*wig > HUman_PBMC_hg19.wig
grep -v '#' HUman_PBMC_hg19.wig > HUman_PBMC_hg19.wig.bed
perl wigAverage.pl HUman_PBMC_hg19.wig.bed > HUman_PBMC_hg19.wig.bed.uni
bedtools sort -i HUman_PBMC_hg19.wig.bed.uni > HUman_PBMC_hg19.wig.bed.uni.sort
echo 'track type=bedGraph name="PBMC" description="PBMC" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > HUman_PBMC_hg19.wig.bed.uni.sort.UCSC
cat  HUman_PBMC_hg19.wig.bed.uni.sort >> HUman_PBMC_hg19.wig.bed.uni.sort.UCSC
head HUman_PBMC_hg19.wig.bed.uni.sort.UCSC
gzip HUman_PBMC_hg19.wig.bed.uni.sort.UCSC

cat Human_CD4T*wig > HUman_CD4T_hg19.wig
grep -v '#' HUman_CD4T_hg19.wig > HUman_CD4T_hg19.wig.bed
perl wigAverage.pl HUman_CD4T_hg19.wig.bed > HUman_CD4T_hg19.wig.bed.uni
bedtools sort -i HUman_CD4T_hg19.wig.bed.uni > HUman_CD4T_hg19.wig.bed.uni.sort
echo 'track type=bedGraph name="CD4+" description="CD4+" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > HUman_CD4T_hg19.wig.bed.uni.sort.UCSC
cat  HUman_CD4T_hg19.wig.bed.uni.sort >> HUman_CD4T_hg19.wig.bed.uni.sort.UCSC
head HUman_CD4T_hg19.wig.bed.uni.sort.UCSC
gzip -f HUman_CD4T_hg19.wig.bed.uni.sort.UCSC

SID="HUman_B_Cell_hg19"
cat Human_BCell*wig > $SID.wig
grep -v '#' $SID.wig > $SID.wig.bed
perl wigAverage.pl $SID.wig.bed > $SID.wig.bed.uni
bedtools sort -i $SID.wig.bed.uni > $SID.wig.bed.uni.sort
echo 'track type=bedGraph name="Human_B_Cell" description="Human_B_Cell" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > $SID.wig.bed.uni.sort.UCSC
cat  $SID.wig.bed.uni.sort >> $SID.wig.bed.uni.sort.UCSC
head $SID.wig.bed.uni.sort.UCSC
gzip -f $SID.wig.bed.uni.sort.UCSC

SID="Human_Neutrophil_hg19"
cat Human_Neut*bw.wig > $SID.wig
grep -v '#' $SID.wig > $SID.wig.bed
perl wigAverage.pl $SID.wig.bed > $SID.wig.bed.uni
bedtools sort -i $SID.wig.bed.uni > $SID.wig.bed.uni.sort
echo 'track type=bedGraph name="Human_Neutrophil" description="Human_Neutrophil" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > $SID.wig.bed.uni.sort.UCSC
cat  $SID.wig.bed.uni.sort >> $SID.wig.bed.uni.sort.UCSC
head $SID.wig.bed.uni.sort.UCSC
gzip -f $SID.wig.bed.uni.sort.UCSC


# Colon cancer (CRC)
wget http://smithdata.usc.edu/methbase/data/Berman-Human-2012/Human_ColonNormal/tracks_hg19/Human_ColonNormal.meth.bw -O Human_ColonNormal_Berman-Human-2012.meth.bw
wget http://smithdata.usc.edu/methbase/data/Hansen-Human-2011/Human_ColonicMucosa/tracks_hg19/Human_ColonicMucosa.meth.bw -O Human_ColonicMucosa_Hansen-Human-2011.meth.bw
wget http://smithdata.usc.edu/methbase/data/Hansen-Human-2011/Human_AdenoPolyp/tracks_hg19/Human_AdenoPolyp.meth.bw -O Human_AdenoPolyp_Hansen-Human-2011.meth.bw

wget http://smithdata.usc.edu/methbase/data/Berman-Human-2012/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw -O Human_ColonCancer_CRC_Berman-Human-2012.meth.bw  
wget http://smithdata.usc.edu/methbase/data/Ziller-Human-2013/Human_Colon_Tumor_Primary/tracks_hg19/Human_Colon_Tumor_Primary.meth.bw -O Human_Colon_Tumor_Primary_CRC_Ziller-Human-2013.meth.bw
wget http://smithdata.usc.edu/methbase/data/Hansen-Human-2011/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw -O Human_ColonCancer_CRC_Hansen-Human-2011.meth.bw
wget http://smithdata.usc.edu/methbase/data/Akalin-Human-2012/Human_HCT116/tracks_hg19/Human_HCT116.meth.bw -O Human_HCT116_CRC_Akalin-Human-2012.meth.bw

cat *_CRC_*.wig > Human_CCRC_Cancer.wig
grep -v '#' Human_CCRC_Cancer.wig > Human_CCRC_Cancer.wig.bed
perl wigAverage.pl Human_CCRC_Cancer.wig.bed > Human_CCRC_Cancer.wig.bed.uni
bedtools sort -i Human_CCRC_Cancer.wig.bed.uni > Human_CCRC_Cancer.wig.bed.uni.sort
echo 'track type=bedGraph name="CRC" description="CRC" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > Human_CCRC_Cancer.wig.bed.uni.sort.UCSC
cat  Human_CCRC_Cancer.wig.bed.uni.sort >> Human_CCRC_Cancer.wig.bed.uni.sort.UCSC
head Human_CCRC_Cancer.wig.bed.uni.sort.UCSC
gzip Human_CCRC_Cancer.wig.bed.uni.sort.UCSC

cat *_CRC_*.wig > Human_CCRC_Cancer.wig
grep -v '#' Human_CCRC_Cancer.wig > Human_CCRC_Cancer.wig.bed
perl wigAverage.pl Human_CCRC_Cancer.wig.bed > Human_CCRC_Cancer.wig.bed.uni
bedtools sort -i Human_CCRC_Cancer.wig.bed.uni > Human_CCRC_Cancer.wig.bed.uni.sort
echo 'track type=bedGraph name="CRC" description="CRC" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > Human_CCRC_Cancer.wig.bed.uni.sort.UCSC
cat  Human_CCRC_Cancer.wig.bed.uni.sort >> Human_CCRC_Cancer.wig.bed.uni.sort.UCSC
head Human_CCRC_Cancer.wig.bed.uni.sort.UCSC
gzip Human_CCRC_Cancer.wig.bed.uni.sort.UCSC

# Lung and Lung cancer
wget http://smithdata.usc.edu/methbase/data/Roadmap-Human-2015/Human_Lung/tracks_hg19/Human_Lung.meth.bw -O Human_Lung_Roadmap_hg19.meth.bw
bigWigToWig Human_Lung_Roadmap_hg19.meth.bw Human_Lung_Roadmap_hg19.meth.bw.wig
SID="Human_Lung_Roadmap_hg19"
cat Human_Lung_Roadmap*bw.wig > $SID.wig
grep -v '#' $SID.wig > $SID.wig.bed
perl wigAverage.pl $SID.wig.bed > $SID.wig.bed.uni
bedtools sort -i $SID.wig.bed.uni > $SID.wig.bed.uni.sort
echo 'track type=bedGraph name="Human_Neutrophil" description="Human_Neutrophil" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > $SID.wig.bed.uni.sort.UCSC
cat  $SID.wig.bed.uni.sort >> $SID.wig.bed.uni.sort.UCSC
head $SID.wig.bed.uni.sort.UCSC
gzip -f $SID.wig.bed.uni.sort.UCSC


# Breast and breast cancer (BRCA)
wget http://smithdata.usc.edu/methbase/data/Lin-Human-2015/Human_HealthyBreast/tracks_hg19/Human_HealthyBreast.meth.bw -O Human_HealthyBreast_Lin-Human-2015.meth.bw

wget http://smithdata.usc.edu/methbase/data/Grimmer-Human-2014/Human_MCF7/tracks_hg19/Human_MCF7.meth.bw -O Human_MCF7_BRCA_Grimmer-Human-2014.meth.bw
wget http://smithdata.usc.edu/methbase/data/Hon-Human-2012/Human_HCC1954/tracks_hg19/Human_HCC1954.meth.bw -O Human_HCC1954_BRCA_Hon-Human-2012.meth.bw
wget http://smithdata.usc.edu/methbase/data/Lin-Human-2015/Human_BreastTumor089/tracks_hg19/Human_BreastTumor089.meth.bw -O Human_BreastTumor089_BRCA_Lin-Human-2015.meth.bw
wget http://smithdata.usc.edu/methbase/data/Lin-Human-2015/Human_BreastTumor126/tracks_hg19/Human_BreastTumor126.meth.bw -O  Human_BreastTumor126_BRCA_Lin-Human-2015.meth.bw
wget http://smithdata.usc.edu/methbase/data/Lin-Human-2015/Human_BreastTumor198/tracks_hg19/Human_BreastTumor198.meth.bw -O Human_BreastTumor198_BRCA_Lin-Human-2015.meth.bw
wget http://smithdata.usc.edu/methbase/data/Menafra-Human-2014/Human_MCF7/tracks_hg19/Human_MCF7.meth.bw -O Human_MCF7_BRCA_Menafra-Human-2014.meth.bw

cat *_MCF7_*.wig > Human_BRCA.wig
grep -v '#' Human_BRCA.wig > Human_BRCA.meth.hg19
bedtools sort -i Human_BRCA.meth.hg19 > Human_BRCA.meth.hg19.wig
echo 'track type=bedGraph name="BRCA" description="BRCA" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > Human_BRCA.meth.hg19.wig.UCSC
cat  BRCA.bedgraph.uni.meth.sort >> Human_BRCA.meth.hg19.wig.UCSC
head Human_BRCA.meth.hg19.wig.UCSC
gzip Human_BRCA.meth.hg19.wig.UCSC

# pancrease cancer
for i in {1..12}
do
wget http://smithdata.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw &
done
wget http://smithdata.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas1/tracks_hg19/Human_NormalPancreas1.meth.bw
wget http://smithdata.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas2/tracks_hg19/Human_NormalPancreas2.meth.bw

cat *_PancreaticCancer*.wig > PANcreaticCancer.wig
grep -v '#' PANcreaticCancer.wig > PANcreaticCancer.wig.bed
perl wigAverage.pl PANcreaticCancer.wig.bed > PANcreaticCancer.wig.bed.uni
bedtools sort -i PANcreaticCancer.wig.bed.uni > PANcreaticCancer.wig.bed.uni.sort
echo 'track type=bedGraph name="PancreaticCancer" description="_PancreaticCancer" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > PANcreaticCancer.wig.bed.uni.sort.UCSC
cat  PANcreaticCancer.wig.bed.uni.sort >> PANcreaticCancer.wig.bed.uni.sort.UCSC
head PANcreaticCancer.wig.bed.uni.sort.UCSC
gzip PANcreaticCancer.wig.bed.uni.sort.UCSC

cat *_NormalPancreas*.wig > NOrmalPancreas.wig
grep -v '#' NOrmalPancreas.wig > NOrmalPancreas.wig.bed
perl wigAverage.pl NOrmalPancreas.wig.bed > NOrmalPancreas.wig.bed.uni
bedtools sort -i NOrmalPancreas.wig.bed.uni > NOrmalPancreas.wig.bed.uni.sort
echo 'track type=bedGraph name="Normal_Pancreas" description="Normal_Pancreas" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > NOrmalPancreas.wig.bed.uni.sort.UCSC
cat  NOrmalPancreas.wig.bed.uni.sort >> NOrmalPancreas.wig.bed.uni.sort.UCSC
head NOrmalPancreas.wig.bed.uni.sort.UCSC
gzip -f NOrmalPancreas.wig.bed.uni.sort.UCSC

# AML 
wget http://smithdata.usc.edu/methbase/data/Akalin-Human-2012/Human_AML/tracks_hg19/Human_AML.meth.bw -O Human_AML_Akalin-Human-2012.meth.bw
wget http://smithdata.usc.edu/methbase/data/Lund-Human-2014/Human_AML3-Control/tracks_hg19/Human_AML3-Control.meth.bw -O Human_AML3_Lund-Human-2014.meth.bw
cat  Human_AML*.wig > HUman_AML_WGBS.wig
grep -v '#' HUman_AML_WGBS.wig > HUman_AML_WGBS.wig.bed
perl wigAverage.pl HUman_AML_WGBS.wig.bed > HUman_AML_WGBS.wig.bed.uni
bedtools sort -i HUman_AML_WGBS.wig.bed.uni > HUman_AML_WGBS.wig.bed.uni.sort
echo 'track type=bedGraph name="AML" description="AML" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > HUman_AML_WGBS.wig.bed.uni.sort.UCSC
cat  HUman_AML_WGBS.wig.bed.uni.sort >> HUman_AML_WGBS.wig.bed.uni.sort.UCSC
head HUman_AML_WGBS.wig.bed.uni.sort.UCSC
gzip -f HUman_AML_WGBS.wig.bed.uni.sort.UCSC

# Human_B_Lymphocyte
wget http://smithdata.usc.edu/methbase/data/Ball-Human-2009/Human_BLymphocyte/tracks_hg19/Human_BLymphocyte.meth.bw -O  Human_BALL_Ball_Human-2009.meth.bw
wget http://smithdata.usc.edu/methbase/data/Lee-Human-2015/Human_BALLhyperdiploid/tracks_hg19/Human_BALLhyperdiploid.meth.bw -O Human_BALL_hyperdiploid_Lee-Human-2015.meth.bw
wget http://smithdata.usc.edu/methbase/data/Lee-Human-2015/Human_BALLETV6RUNX1/tracks_hg19/Human_BALLETV6RUNX1.meth.bw -O Human_BALL_ETV6RUNX1-Human-2015.meth.bw
cat  Human_BALL_*.wig > HUman_BALL_WGBS.wig
grep -v '#' HUman_BALL_WGBS.wig > HUman_BALL_WGBS.wig.bed
perl wigAverage.pl HUman_BALL_WGBS.wig.bed > HUman_BALL_WGBS.wig.bed.uni
bedtools sort -i HUman_BALL_WGBS.wig.bed.uni > HUman_BALL_WGBS.wig.bed.uni.sort
echo 'track type=bedGraph name="B-ALL" description="B-ALL" visibility=full color=227,207,87 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > HUman_BALL_WGBS.wig.bed.uni.sort.UCSC
cat  HUman_BALL_WGBS.wig.bed.uni.sort >> HUman_BALL_WGBS.wig.bed.uni.sort.UCSC
head HUman_BALL_WGBS.wig.bed.uni.sort.UCSC
gzip HUman_BALL_WGBS.wig.bed.uni.sort.UCSC

# Non-small cell lung cancer cells, IMR-90 (normal lung cell line fibroblast), A549
wget http://smithdata.usc.edu/methbase/data/Brocks-Human-2017/Human_DMSO/tracks_hg19/Human_DMSO.meth.bw -O Human_H1299_Brocks-Human-2017.meth.bw
bigWigToWig Human_H1299_Brocks-Human-2017.meth.bw Human_H1299_Brocks-Human-2017.meth.bw.wig
grep -v '#' Human_H1299_Brocks-Human-2017.meth.bw.wig > Human_H1299_Brocks-Human-2017.meth.wig &
echo 'track type=bedGraph name="H1299" description="H1299" visibility=full color=0,100,200 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > Human_H1299_Brocks-Human-2017.meth.wig.UCSC
cat  Human_H1299_Brocks-Human-2017.meth.wig >> Human_H1299_Brocks-Human-2017.meth.wig.UCSC
head Human_H1299_Brocks-Human-2017.meth.wig.UCSC
gzip Human_H1299_Brocks-Human-2017.meth.wig.UCSC

awk '{print $1,$2,$3,$5/1000}' GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.bed > GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.wig
head GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.wig
echo 'track type=bedGraph name="HepG2" description="HepG2" visibility=full color=0,100,200 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.UCSC.wig
cat  GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.wig >> GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.UCSC.wig
head GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.UCSC.wig
gzip GSM1204464_BiSeq_cpgMethylation_BioSam_1502_IMR_90_304071.BiSeq.UCSC.wig


# prostate cancer
wget http://smithdata.usc.edu/methbase/data/Pidsley-Human-2016/Human_LNCaP/tracks_hg19/Human_LNCaP.meth.bw -O Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.bw
grep -v '#' Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.bw.wig > Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.wig &
echo 'track type=bedGraph name="LNCaP_ProstateCancer" description="LNCaP_ProstateCancer" visibility=full color=0,100,200 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full' > Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.wig.UCSC
cat  Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.wig >> Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.wig.UCSC
head Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.wig.UCSC
gzip Human_LNCaP_ProstateCancer_Pidsley-Human-2016.meth.wig.UCSC

# liver cancer, 
wget -r https://genome.ucsc.edu/gbdb/methylome/public/Ziller-Human-2013/Human_HepG2/tracks_hg19/
track type=bedGraph name="HepG2" description="HepG2" visibility=full color=0,100,200 altColor=200,100,0 priority=20 maxHeightPixels=128:64:32 visibility=full
awk '{print $1,$2,$3,$5/1000}' GSM1204463_BiSeq_cpgMethylation_BioSam_1500_HepG2_304072.BiSeq.bed > GSM1204463_BiSeq_cpgMethylation_BioSam_1500_HepG2_304072.BiSeq.wig
gzip GSM1204463_BiSeq_cpgMethylation_BioSam_1500_HepG2_304072.BiSeq.wig


# K562

# cancer associated fibroblasts and non-malignant prostate fibroblasts (NPFs): GSE85609

# WGBS of esophagus squamous epithelium/esophagus muscularis mucosa


# Lung and lung cancer (450K)
GSE39279: 444 lung cancer 450K array
GSE92843: 3 lung cancer cell lines (A549, A427 and H322), normal bronchial ephitelial cells (NHBEC)
GSE63704: 17 LC, 37 idiopathic lung fibrosis,32 patients suffering from chronic obstructive pulmonary disease and 43 DNA samples derived from healthy-lungs 
GSE63940: 36 lung adenocarcinoma cell lines was bisulfite treated and analyzed on the Illumina Infinium HumanMethylation27K 
GSE85566: isolated airway epithelial cells of asthmatics and non-asthmatics
# lung cancer (27K)
GSE62950: 28 pair LC-Normal(27K and mRNA)
GSE32867: 59 pair LC-Normal(27K and mRNA)
# lung cancer (27K)


wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39279/suppl/GSE39279_RAW.tar
mkdir GSE39279
mv GSE39279_RAW.tar GSE39279



perl VEPmerge.pl > BLM-VEP.hg19.bed

for i in `ls *multianno.txt.bed`
do
bedtools intersect -wo -a $i -b BLM-VEP.hg19.bed > $i.vep.bed
done

awk '{print $2,$4}' OFS="\t" BLM-VEP.hg19.txt |  awk -F "[:-]" '{print "chr"$1,$2,$3,$4}' OFS="\t" | grep -v chrLocation | sort -u > BLM-VEP.hg19.bed

data<-read.table("LBM.MutationProfile.txt",head=T,sep="\t",row.names=1)
target<-c("TP53","KRAS","FAT4","STK11","EGFR","KMT2C","CHEK2P2","MIR4436A","BAGE2","BAGE3","BAGE4","BAGE5","LOC102723769","MIR6077","FER1L4","FRG1BP","FRG1DP","LOC100996724","AHNAK2","FRG1BP","HMCN2","TPTE")
input<-data[,na.omit(match(target,colnames(data)))]

perl VEPmerge.pl > BLM-VEP.hg19.bed
perl bed2mutationprofile.pl > LBM.MutationProfile.txt
cp LBM.MutationProfile.txt LBM.MutationProfile.test.txt
perl -p -i -e 's/3_prime_UTR_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/5_prime_UTR_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/coding_sequence_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/downstream_gene_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/frameshift_variant/frameshift_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/inframe_deletion/frameshift_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/intergenic_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/intron_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/mature_miRNA_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/missense_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/NMD_transcript_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/non_coding_transcript_exon_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/non_coding_transcript_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/protein_altering_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/regulatory_region_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_acceptor_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_donor_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_region_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/start_lost/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/stop_gained/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/stop_lost/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/synonymous_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/TF_binding_site_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/upstream_gene_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt




awk '$8 !="." & $8 != "synonymous"{print FILENAME}' *hg19_multianno.txt.bed

grep -v "#" *.vcf | awk '{print $1,$2,$3,$4,%5}' 

awk '{print $1,$2,".",$4,$5,".",FILENAME,"."}' OFS="\t" *.bed > LungBrainMetasis.hg19.vcf
perl -p -i -e 's/chr//g' LungBrainMetasis.hg19.vcf
sort -u LungBrainMetasis.hg19.vcf > LungBrainMetasis.uni.hg19.vcf

dss2bedgraph

mkdir temp
for i in `ls *dss`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job 
echo awk \'\$3\>0{print \$1,\$2-1,\$2,\$4\/\$3}\' OFS=\"\\t\" $i \> $i.bed >> ./temp/$i.job
qsub ./temp/$i.job
done

bedtools sort -i GSE112658.RA.hg19.bedgraph > GSE112658.RA.sort.hg19.bedgraph &
bedtools sort -i GSE112658.OA.hg19.bedgraph > GSE112658.OA.sort.hg19.bedgraph &

cat GSE112658.RA.sort.hg19.bedgraph >> RA.hg19.bedgraph
cat GSE112658.OA.sort.hg19.bedgraph >> OA.hg19.bedgraph

gzip -c RA-FLS.hg19.bedgraph > RA-FLS.hg19.bedgraph.gz &
gzip -c OA-FLS.hg19.bedgraph > OA-FLS.hg19.bedgraph.gz &

bedtools intersect -wo -a RA-Naturecommnunication.DMR.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed > RA-Naturecommnunication.DMR.RefGene.hg19.bed
awk '{print $15"."$3,$1,$2,$3,$15"."$3}' OFS="\t" RA-Naturecommnunication.DMR.RefGene.hg19.bed | sort -u > RA-Naturecommnunication.DMR.RefGene.hg19.uni.bed



Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install pysam and methpipe 
# download the lastest pysam (older version will report error)
# Go to: https://pypi.python.org/pypi/pysam
tar xzvf pysam-0.9.1.4.tar.gz
cd pysam-0.9.1.4
python setup.py build
python setup.py install --user

wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz


Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex




Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\mh450
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-03-09.txt

Set-Location //mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancancer/RNA-seq
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-03-09.txt

Pancancer-TCGA-meth450-8893-clinical.tsv

data1<-read.table("/gpfs/home/guosa/hpc/db/hg19/H3k27ac/VIP.genelist.txt",head=F)
gene1<-unique(data[,1])
data2<-read.table("/gpfs/home/guosa/hpc/GWAS_Catalog/immnue.Gene.txt")
gene2<-names(which(table(data2[,1])>10))

db<-read.table("~/hpc/db/hg19/refGeneV2.hg19.bed",head=F)
newdb<-subset(db,V8=="Exon1"| V8=="Enhancer" | V8=="Promoter")

output1<-newdb[newdb$V6 %in% gene1,]
output2<-newdb[newdb$V6 %in% gene2,]

output<-unique(rbind(output1,output2))
unique(output$V6)
write.table(output,file="VIP.PID.Regulatory.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)

bedtools sort -i VIP.PID.Regulatory.hg19.bed > VIP.PID.Regulatory.hg19.sort.bed
mv VIP.PID.Regulatory.hg19.sort.bed VIP.PID.Regulatory.hg19.bed

bedtools intersect -b VIP.PID.Regulatory.hg19.bed -a wgEncodeRegTfbsClusteredWithCellsV3.bed | sort -u > VIP.TFBS.hg19.bed
bedtools intersect -b VIP.PID.Regulatory.hg19.bed -a wgEncodeRegDnaseClusteredV3.bed | sort -u > VIP.Dnase.hg19.bed
bedtools intersect -a VIP.Dnase.hg19.bed -b VIP.TFBS.hg19.bed | sort -u > VIP.TFBS.Dnase.hg19.bed
bedtools sort -i VIP.TFBS.Dnase.hg19.bed >  VIP.TFBS.Dnase.hg19.sort.bed
bedtools merge -d 500 -i VIP.TFBS.Dnase.hg19.sort.bed  | awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" > VIP.TFBS.Dnase.hg19.merge.bed

bedtools intersect -wo -a VIP.TFBS.Dnase.hg19.merge.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed | grep PTPN

# bw2tab
mkdir temp
for i in `ls *bigWig`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job 
echo bigWigAverageOverBed $i VIP.TFBS.Dnase.hg19.merge.bed $i.tab >> ./temp/$i.job
qsub ./temp/$i.job
done


tab2matrix.pl > Matrix.txt

data<-read.table("Matrix.txt")
data$RowMean<-rowMeans(data)
write.table(data,file="Matrix.RowMean.bed",sep="\t",quote=F,row.names=T,col.names=NA)

awk '{print $1,$23}' OFS="\t" Matrix.RowMean.bed | grep -v Encod > Matrix.RowMean.bid
perl -p -i -e 's/:/\t/' Matrix.RowMean.bid
perl -p -i -e 's/-/\t/' Matrix.RowMean.bid

bedtools intersect -wo -a Matrix.RowMean.bid -b ~/hpc/db/hg19/refGeneV2.hg19.bed > Matrix.RowMean.Symbol.bid

data<-read.table("Matrix.RowMean.Symbol.bid",head=F)
newdata<-subset(data,V11=="Exon1"| V11=="Enhancer" | V11=="Promoter")
dim(newdata)

output<-c()
for(i in unique(data$V10)){
input<-unique(subset(data,V10==i))
new<-unique(input[,1:4])
if(nrow(new)>1){
Max<-which.max(new$V4)
output<-rbind(output,input[Max,])
}else{
output<-rbind(output,input)
}
}
rlt<-unique(output[,1:4])
write.table(rlt,file="VIP.PID.Regulatory.MaxScore-Region.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(output,file="VIP.PID.Regulatory.MaxScore-Gene.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneH1hescH3k4me3StdSig.bigWig  
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneHsmmH3k4me3StdSig.bigWig     
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneHuvecH3k4me3StdSig.bigWig    
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig    
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneNhekH3k4me3StdSig.bigWig    
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/wgEncodeBroadHistoneNhlfH3k4me3StdSig.bigWig    
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/md5sum.txt

wget -r http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k27ac/

wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k27ac/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me1/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/

tabix -p vcf /gpfs/home/guosa/hpc/db/hg19/allSNP151.vcf.gz -B GWAS.hg19.txt > GWAS.allSNP151.vcf

for i in `ls *dss`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job
echo awk \'\$3\>4{print \$1,\$2-1,\$2,\$4\/\$3}\' OFS\=\"\\t\" $i \> $i.bedgraph >> ./temp/$i.job
qsub ./temp/$i.job
done

for i in `ls *.bedgraph`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job
echo bedtools sort -i $i \> $i.sort >> ./temp/$i.job
echo bedGraphToBigWig $i.sort /gpfs/home/guosa/hpc/db/hg19/hg19.chrom.sizes $i.bw >> ./temp/$i.job
qsub ./temp/$i.job
done



for i in `ls *.bedgraph`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedmethyl >> $i.job
echo j=\"${i//.bed.bedgraph/}\" >> $i.job
echo wigToBigWig $i ~/hpc/db/hg38/hg38.chrome.sizes \$j.bw >> $i.job
# qsub $i.job
done



Introduce an novel disease gene mapping approach and show an example to identify new hemochromatosis gene, FGF6, with PMRP dataset. Understanding the mechanisms of iron metabolism dysfunction can aid the development of therapeutic interventions for patients with hemochromatosis and possibly anemia



chr10   128946873       128947161       DOCK1
chr10   52835666        52835788        PRKG1
chr1    207497041       207497487       CD55
chr1    214601110       214601543       PTPN14
chr14   69065759        69065829        RAD51B
chr18   463850  464312  COLEC12
chr2    113885018       113885168       IL1RN
chr21   36421047        36421266        RUNX1
chr3    57584440        57584495        PDE12
chr5    56113778        56114029        MAP3K1
chr7    19152486        19152738        TWIST1
chr7    19154311        19154497        TWIST1


bedtools intersect -wo -a WGBS.hg19.bed -b 2133ImmuGene.hg19.bed | awk '{print $1,$2,$3,$10}' OFS="\t" > CandidateList1.hg19.bed
cat CandidateList1.hg19.bed WGBS.GREAT.Gene.hg19.bed  | awk '{print $1,$2,$3}' OFS="\t" | sort -u | bedtools sort > CandidateList2.hg19.bed
bedtools closest -d -wo -a CandidateList2.hg19.bed -b GWAS-SNP-359.GRCH37.bed > CandidateList3.hg19.bed
bedtools closest -d -wo -a CandidateList3.hg19.bed -b ~/hpc/db/hg19/refGene.hg19.bed | awk '{print $1,$2,$3,$7,$8,$14}' OFS="\t" | sort -u > CandidateList4.hg19.bed

library(DSS)
require(bsseq)
path <- file.path(system.file(package="DSS"), "extdata")
dat1.1 <- read.table(file.path(path, "cond1_1.txt"), header=TRUE)
dat1.2 <- read.table(file.path(path, "cond1_2.txt"), header=TRUE)
dat2.1 <- read.table(file.path(path, "cond2_1.txt"), header=TRUE)
dat2.2 <- read.table(file.path(path, "cond2_2.txt"), header=TRUE)
BSobj <- makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2),c("C1","C2", "N1", "N2") )[1:1000,]
BSobj
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))
dmlTest.sm <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"), smoothing=TRUE)
dmls <- callDML(dmlTest, p.threshold=0.001)
dmls2 <- callDML(dmlTest, delta=0.1, p.threshold=0.001)
dmrs <- callDMR(dmlTest, p.threshold=0.01)
dmrs2 <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
showOneDMR(dmrs[1,], BSobj)


/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication




Transfer miR_Target database to bed format
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/All_Target_Locations.hg19.bed
for i in `ls *.bed`; 
do 
perl format.pl $i > $i.txt & 
done
mirTarget2bed.pl

#################################################################################################################################
grep PTPN ~/hpc/db/hg19/refGene.hg19.bed >> RAcandidate.hg19.bed
grep PADI ~/hpc/db/hg19/refGene.hg19.bed >> RAcandidate.hg19.bed
perl -lane "{print if /\s+IL\d+/}"  ~/hpc/db/hg19/refGene.hg19.bed >>  RAcandidate.hg19.bed
perl -lane "{print if /\s+MUC\d+/}"  ~/hpc/db/hg19/refGene.hg19.bed >>  RAcandidate.hg19.bed
perl -p -i -e "s/chr//" RAcandidate.hg19.bed
sort -u RAcandidate.hg19.bed > RAcandidate.uni.hg19.bed
mv RAcandidate.uni.hg19.bed  RAcandidate.hg19.bed
panel="refGene"
mkdir temp
perl -p -i -e "{s/chr//}" $panel.hg19.bed
## Function Variants in Exom Regions
for i in {1..22} X Y
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
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.01 \& INFO/AF_eas\<0.99 \& \(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' -R $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
#################################################################################################################################


#################################################################################################################################
grep PTPN ~/hpc/db/hg19/refGene.hg19.bed >> RAcandidate.hg19.bed
grep PADI ~/hpc/db/hg19/refGene.hg19.bed >> RAcandidate.hg19.bed
perl -lane "{print if /\s+IL\d+/}"  ~/hpc/db/hg19/refGene.hg19.bed >>  RAcandidate.hg19.bed
perl -lane "{print if /\s+MUC\d+/}"  ~/hpc/db/hg19/refGene.hg19.bed >>  RAcandidate.hg19.bed
perl -p -i -e "s/chr//" RAcandidate.hg19.bed
sort -u RAcandidate.hg19.bed > RAcandidate.uni.hg19.bed
mv RAcandidate.uni.hg19.bed  RAcandidate.hg19.bed


panel="MRCI"
mkdir temp
perl -p -i -e "{s/chr//}" $panel.hg19.bed

## Function Variants in Exom Regions
for i in {1..22} X Y
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
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO/AF_eas\<0.999 \& \(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' -R $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
#################################################################################################################################


# Since I want to see the haplotype, I increased the AF
#####################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/UTR3SNP
perl -p -i -e 's/chr//g' UTR3-miRNA.hg19.bed
panel="UTR3miRNAsNP"
mkdir temp
## Function Variants in Genome Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \# bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \# tabix -p vcf gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.01 \&INFO/AF_eas\<0.99\' -T $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *genomes*rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf | awk '{print "chr"$1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
wc -l gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
#####################################################################

awk '{print "chr"$2,$3-1,$3,$4}' OFS="\t"  ASACHIA_rsID.txt | grep -v 'CHROM' | sort -u | bedtools sort > ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.4728.rec.ReactomePathWay.immnueGene.hg19.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a GWAS-Meta-128-SNPs.20190208.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed
bedtools intersect -v -a 48-SNPs-Genetic-Risk-Score-EUR.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.HLA.hg19.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 1375.gnomad.exomes.r2.1.sites.rec.RA-GWAS-Cytoband.hg19.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a GWAS-immnue-3325_SNP.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.40KGWASAID.merge.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.InnateDB.merge.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.genomes.r2.1.sites.rec.innateDbUTR3.merge.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 15446.MRCI.ASA.eQTL.hg19.MAF0.001.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.eQTL.set2.hg19.vcf.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a cpgSNPisland.AID.GWAS.SNP.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.RAcandidate.hg19.vcf.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a T325.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a hsa-miRNALD.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 63.miRNA.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.8501.rec.GHRA_ASA.hg19.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.genomes.r2.1.sites.rec.GWASCatalog.merge.vcf.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.genomes.r2.1.sites.rec.ExomGenome.merge.vcf.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

grep -v 'CHROM'  ASA.hg19.bed | sort -u > ASA.hg19.P2.bed
wc -l ASA.hg19.P2.bed
awk '{print "chr"$2,$3-1,$3,$4}' OFS="\t" ASACHIA_rsID.txt | grep -v CHROM | sort -u > ASA.hg19.P1.bed
wc -l ASA.hg19.P1.bed
bedtools intersect -v -a ASA.hg19.P2.bed  -b ASA.hg19.P1.bed| awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > MRCI.hg19.bed
wc -l MRCI.hg19.bed




bedtools intersect -v -a gnomad.genomes.r2.1.sites.rec.ExomGenome.merge.vcf.hg19.bed -b ../ASA.hg19.P1.bed| awk '{print $1,$2,$3,$4}' OFS="\t" | wc -l
bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.TotalCandidateGene.hg19.vcf.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed






cp /gpfs/home/guosa/hpc/rheumatology/RA/ASA/cpgSNPisland/cpgSNPisland.AID.GWAS.SNP.hg19.bed ./

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/2665-GWAS-LD0.8-ASN.hg19.rsList.txt
awk '{print "chr"$1,$2-1,$2,$3,$4,$5}' OFS="\t" 2665-GWAS-LD0.8-ASN.hg19.rsList.txt > 2665-GWAS-LD0.8-ASN.hg19.rsList.hg19.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/GWAS-immnue-3325_SNP.hg19.bed

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/gnomad.exomes.r2.1.sites.rec.40KGWASAID.merge.vcf.bed
awk '{print "chr"$1,$2-1,$2,$3,$4,$5}' OFS="\t" gnomad.exomes.r2.1.sites.rec.40KGWASAID.merge.vcf.bed > gnomad.exomes.r2.1.sites.rec.40KGWASAID.merge.hg19.bed

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf.bed
awk '{print "chr"$1,$2-1,$2,$3,$4,$5}' OFS="\t" gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf.bed > gnomad.exomes.r2.1.sites.rec.InnateDB.merge.hg19.bed

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/gnomad.genomes.r2.1.sites.rec.innateDbUTR3.merge.vcf.bed
awk '{print "chr"$1,$2-1,$2,$3,$4,$5}' OFS="\t" gnomad.genomes.r2.1.sites.rec.innateDbUTR3.merge.vcf.bed > gnomad.genomes.r2.1.sites.rec.innateDbUTR3.merge.hg19.bed

1375.gnomad.exomes.r2.1.sites.rec.RA-GWAS-Cytoband.hg19.vcf.bed



# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz



cd 
cp ~/hpc/db/hg19/cpgCSNP.hg19.bed ./
awk '{print "chr"$1,$2-1,$2,$3}' OFS="\t" gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.bed > gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.hg19.bed
bedtools intersect -wo -a gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.hg19.bed -b cpgCSNP.hg19.bed | sort -u >  gnomad.genomes.eQTL.cpgSNP.uni.hg19.bed

bedtools intersect -wo -a /gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP/gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed | sort -u >  gnomad.genomes.eQTL.TFBS.uni.hg19.bed


bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.hg19.bed -b  | sort -u > GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed






#####################################################################
## CpG-SNP counts in 2000bp Genomic Windows
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/cpgSNPisland
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=16 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >>./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job
echo bedtools intersect -c -a ~/hpc/db/hg19/window2000/hg19.chr$i.win2K.bed -b ~/hpc/db/hg19/cpgCSNP.hg19.bed \>chr$i.count.txt >> ./temp/$i.job
qsub ./temp/$i.job
done
#########################################################################

/gpfs/home/guosa/hpc/db/hg38/fa/chroms/cpgSNP.hg38.bed

/gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP


for i in {21..43}
do
qdel 5065$i.bright
done

# Since I want to see the haplotype, I increased the AF
#####################################################################
cd /gpfs/home/guosa/hpc/cpgSNP/Asian
cp /gpfs/home/guosa/hpc/db/hg19/cpgSNP.hg19.bed ./
perl -p -i -e 's/chr//g' cpgSNP.hg19.bed
panel="cpgSNPs"
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
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.3 \&INFO/AF_eas\<0.7\' -T $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
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


#####################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP
perl -p -i -e 's/chr//g' eQTL.hg19.bed

panel="eQTL"
bed="eQTL.hg19.bed"
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


awk '{print "dbsnp",$1}' OFS="\t" vip.eqtl.snp.txt | sort -u > eQTL.SnpnexusInput.txt



wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/EyreS_23143596_GCST005569/eyre_2012_23143596_ra_efo0000685_1_ichip.sumstats.tsv.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/EyreS_23143596_GCST005569/harmonised/23143596-GCST005569-EFO_0000685-Build37.f.tsv.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/EyreS_23143596_GCST005569/harmonised/23143596-GCST005569-EFO_0000685.h.tsv.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/EyreS_23143596_GCST005569/harmonised/readme.txt


/gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/GTEx_Analysis_v7_eQTL/snp2bed



grep -fF rsIDs.txt /gpfs/home/guosa/hpc/db/hg19/allSNP151.hg19.bed > rsIDs.hg19.bed


for R in rs9955526 rs9955995 rs9958080 rs9958149  rs9958799
do
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -Dhg19 -N -e "select concat(chrom,':',chromStart+1,'-',chromEnd) from snp151 where name='$R';" | xargs tabix my.vcf.gz
done


wget -qO- ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz | gunzip -c | convert2bed --input=vcf --output=bed --sort-tmpdir=${PWD} - > hg19.snp151.bed

bcftools view gnomad.exomes.r2.1.sites.vcf.bgz | grep -v '#' | awk 'print $1,$2-1,$2,$3,$4,$5' OFS="\t" > gnomad.hg19.bed
grep -F rs554008981 gnomad.hg19.bed
grep -fF rsIDs.txt gnomad.hg19.bed > rsIDs.hg19.bed



LC_ALL=C
wget -qO- https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz \
   | gunzip -c \
   | vcf2bed --sort-tmpdir=${PWD} - \
   | awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $4 }' \
   | sort -k4,4 \
   > hg19.dbSNP151.sRsID.hg19.bed

   
 

split -l 100 snp.txt
mkdir temp
for i in `ls x*`
do
echo \#PBS -N $i  > A.$i.job
echo \#PBS -l nodes=1:ppn=1 >> A.$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> A.$i.job
echo \#PBS -m abe  >> A.$i.job
echo \#PBS -o $(pwd)/temp/ >>A.$i.job
echo \#PBS -e $(pwd)/temp/ >>A.$i.job
echo cd $(pwd) >> A.$i.job
echo 'bcftools view -i '\'ID=@$i\'' ~/hpc/db/hg19/gnomad.genomes.r2.1.sites.vcf.bgz | grep -v '\'#\'' | awk '\'{print \$1,\$2-1,\$2,\$3,\$4,\$5}\'' OFS="\t"' \>\> output.$i >> A.$i.job
qsub A.$i.job
done



wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL

qval_threshold=0.05
data1<-subset(read.table("Whole_Blood.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data2<-subset(read.table("Liver.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data3<-subset(read.table("Small_Intestine_Terminal_Ileum.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data4<-subset(read.table("Stomach.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(as.character(data1[,19]),as.character(data2[,19]),as.character(data3[,19]),as.character(data4[,19]))
length(table(eqtl))
vip.eqtl.snp<-names(table(eqtl))
write.table(vip.eqtl.snp,file="vip.eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)


##################################################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP
panel="eQTL"
bed="eQTL.hg19.bed"
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
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \&INFO/AF_eas\<0.999\' -R $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
##################################################################################################

gnomad.genomes.r2.1.sites.chr11.rec.innateDbUTR3.sort.rmdup.biallelic.vcf.bgz


wget http://www.targetscan.org/vert_72/vert_72_data_download/Predicted_Target_Locations.default_predictions.hg19.bed.zip
unzip Predicted_Target_Locations.default_predictions.hg19.bed.zip

bcftools view -i 'ID=@InnateDB.UTR3.rsid.txt' 

cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/innateDB
/home/local/MFLDCLIN/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr19.rec.vcf.bgz


InnateDB.UTR3.rsid.txt


bcftools view -i 'ID=@InnateDB.UTR3.rsid.txt' /gpfs/home/guosa/hpc/db/hg19/All_20180423.vcf.gz | grep -v '#' | awk '{print $1,$2,$3,$4,$5}' OFS="\t"


wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz.tbi &
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz &
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz.tbi &

wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

awk '{print $4"\t"$5}' miRNA_targets_gain_by_SNPs_in_seed_regions.txt | sort -u | grep rs > miRNA_seed_SNPs.list.txt
awk '{print $4"\t"$5}' miRNA_targets_loss_by_SNPs_in_seed_regions.txt | sort -u | grep rs >> miRNA_seed_SNPs.list.txt
awk '{print $5}' miRNA_targets_gain_by_SNPs_in_seed_regions.txt | sort -u | grep rs > miRNA_seed_SNPs.snponly.list.txt
awk '{print $5}' miRNA_targets_loss_by_SNPs_in_seed_regions.txt | sort -u | grep rs >> miRNA_seed_SNPs.snponly.list.txt

#####################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/innateDB/utrSNP
perl -p -i -e 's/chr//g' innateDbUTR3.hg19.bed

panel="innateDbUTR3"
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

# 2019-02-11
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
VIP.biallelic.MFS.sh

for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
qsub $i.job
done

#########################################################################################################
panel="TotalCandidateGene"
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
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \&INFO/AF_eas\<0.999 \& \(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"missense\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' -R $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
#########################################################################################################

 
 my $cytoband=shift @ARGV;
 open F,$cytoband;
 while(<F>){
 if(/chr(\d+)/){
 chomp;
 my ($chr,$start,$end,$cytoband,$tp)=split/\s+/;
 print "$chr\t$start\t$end\t$1$cytoband\t$tp\n";
 }
 }
 
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


wget -q https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.hg19.bed
grep -v "[;|#|*|SNP|HLA]" immuneGWAS.hg19.bed | grep rs | awk '$2>0{print $1"\t"$2-250000"\t"$3+250000"\t"$6}' > immuneGWAS.hg19.uni.bed
bedtools intersect -wa -a ~/hpc/db/hg19/refGene.hg19.bed -b immuneGWAS.hg19.uni.bed > immuneGWAS.Gene.hg19.bed
awk '{print $5}' immuneGWAS.Gene.hg19.bed | sort -u | wc -l

panel="ReactomePathWay.immnueGene"

wget -q https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/Epigene.hg19.bed
wget -q https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/VIP.hg19.bed
wget -q https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/InnateDB.hg19.bed
wget -q https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.Gene.hg19.bed

grep PTPN ~/hpc/db/hg19/refGene.hg19.bed > $panel.hg19.bed
grep PADI ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed
grep MUC ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed
grep HLA ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed
for i in TLR8 UBASH3A IRF8 CREM PAX5 BC017643 KLRC3 MEF2C FAIM3 CXCR4 ID2 IL2 FOXP3 CD55
do
grep $i ~/hpc/db/hg19/refGene.hg19.bed >> $panel.hg19.bed
done

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
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed

bedtools intersect -c -a ~/hpc/db/hg19/refGene.hg19.bed -b gnomad.exomes.r2.1.sites.rec.GHRA_ASA.hg19.vcf.bed | awk '$6>0'


### Fine-mapping to GWAS 
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/immuneGWAS.hg19.bed

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





 
\\mcrfnas2\bigdata\Genetic\Projects\shg047\hemochromatosis\2LOF
#!/bin/bash 
for i in `seq $4 250 $3`
do
j=$((i+500))
tabix -h -p vcf ALL.chr$2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $2:$i-$j > ./hap500/temp.chr$2.$1.$i.vcf
bcftools view -H ./hap500/temp.chr$2.$1.$i.vcf -S $1.1000G.S1.txt > ./hap500/temp.chr$2.$1.$i.vcf.txt
rm ./hap500/temp.chr$2.$1.$i.vcf
cut -f 10- ./hap500/temp.chr$2.$1.$i.vcf.txt > ./hap500/temp.chr$2.$1.$i.txt
rm ./hap500/temp.chr$2.$1.$i.vcf.txt
perl -p -i -e "s/\|/\t/g" ./hap500/temp.chr$2.$1.$i.txt
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' ./hap500/temp.chr$2.$1.$i.txt > ./hap500/temp.chr$2.$1.$i.trans.txt
rm ./hap500/temp.chr$2.$1.$i.txt
sort -u ./hap500/temp.chr$2.$1.$i.trans.txt | wc -l >> chr$2.$1.$1.hap.txt
rm ./hap500/temp.chr$2.$1.$i.trans.txt
done


#!/bin/bash 
for i in `seq 63854 250 171114067`
do
j=$((i+500))
tabix -h -p vcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 4:$i-$j > ./hap500/temp.chr6.$1.$i.vcf
bcftools view -H ./hap500/temp.chr6.$1.$i.vcf -S $1.1000G.S1.txt > ./hap500/temp.chr6.$1.$i.vcf.txt
rm ./hap500/temp.chr6.$1.$i.vcf
cut -f 10- ./hap500/temp.chr6.$1.$i.vcf.txt > ./hap500/temp.chr6.$1.$i.txt
rm ./hap500/temp.chr6.$1.$i.vcf.txt
perl -p -i -e "s/\|/\t/g" ./hap500/temp.chr6.$i.txt
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' ./hap500/temp.chr6.$1.$i.txt > ./hap500/temp.chr6.$1.$i.trans.txt
rm ./hap500/temp.chr6.$1.$i.txt
sort -u ./hap500/temp.chr6.$1.$i.trans.txt | wc -l >> chr6.$z.hap.txt
rm ./hap500/temp.chr6.$1.$i.trans.txt
done


for i in  ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI

for i in  CEU CHB YRI
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo sh ./hap500.sh $i Y 48129895 2655180 >>$i.job
qsub $i.job
done





chr1	249250621	10177
chr2	243199373	10179
chr3	198022430	60069
chr4	191154276	10005
chr5	180915260	10043
chr6	171115067	63854
chr7	159138663	14808
chr8	155270560	11740
chr9	146364022	10163
chr10	141213431	60494
chr11	135534747	61395
chr12	135006516	60181
chr13	133851895	19020047
chr14	115169878	19000017
chr15	107349540	20000041
chr16	102531392	60086
chr17	90354753	52
chr18	81195210	10083
chr19	78077248	60842
chr20	63025520	60343
chr21	59373566	9411239
chr22	59128983	16050075
chrX	51304566	60020
chrY	48129895	2655180




for i in {1..22} X Y
do
zcat ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | head -n 1 | awk '{print $2}'
done


/gpfs/home/guosa/hpc/project/cpgSNP

vcf2tai
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> $i.job
qsub $i.job
done


tabix -p vcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 4:10005-10494


ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf


# 2019-02-09
cd /gpfs/home/guosa/hpc/rheumatology/RA/TFBS_GWAS_RA_SNP


for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip wgEncodeRegTfbsClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3"\t"$4}' wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClusteredV3.hg19.bed

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz
gunzip wgEncodeRegDnaseClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3}' wgEncodeRegDnaseClusteredV3.bed > wgEncodeRegDnaseClusteredV3.hg19.bed

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/TFBS-GWAS-SNP/GWAS-RA-792.hg19.bed

bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed > GWAS-RA-R2.6.tfbs.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.tfbs.hg19.bed -b wgEncodeRegDnaseClusteredV3.bed | sort -u > GWAS-RA-R2.6.tfbs.DNase.hg19.bed
bedtools intersect -wa -a GWAS-RA-R2.6.tfbs.DNase.hg19.bed -b ~/hpc/db/hg19/CpGI.hg19.bed > GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.bed
bedtools sort -i GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.bed > GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.bed
bedtools merge -d 2000 -i GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.bed > GWAS-RA-R2.6.tfbs.DNase.CpGI.129.hg19.sort.merge.hg19.bed



awk '{print "chr"$1"\t"$2-1"\t"$2+1"\t"$3"\t"$4"\t"$5}' gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.bed  > gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.hg19.bed
bedtools intersect -a gnomad.genomes.r2.1.sites.rec.CpGI.merge.vcf.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed > CpGI.TFBS.SNP.hg19.txt 
bedtools intersect -wa -a CpGI.TFBS.SNP.hg19.sort.uni.txt -b wgEncodeRegDnaseClusteredV3.bed | sort -u > CpGI.TFBS.DNase.SNP.hg19.sort.uni.bed
bedtools intersect -wa -a CpGI.TFBS.DNase.SNP.hg19.sort.uni.bed -b /gpfs/home/guosa/hpc/db/hg19/BUR.GRCH37.hg19.bed | sort -u > CpGI.TFBS.DNase.BUR.hg19.bed



cp /gpfs/home/guosa/hpc/db/hg19/CpGI.hg19.bed ./
perl -p -i -e "s/chr//" CpGI.hg19.bed



cp GWAS-RA-792.R2.6.rsSNP.hg19.bed > GWAS-RA-792.R2.6.rsSNP.hg19.bed.input
perl -p -i -e "s/chr//" GWAS-RA-792.R2.6.rsSNP.hg19.bed.input





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





wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz

setwd("GTEx_Analysis_v7_eQTL")
qval_threshold=0.05
eqtl<-c()
for(i in c("Whole_Blood","Whole_Blood","Liver","Small_Intestine_Terminal_Ileum","Stomach","Colon_Sigmoid","Lung","Spleen","Ovary")){
data<-subset(read.table(paste(i,".v7.egenes.txt",sep=""),head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(eqtl,as.character(data[,19]))
}
length(table(eqtl))
eqtl.snp<-names(table(eqtl))
write.table(eqtl.snp,file="eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
input<-read.table("../GWAS-RA-792.R2.6.rsSNP.hg19.bed",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.PRA.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)


qval_threshold=0.05
file=list.files(pattern="*.v7.egenes.txt")
qval_threshold=0.05
eqtl<-c()
for(i in file){
data<-subset(read.table(i,head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(eqtl,as.character(data[,19]))
}
length(table(eqtl))
eqtl.snp<-names(table(eqtl))
write.table(eqtl.snp,file="eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
input<-read.table("../GWAS-RA-792.R2.6.rsSNP.hg19.bed",head=F)
output<-input[input[,4] %in% eqtl.snp,]
dim(output)
write.table(output,file="../GWAS-RA-792.R2.6.rsSNP.FullRA.eQTL.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)






# 2019-02-08
cp GWAS-RA-792.R2.6.rsSNP.hg19.bed GWAS-RA-792.R2.6.rsSNP.input.hg19.bed
perl -p -i -e "s/chr//" GWAS-RA-792.R2.6.rsSNP.input.hg19.bed

panel="GWAS-RA-792.R2.6.rsSNP.input"
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
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed

# 2019-02-08
cd /gpfs/home/guosa/hpc/db/hg19/CpGI
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



cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/PTPN_PADI
grep PTPN ~/hpc/db/hg19/refGene.hg19.bed > PTPN_PADI.hg19.bed
grep PADI ~/hpc/db/hg19/refGene.hg19.bed >> PTPN_PADI.hg19.bed
perl -p -i -e "s/chr//" PTPN_PADI.hg19.bed
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
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R PTPN_PADI.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.PTPN_PADI.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.5KGWASAID.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.5KGWASAID.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.5KGWASAID.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.5KGWASAID.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.5KGWASAID.merge.vcf.bed

awk '{print $1"\t"$2-20000"\t"$3+20000"\t"$4}' GWAS-immnue-3325_SNP.hg19.bed > 5KGWASAID.hg19.bed
awk '{print $1"\t"$2-10000"\t"$3+10000"\t"$4}' GWAS-immnue-3325_SNP.hg19.bed > GWAS_Catalog_5K_3000_Autoimmnue.hg19.bed
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R InnateDB.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.InnateDB.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf.bed

wget http://210.46.85.180:8080/SNP_linc_tfbs/download/linc_tfbs_chipseq_snp.txt
wget http://210.46.85.180:8080/SNP_linc_tfbs/download/linc_tfbs_snp.txt
awk '{print $2"\t"$17-1"\t"$17+1"\t"$16}' linc_tfbs_snp.txt | grep rs |sort -u > linc_tfbs_snp.hg19.bed
awk '{print $2"\t"$19-1"\t"$19+1"\t"$18}' linc_tfbs_chipseq_snp.txt |grep rs | sort -u > linc_tfbs_chipseq_snp.hg19.bed


cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/PTPN_PADI
grep ~/hpc/db

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/ASA/GWAS-immnue-3325_SNP.hg19.bed
perl -p -i -e "s/chr//" GWAS-immnue-3325_SNP.hg19.bed
awk '{print $1"\t"$2-10000"\t"$3+10000"\t"$4}' GWAS-immnue-3325_SNP.hg19.bed > GWAS_Catalog_5K_3000_Autoimmnue.hg19.bed
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \#bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R InnateDB.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.InnateDB.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.InnateDB.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.InnateDB.merge.vcf.bed




####### 2019-02-06

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




####### 2019-02-05

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




wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

awk '{print $4"\t"$5}' miRNA_targets_gain_by_SNPs_in_seed_regions.txt | sort -u | grep rs > miRNA-SNPs.list.txt
awk '{print $4"\t"$5}' miRNA_targets_loss_by_SNPs_in_seed_regions.txt | sort -u | grep rs >> miRNA-SNPs.list.txt
awk '{print $5}' miRNA_targets_gain_by_SNPs_in_seed_regions.txt | sort -u | grep rs > miRNA-SNPs.snponly.list.txt
awk '{print $5}' miRNA_targets_loss_by_SNPs_in_seed_regions.txt | sort -u | grep rs >> miRNA-SNPs.snponly.list.txt


perl -F"\s" -lane  "{print @F[6]}" InnateDB_genes.txt | grep -v name > InnateDB_genes.symbol.txt
genesymbol<-read.table("InnateDB_genes.symbol.txt")
db<-read.table("~/hpc/db/hg19/refGene.hg19.bed",head=F)
rlt<-db[db[,5] %in% genesymbol[,1],]
write.table(rlt,file="InnateDB.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

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


# 02/04/2019
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
grep -v "#" gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | grep "rs" > gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf.bed

wget VIP.hg19.bed
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo tabix -p vcf gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \| INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R VIP.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.vip.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.vip.merge.vcf

# Update 3
for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \| INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R Epigene.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Oz -o  gnomad.exomes.r2.1.sites.chr$i.rec.Epi.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.Epi.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.Epi.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.Epi.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.Epi.merge.vcf




 
 

rlt<-db[db[,5] %in% ef[,1],]
write.table(rlt,file="Epigene.hg19.bed",sep="\t",col.names=F,row.names=F,quote=F)

bedtools intersect -wao -a 001A_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt -b 001A_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt

mkdir bed
for i in `ls *.txt`
do
grep -v Start $i | awk '$1"\t"$1"\t"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10}' > ./bed/$i.bed
done

bedtools intersect -wa -a 001E_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed -b 001E_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
wc -l 001A_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
wc -l 001A_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed

rm *.snp.*

mv S001B_Tumor1-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed  001B_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001B_Tumor2-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001B_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001D_Tumor1_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001D_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001D_Tumor2_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001D_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001L_Tumor1-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001L_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001L_Tumor2-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001L_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001N_Tumor1_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001N_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed
mv S001N_Tumor2_mem.merged.HighConf.snpEff_ann.hg19_multianno.hg19_multianno.txt.bed 001N_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed

touch -m *
# D and N do not have normal control. 
mkdir SDN
mv 001D* SDN
mv 001N* SDN

cd /home/local/MFLDCLIN/guosa/hpc/project/LungBrainMetastasis/bed
for i in A B C E F G H I J K L M O
do
wc -l 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed > $i.txt
wc -l 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed >> $i.txt
bedtools intersect -wa -a 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed -b 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed | wc -l >> $i.txt
done

use strict;
my @file=glob("*.txt");
foreach my $file(@file){
open F,$file;
print "$file";
while(<F>){
my ($num)=split/\s+/;
print "\t$num";
}
print "\n"
}

use strict;
my @file=glob("*.bed");
my %cat;
my %type;
foreach my $file(@file){
open F,$file;
my($sam)=split/\_/,$file;
while(<F>){
my @line=split/\s+/;
$cat{$sam}{$line[6]}++;
$type{$line[6]}=$line[6];
}
}

foreach my $sam(sort keys %cat){
print "$sam";
 foreach my $type(sort keys %type){
 print "\t$cat{$sam}{$type}";
 }
 print "\n";
}


cd /home/local/MFLDCLIN/guosa/hpc/project/LungBrainMetastasis/bed
for i in A B C E F G H I J K L M O
do
bedtools intersect -wa -a 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed -b 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed >$i.overlap.bed
done



dim(GSE125362)
dim(GSE32413)
dim(GSE76885)
dim(GSE76807)
dim(GSE68698)
dim(GSE59785)
dim(GSE45485)
dim(GSE19617)

GSE76809

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/FGF6/Iron.Related.Gene.txt
wget https://raw.githubusercontent.com/Shicheng-Guo/GEO/master/SSc/mRNA/GEO-SSc-Vs-Normal-GPL6480.txt

iron.gene<-read.table("Iron.Related.Gene.txt")
SSc<-read.table("GEO-SSc-Vs-Normal-GPL6480.txt",head=T,sep="\t")
iron.SSc<-SSc[match(iron.gene[,1],SSc[,7]),]
hkg<-read.table("Human_Housekeeping_Gene_List.txt")
hkg.SSc<-SSc[match(hkg[,1],SSc[,7]),]
library("Haplin")
png("qqplot.iron.png")
pQQ(na.omit(iron.SSc$P.Value), nlabs =nrow(na.omit(iron.SSc)), conf = 0.95)
dev.off()
png("qqplot.hkg.png")
pQQ(na.omit(hkg.SSc$P.Value), nlabs =nrow(na.omit(hkg.SSc)), conf = 0.95)
dev.off()



cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/miRNA-SNP

wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt

data<-read.table("miRNA_targets_loss_by_SNPs_in_seed_regions.txt",head=T,sep="\t")
data<-read.table("miRNA_targets_gain_by_SNPs_in_seed_regions.txt",head=T,sep="\t")



02/01/2019
#PBS -N 22
#PBS -l nodes=1:ppn=1
cd /gpfs/home/guosa/hpc/db/Gnomad/vcf
bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz
bcftools view -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.001 & INFO/AF_eas>0.001 & INFO/vep ~ "missense_variant"' -R VIP.hg19.bed gnomad.exomes.r2.1.sites.chr22.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vip.vcf.bgz
bcftools sort gnomad.exomes.r2.1.sites.chr22.rec.vip.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.vcf.bgz
bcftools norm -d all gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.rmdup.vcf.bgz
bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.rmdup.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr22.rec.vip.sort.rmdup.biallelic.vcf.gz

time(bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.sort.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz)
time(bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.sort.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz)

bcftools sort check.temp.17.vcf > check.17.vcf
bcftools norm  -m + check.17.vcf > check.17.recover.vcf
bcftools norm -d both check.17.recover.vcf > check.17.recover.2.vcf
bcftools view -m2 -M2 -v snps check.17.recover.2.vcf > check.17.trim.vcf
grep rs1799966 check.17.trim.vcf | less -S 

bcftools norm -d both --threads=32 ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O z  -o chr1.vcf.gz

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=finished.pdf Supplemental_29Jan2019.pdf   Supplementary_Table_1.pdf Supplementary_Table_2_.pdf Supplementary_Table_3.pdf Supplementary_Figure_1.pdf Supplementary_Figure_2.pdf Supplementary_Figure_3.pdf Supplementary_Figure_4_.pdf Supplementary_Figure_5.pdf Supplementary_Figure_6.pdf supplementary_figure_7_.pdf supplementary_figure_8.pdf supplementary_figure_9.pdf Supplementary_Figure_10.pdf Supplementary_Figure_11.pdf 

tabix -p vcf -h gnomad.exomes.r2.1.sites.chr4.vcf.bgz 4:89011417-89152475 > ABCB1.vcf
cd ~/hpc/db/Gnomad/vcf

for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done


for i in {1..22} X 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools norm -m \+ gnomad.exomes.r2.1.sites.chr$i.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps -f PASS -i \'INFO\/AF[0] \> 0.001 \| INFO/AF_eas\>0.001 \& INFO\/vep \~ \"missense_variant\"\' -R VIP.hg19.bed  gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.vip.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.vip.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.vip.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.vip.merge.vcf

 


perl -lane 'print $_ if $_ !~/#/' VIP*vcf | awk '{print $1,$2,$3,$4,$5}' > VIP.snp.hg19.2.txt
bcftools view -m2 -M2 -v snps -f PASS -R VIP.hg19.bed gnomad.exomes.r2.1.sites.chr17.vcf.bgz > VIP.chr17.vcf
bcftools view -f PASS -R VIP.hg19.bed gnomad.exomes.r2.1.sites.chr17.vcf.bgz > VIP.chr17.temp.vcf
bcftools view -m2 -M2 -v snps VIP.chr17.temp.vcf | grep rs1799966
bgzip X.vcf
tabix X.vcf.gz
bcftools query -f PASS -i 'INFO/AF[0] > 0.005 | INFO/AF_eas>0.005' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\t%INFO/vep\n' X.vcf.gz  | grep -E "(stop gained)" 
bcftools view -f PASS ABCB1.vcf -S VIP.bed 'INFO/AF[0] > 0.01'  > ABCB1.pass.vcf
bcftools view -f PASS -i 'INFO/AF[0] > 0.1 & INFO/AF_eas>0.01' ABCB1.vcf.gz | less -S 
bcftools query -f PASS -i 'INFO/AF[0] > 0.005 | INFO/AF_eas>0.005' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\t%INFO/vep\n' ABCB1.vcf.gz  | grep -E "(missense_variant|frameshift|stop_gained|stop_lost|inframe_deletion)"

bcftools query  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.1 & INFO/vep ~ "[missense_variant|frameshift|stop_gained|stop_lost|inframe_deletion]"' -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\t%INFO/vep\n' gnomad.exomes.r2.1.sites.chr13.vcf.bgz | less -S 


bcftools view  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.1 & INFO/AF_eas>0.01 & INFO/vep ~ "[missense_variant|frameshift|stop_gained|stop_lost|inframe_deletion]"' ABCB1.vcf.gz | wc -l
bcftools view  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.001 & INFO/AF_eas>0.0001 & INFO/vep ~ "missense_variant"' ABCB1.vcf.gz | wc -l
bcftools view  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.0001 & INFO/AF_eas>0.00001 & INFO/vep ~ "frameshift"' ABCB1.vcf.gz | wc -l
bcftools view  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.0001 & INFO/AF_eas>0.00001 & INFO/vep ~ "stop_gained"' ABCB1.vcf.gz | wc -l
bcftools view  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.0001 & INFO/AF_eas>0.00001 & INFO/vep ~ "stop_lost"' ABCB1.vcf.gz | wc -l
bcftools view  -m2 -M2 -v snps -f PASS -i 'INFO/AF[0] > 0.0001 & INFO/AF_eas>0.00001 & INFO/vep ~ "inframe_deletion"' ABCB1.vcf.gz | wc -l




vep=A|missense_variant|MODERATE|ABCB1|ENSG00000085563|Transcript|ENST00000265724|protein_coding|29/29||ENST00000265724.3:c.3835C>T|ENSP00000265724.3:p.Arg1279Cys|4253|3835|1279|R/C|Cgc/Tgc|rs137996914|1||



grep IL7 IL-7_RA.tsv > output.txt
grep CD4 IL-7_RA.tsv >> output.txt
awk '{print $1}' output.txt > output2.txt
awk '{print $2}' output.txt >> output2.txt
sort -u output2.txt > output3.txt

x1<-"/gpfs/home/guosa/hpc/rheumatology/RA/ASA/miRNA-SNP/miRNA-UTR-double-pair.txt"
x2<-"/gpfs/home/guosa/hpc/rheumatology/RA/InnateDB_genes.txt"
data1<-read.table(x1,sep="\t")
data2<-read.table(x2,sep="\t")
newdata<-data1[which(data1[,1]%in% data2[,1]),]
write.table(newdata,file="miRNA-UTR-double-pair-RA-immune-Gene-GWAS.txt",sep="\t",quote=F,col.names=F,row.names=T)


cat miRNA_targets_gain_by_SNPs_in_seed_regions.txt  miRNA_targets_loss_by_SNPs_in_seed_regions.txt >> miRNA.txt
cat miRNA_gain_by_SNPs_in_gene_3utr.txt miRNA_loss_by_SNPs_in_gene_3utr.txt >> UTR.txt
awk '{print $2,"\t",$4,"\t",$5,"\t","miRNA"}' miRNA.txt | sort -u > miRNA.rs.txt
awk '{print $2,"\t",$6,"\t",$4,"\t","UTR"}' UTR.txt | sort -u > UTR.rs.txt
cat miRNA.rs.txt UTR.rs.txt > miRNA.UTR.rs.txt

wget ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz
install.packages("ipcid")
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt


cd /gpfs/home/guosa/hpc/tcga/xiong/hnsc 

mkdir mh450
mkdir fpkm-uq
mkdir miRNA
mkdir mCNS
mkdir clinic
mkdir SNP

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\mh450
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\fpkm-uq
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\miRNA
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\mCNS
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\clinic
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\xiong\hnsc\SNP
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\meth450Pancancer
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-20.txt

Set-Location \\mcrfnas2\bigdata\Genetic\Projects\shg047\tcga\rnaseqPancancer
C:\Admin\gdc-client.exe download --manifest gdc_manifest.2019-01-23.txt



plink --bfile /gpfs/home/guosa/hpc/db/Hapmap/hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode --extract snplist.txt --make-bed PGRA
plink --bfile /gpfs/home/guosa/hpc/db/Hapmap/hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode --extract snplist.txt --freq --make-bed --out PGRA 

awk '{print $7}' relationships_w_pops_051208.txt | sort -u | grep -v pop
for i in ASW CEU CHB CHD GIH JPT LWK MEX MKK TSI YRI
do
plink --file /gpfs/home/guosa/hpc/db/Hapmap/hapmap3/hapmap3_r1_b36_fwd.$i.qc.poly.recode --extract snplist.txt --freq --make-bed --out hapmap3_r1_b36_fwd.$i.qc.poly.recode.pgx
done



perl -lane '{print $1 if $_=~/(rs\d+)/}' pubmed_result.nejm.txt
C. Todd Stewart Award for Clinical Excellence

mutect2 
mutect1
scalpel
vardit

You can find Tai-Hsien Ou Yang with google and he is working in Columbia University. He graduated from BSEng, National Taiwan University which is No 1 university in Taiwan.  


First, relatives defined based on identity-by-descent and samples with sex inconstancies were excluded through the PLINK options of “--genome” and “--check-sex”, respectively. Also, samples with low call rate (<99%) were thrown out by option “--mind”. Next, PLINK option “--geno 0.95 (call rate < 95%)” was used to remove variants with completely missing data or having a low call rate. Also, additional variants were excluded based on deviation from Hardy-Weinberg Equilibrium (HWE), a below-threshold P value (HWE P < 0.000001), and minor allele count (MAC) < 2. Duplicated variants were removed.


plink --bfile ../FinalRelease_QC_20140311_Team1_Marshfield --mind 0.05 --geno 0.95 --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO --freq --out plink1
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO --extract fSNP.bed --range --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO.fSNP
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO.fSNP --freq --out plink2



Monomorphic variants

wget -r -l 1 -nd -e robots=off --reject jpg,html https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16256/suppl/

library("GEOquery")
GSE16256 <- getGEO("GSE16256")
data <- as.data.frame(exprs(GSE16256[[1]]))
phen <- pData(phenoData(GSE16256[[1]]))

phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  # status 1:control, 2:scz
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"schizophrenia"
phen2<-sapply(strsplit(as.character(phen$characteristics_ch1),"[:]"),function(x) (unlist(x)[2]))  # gender

data1=na.omit(data)
PCAPlot(t(data1),phen1,output="GSE41169.scz.normal.pdf",multifigure=T)  # status
PCAPlot(t(data1),phen2,output="GSE41169.gender.pdf",multifigure=T)  # gender

for i in {1..22}
do
perl -lane '{print $_ if @F[6]=~/rs/}' gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.annovar.txt > gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.annovar.hg19.bed
done


cd ~/gpfs/home/guosa/hpc/project/pmrp/phase1/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'R2\>0.8\' $i.dose.vcf.gz -Oz -o $i.dose.filter.vcf.gz >> $i.job
echo tabix -p vcf $i.dose.filter.vcf.gz >> $i.job
echo zcat $i.dose.filter.vcf.gz \| awk \'{print \$1,\$2,\$2,\$4,\$5}\' OFS=\"\\t\" \| grep -v \'#\' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl $i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ../annovar/$i -remove -protocol refGene,dbnsfp35c,gwasCatalog -operation gx,f,r -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub $i.job
done

table_annovar.pl chr12.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out chr12 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt

cd /gpfs/home/guosa/hpc/hemochromatosis/FGF6-imputation
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo sort -k1,1 -k2,2n $i.anno.txt \> $i.sort.anno.txt >> $i.job
echo bgzip -c $i.sort.anno.txt \> $i.anno.gz  >> $i.job
echo tabix -s1 -b2 -e2 $i.anno.gz >> $i.job
qsub $i.job
done


echo "##fileformat=VCFv4.1" > head.hdr
echo "##INFO=<ID=TAG1,Number=1,Type=String,Description="Yet another header line">" >> head.hdr
echo "##INFO=<ID=TAG2,Number=1,Type=String,Description="Yet another header line">" >> head.hdr
echo "##INFO=<ID=TAG3,Number=1,Type=String,Description="Yet another header line">" >> head.hdr


for i in chr{11..12}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view $i.dose.filter.vcf.gz -R $i.anno.txt -Oz -o $i.2LOF.vcf.gz >> $i.job
echo bcftools annotate -x INFO,FILTER,FORMAT/DS,FORMAT/GP $i.2LOF.vcf.gz -Oz -o $i.2LOF.vcf.tmp.gz >> $i.job
echo bcftools annotate -a $i.anno.gz -h head.hdr -c CHROM,POS,REF,ALT,-,-,TAG3 $i.2LOF.vcf.tmp.gz -Oz -o $i.update.vcf.gz >> $i.job
echo rm $i.2LOF.vcf.tmp.gz >> $i.job
qsub $i.job
done



#PBS -N chr12
#PBS -l nodes=1:ppn=1
cd /gpfs/home/guosa/hpc/hemochromatosis/FGF6-imputation
bcftools view chr12.dose.filter.vcf.gz -R chr12.anno.txt -Oz -o chr12.2LOF.vcf.gz
bcftools annotate -x INFO,FILTER,^FORMAT/GT chr12.2LOF.vcf.gz -Oz -o chr12.2LOF.vcf.tmp.gz
bcftools annotate -a chr12.anno.gz -h head.hdr -c CHROM,POS,REF,ALT,-,-,TAG3 chr12.2LOF.vcf.tmp.gz -Oz -o chr12.update.vcf.gz
rm chr12.2LOF.vcf.tmp.gz
bcftools annotate -x INFO,FILTER,FORMAT/PL chr12.2LOF.vcf.gz -Oz -o chr12.2LOF.vcf.tmp.gz


Rscript --vanilla 2LOF.R  All_samples_Exome_QC.chr19.vcf.DR2L0.8.lof.vcf.anno.gz All_samples_Exome_QC.phen


cp ~/hpc/project/pmrp/Exom2/2LOF/chr22.update.vcf  ./
cp ~/hpc/project/pmrp/Exom2/2LOF/2LOF.R  ./


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF

for i in {1..22}
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp10_Obesity_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp1_RA_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp7_Iron_C1_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp7_Iron_C2_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp3_PA_rev2_SampleIDs.Michigen.txt >> $i.job 
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp5_Thyroid_C1_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp5_Thyroid_C2_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_SSc_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ANA_rev2_SampleIDs.Michigen.txt >> $i.job
echo  Rscript --vanilla 2LOF.R chr$i.update.vcf /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ENA_rev2_SampleIDs.Michigen.txt >> $i.job
qsub $i.job
done

chr12.dose.filter.vcf.gz

table_annovar.pl ~/hpc/db/hg19/cpgSNP.hg19.bed.avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out cpgSNP.hg19.bed.avinput -remove -protocol Adipose_Subcutaneous.eqtl,Adipose_Visceral_Omentum.eqtl,Adrenal_Gland.eqtl,Artery_Aorta.eqtl,Artery_Coronary.eqtl,Artery_Tibial.eqtl,Brain_Amygdala.eqtl,Brain_Anterior_cingulate_cortex_BA24.eqtl,Brain_Caudate_basal_ganglia.eqtl,Brain_Cerebellar_Hemisphere.eqtl,Brain_Cerebellum.eqtl,Brain_Cortex.eqtl,Brain_Frontal_Cortex_BA9.eqtl,Brain_Hippocampus.eqtl,Brain_Hypothalamus.eqtl,Brain_Nucleus_accumbens_basal_ganglia.eqtl,Brain_Putamen_basal_ganglia.eqtl,Brain_Spinal_cord_cervical_c-1.eqtl,Brain_Substantia_nigra.eqtl,Breast_Mammary_Tissue.eqtl,Cells_EBV-transformed_lymphocytes.eqtl,Cells_Transformed_fibroblasts.eqtl,Colon_Sigmoid.eqtl,Colon_Transverse.eqtl,Esophagus_Gastroesophageal_Junction.eqtl,Esophagus_Mucosa.eqtl,Esophagus_Muscularis.eqtl,Heart_Atrial_Appendage.eqtl,Heart_Left_Ventricle.eqtl,Liver.eqtl,Lung.eqtl,Minor_Salivary_Gland.eqtl,Muscle_Skeletal.eqtl,Nerve_Tibial.eqtl,Ovary.eqtl,Pancreas.eqtl,Pituitary.eqtl,Prostate.eqtl,Skin_Not_Sun_Exposed_Suprapubic.eqtl,Skin_Sun_Exposed_Lower_leg.eqtl,Small_Intestine_Terminal_Ileum.eqtl,Spleen.eqtl,Stomach.eqtl,Testis.eqtl,Thyroid.eqtl,Uterus.eqtl,Vagina.eqtl,Whole_Blood.eqtl,bed -operation f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r --bedfile hg19_wgEncodeRegTfbsClusteredWithCellsV3.bed -argument '-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4' -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 


table_annovar.pl avinput ~/hpc/tools/annovar/humandb/  -buildver hg19 --csvout -out cpgSNP.hg19.bed.avinput -remove -protocol



bcftools view -i 'R2>0.9' chr12.dose.vcf.gz -Oz -o chr12.dose.filter.vcf.gz 


tabix -h chr12.dose.filter.vcf.gz  12:4543445-4554651 > chr12.dose.filter.FGF6.vcf

tabix -h FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz  12:4543445-4554651 > FinalRelease_QC_20140311_Team1_Marshfield.FGF6.vcf

zcat chr12.dose.filter.vcf.gz | awk '{print $1,$2,$3,$4,$5}' OFS="\t" | grep -v '#' > chr12.vcf.avinput


cat chr12.dose.filter.FGF6.vcf | awk '{print $1,$2,$2,$4,$5}' OFS="\t" | grep -v '#' > chr12.vcf.FGF6.avinput
table_annovar.pl chr12.vcf.FGF6.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out chr12 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt

tabix -h chr12.dose.vcf.gz  12:4543445-4554652 > chr12.dose.FGF6.vcf
cat chr12.dose.FGF6.vcf | awk '{print $1,$2,$2,$4,$5}' OFS="\t" | grep -v '#' > chr12.vcf.FGF6.avinput
table_annovar.pl chr12.vcf.FGF6.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out chr12 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt





~/hpc/hemochromatosis/plink
two_alof_all_combed_v4.phen
plink --file FGF6-C11 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out  FGF6-C11
--allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov --linear --out test


cd /gpfs/home/guosa/hpc/hemochromatosis/FGF6-imputation
for i in chr{11..12} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'R2\>0.6\' $i.dose.vcf.gz -Oz -o $i.dose.filter.vcf.gz >> $i.job
echo tabix -p vcf $i.dose.filter.vcf.gz >> $i.job
echo zcat $i.dose.filter.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl ./annovar/$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ./annovar/$i -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub  $i.job
done


for i in {6..22}
do
7za x chr_$i.zip -ppydG2EucK5bIhG &
done

cd /gpfs/home/guosa/hpc/project/pmrp/phase1/imputation
for i in {6..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo 7za x chr_$i.zip -ppydG2EucK5bIhG >>chr$i.job
qsub chr$i.job
done



 
cd ~/hpc/tools/annovar
annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
# just for allele frequency
annotate_variation.pl -downdb -webfrom annovar exac03 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2_all humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar gnomad_exome humandb -buildver hg38 &
# whole-exome data
annotate_variation.pl -downdb -webfrom annovar 1000g2015aug humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar kaviar_20150923 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar hrcr1 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cg69 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar gnomad_genome humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar dbnsfp30a humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg38 &
annotate_variation.pl -downdb esp6500siv2 humandb -buildver hg38 &
#  whole-genome data
annotate_variation.pl -downdb -webfrom annovar gerp++ humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cadd humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cadd13 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar fathmm humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar eigen humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar gwava humandb -buildver hg38  &
# for CNV
annotate_variation.pl -downdb -webfrom annovar dbscsnv11 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar spidex humandb -buildver hg38  &
# disease-specific variants
annotate_variation.pl -downdb -webfrom annovar clinvar_20160302 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cosmic70 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar icgc21 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar nci60 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar dbnsfp35c humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar dann humandb --buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar dann humandb --buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_all humandb --buildver hg19&

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gwasCatalog humandb/ &
annotate_variation.pl -buildver hg38 -downdb  tfbsConsSites humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar wgRna humandb/ &
annotate_variation.pl -buildver hg38 -downdb targetScanS humandb/ 

annotate_variation.pl -downdb -webfrom annovar dann humandb --buildver hg38 &

annotate_variation.pl -buildver hg38 -downdb  tfbsConsSites humandb/ 
annotate_variation.pl -buildver hg38 -downdb  tfbsConsSites humandb/ 


annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom ucsc gnomad_genome humandb/ 




for i in {1..22} X Y
do
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done

cd /gpfs/home/guosa/hpc/db/hg38/1000genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo vcftools --gzvcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --chr $i --recode --out ./annovar/chr$i >>chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo chr$i.recode.vcf  \> ../annovar/chr$i.vcf.avinput >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg38 --csvout -out ../annovar/chr$i.hg38 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done
rm *.e*
rm *.o*
findfncVar.pl

retrieve_seq_from_fasta.pl -format refGene -seqdir hg38_seq/ -outfile example.fa hg38_refGene.txt

cd /gpfs/home/guosa/hpc/tools/annovar/humandb/hg38_seq
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
gunzip http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz

cd ~/hpc/db/hg38/1000genome/annovar
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg38 --csvout -out ../annovar/chr$i.hg38 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done
rm *.e*
rm *.o*


echo convert2annovar.pl -format vcf4 chr$i.dose.filter.vcf.gz  \> ../annovar/chr$i.dose.vcf.avinput >> chr$i.job
/home/local/MFLDCLIN/guosa/hpc/db/hg38/1000genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

22      16053658        16053659        rs915675        A       C       0.857628
22      16055941        16055942        rs738830        C       T       0.819289
22      16058069        16058070        rs2843238       A       G       0.995008
22      16058882        16058883        rs2844897       A       G       0.649561
22      16060512        16060513        rs202125193     T       C       0.523562

http://faculty.washington.edu/tathornt/SISG2015/exercises/Lecture_3_PLINK_Code.txt

plink --bfile ~/hpc/db/hg19/1000Genome/chr5 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --extract two.test --recode fastphase --out two 


pdf("MRCI.PMRP.MAF.distribution.pdf")
hist(log(newdata$MAF,10),xlim=c(-5,0),main="Histogram of Log(MAF, 10)",xlab="Log(MAF, 10)",col="blue")
dev.off()

pdf("MRCI.PMRP.Hardy.P.distribution.pdf")
hist(log(data$P,10),main="Histogram of Hardy(HEW, 10)",xlab="Log(HWE, 10)",col="blue")
dev.off()

library("SKAT")
File.Bed=""
File.Bim=""
File.Fam=""
File.SetID=""
File.SSD=""
File.Info=""
SKAT.input<-Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)

--mpheno $phenotypeColumnNumber

rvtest --inVcf input.vcf --pheno phenotype.ped --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac

rvtest --inVcf input.vcf --pheno phenotype.ped --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac


rvtest --noweb --inVcf FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --pheno pmrp.exom1.2ALOF.phen --mpheno 2 --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac
rvtest --noweb --inVcf FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --pheno pmrp.exom1.2ALOF.phen --gene FGF6 --pheno-name PheTyp7_Iron_C1 --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac


data<-read.csv("two_alof_all_combed_v4.csv")
colnames(data)[1:4]<-c("FID","IID","LID","ESID")
write.table(data,file="pmrp.exom1.2ALOF.phen",sep="\t",quote=F,col.names=T,row.names=F)



sum(newdata$MAF>0.01)/nrow(newdata)

cd /gpfs/home/guosa/hpc/project/pmrp/phase1
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --maf 0.0001 --geno --freq counts --allow-no-sex --out ./FGF6/FGF6
cd /gpfs/home/guosa/hpc/project/pmrp/phase1/FGF6
wc -l FGF6.bim
plink --bfile FGF6 --recodeAD --out FGF6
less FGF6.raw


library(SKAT)
data(SKAT.example)
names(SKAT.example)
attach(SKAT.example)

bgzip -c FGF6.vcf > FGF6.vcf.gz
tabix -p vcf FGF6.vcf.gz


ped<-read.table()
fam<-read.table("FGF6.fam")
sam<-read.csv("../two_alof_all_combed_v4.csv")

remove=which(sam[match(fam[,1],sam[,2]),]$PheTyp7_Iron_C1=="-9")


obj<-SKAT_Null_Model(y.b ~ X, out_type="D")
SKAT(Z, obj)$p.value
SKAT(Z, obj, method="SKATO")$p.value
SKAT_CommonRare(Z, obj)$p.value
SKAT_CommonRare(Z, obj, method="A")$p.value


PTPN11: rs3750050, rs9640663
SNP1="rs3750050"
SNP2="rs9640663"
Genesymbol="PTPN11"
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 7:67256713-87256713 > $Genesymbol.vcf
plink --vcf $Genesymbol.vcf --maf 0.01 --make-bed --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out $Genesymbol
plink --bfile $Genesymbol --ld $SNP1 $SNP2
plink --bfile $Genesymbol --ld $SNP1 $SNP2 --freq
grep $SNP1 plink.frq
grep $SNP2 plink.frq


PTPN11: rs3750050, rs9640663
SNP1="rs3750050"
SNP2="rs9640663"
Genesymbol="PTPN11"
grep ~/hpc/db/hg19/PTPN11
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 7:67256713-87256713 > $Genesymbol.vcf
plink --vcf $Genesymbol.vcf --maf 0.01 --make-bed --keep ~/hpc/db/hg19/1000Genome/CHS.txt --freq --allow-no-sex --out $Genesymbol
grep $SNP1 $Genesymbol.frq
grep $SNP2 $Genesymbol.frq
plink --bfile $Genesymbol --ld $SNP1 $SNP2
plink --bfile $Genesymbol --ld $SNP1 $SNP2 
plink --vcf $Genesymbol.vcf --maf 0.01 --make-bed --keep ~/hpc/db/hg19/1000Genome/CHB.txt --freq --out $Genesymbol
grep $SNP1 $Genesymbol.frq
grep $SNP2 $Genesymbol.frq

plink --bfile $Genesymbol --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out $Genesymbol
plink --bfile $Genesymbol --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --recode HV --out $Genesymbol
plink --bfile delete --ld $SNP1 $SNP2
plink --bfile delete --r2 --ld-snp-list list.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ld_results
plink --bfile delete --r2 --ld-snp-list list.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --recode HV --out haploview

7:67256713-87256713

for i in {2011..2018}
do
cp $i.txt CCP$i.txt
done

cp GHCC.RA.CCP.uni GHCC.RA.CCP.uni.R
perl -p -i -e 's/\>//g' GHCC.RA.CCP.uni.R
perl -p -i -e 's/\<//g' GHCC.RA.CCP.uni.R
cd /home/local/MFLDCLIN/guosa/hpc/rheumatology/biobank/clinicaldata

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/biobank/clinicaldata")
library(Rtsne)  
data<-read.table("GHCC.RA.CCP.uni.R",head=T,sep="\t",row.names=1)
d <- dist(data) # euclidean distances between the rows
rlt<-Rtsne(d, is_distance=T)
pdf("tsne_RA_clinical.pdf")
d_tsne_1 = as.data.frame(rlt$Y)  
plot(d_tsne_1$Y,pch=16,cex=0.5)
dev.off()


input<-fit$points
input<-input[input[,1]> -5096,]
input<-input[input[,2]> -2596,]
input<-input[input[,2]< 2596,]
temp<-data[match(rownames(input),rownames(data)),]
pdf("MDSplot-chg4.pdf") 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=); 
dev.off()

library("tsne")
data<-read.table("GHCC.RA.CCP.uni.R",head=T,sep="\t",row.names=1)
tsne_iris = tsne(data, perplexity=50)
km <- kmeans(data,5,10000)
pc<-prcomp(data)

data<-read.table("len.tim.full.txt")
pdf("liveid-2.pdf")
par(mfrow=c(2,4))
barplot(table(data$V4),col="red",main="Year-Full")
pie(table(data$V4),main="Year-Full")
barplot(table(data$V5),col="red",main="Record-Full")
boxplot(data$V5,outline=F,main="Record-Full")
data=subset(data,V6=="live")
barplot(table(data$V4),col="red",main="Year-Live")
pie(table(data$V4),main="Year-Full")
barplot(table(data$V5),col="red",main="Record-Live")
boxplot(data$V5,outline=F,main="Record-Live")
dev.off()



cd ~/hpc/project/pmrp/Exom2/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
#echo bcftools view -i \'R2\>0.9\' $i.dose.vcf.gz -Oz -o $i.dose.filter.9.vcf.gz >> $i.job
#echo tabix -p vcf $i.dose.trim.9.vcf.gz >> $i.job
#echo zcat $i.dose.trim.9.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl ../annovar/$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ../annovar/$i -remove -protocol refGene,dbnsfp35c,gwasCatalog -operation gx,f,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub $i.job
done


# add INFO annotation to LOSS ALLLELS
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla readanno.R $i \2>> $i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo sort -k1,1 -k2,2n $i.mPred.2.anno.txt \> $i.sort.mPred.2.anno.txt >> $i.job
echo bgzip -c $i.sort.mPred.2.anno.txt \> $i.mPred.2.anno.gz  >> $i.job
echo tabix -s1 -b2 -e2 $i.mPred.2.anno.gz >> $i.job
qsub $i.job
done


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF/
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view ../imputation/$i.dose.trim.9.vcf.gz -R ../annovar/$i.sort.mPred.2.anno.txt -Oz -o ../2LOF/$i.vcf.gz >> $i.job
echo bcftools annotate -a ../annovar/$i.mPred.2.anno.gz -h ../annovar/head.hdr -c CHROM,POS,REF,ALT,-,-,TAG1 $i.vcf.gz -Ov -o $i.update.vcf >> $i.job
qsub $i.job
done




/gpfs/home/guosa/hpc/project/pmrp/phase1/Michigan/HLA
bcftools view -i 'R2>0.9' chr6.dose.vcf.gz  -Oz -o chr6.dose.R2.0.9.vcf.gz


for i in {1..22} X Y  
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo "zcat chr3.dose.vcf.gz | perl -lane '{next if /^#/;/R2=(\d+.\d+)/;print \$1}' > chr$i.R2.txt" >> $i.job
qsub $i.job
done


echo zcat chr3.dose.vcf.gz | perl -lane '{next if /^#/; /R2=(\d+.\d+)/;print $1}'


cd /gpfs/home/guosa/hpc/project/pmrp/phase1/Michigan
for i in {1..22} X Y  
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo vcftools --gzvcf FinalRelease_phase1.vcf.gz --chr $i --recode --out FinalRelease_phase1.chr$i.vcf >> $i.job
qsub $i.job
done


HLA-B chr6:31323299-31324,734
MICA chr6:31367561-31383090
NOTCH4 chr6:32163725-32165371

convert_bim_allele.pl: convert SNP allele coding between Illumina A/B alleles, Illumina 1/2 alleles, Illumina TOP strand alleles and dbSNP forward strand alleles.
/gpfs/home/guosa/hpc/tools/GenGen-1.0.1
convert_bim_allele.pl old.bim hh550_610.snptable -outfile new.bim -intype ilmnab -outtype top
convert_bim_allele.pl old.bim hh550_610.snptable -outfile new.bim -intype top -outtype dbsnp


tabix FinalRelease_phase1.vcf.gz 6:31323299-32165371


1       exm9324 0       6504575 0       G
1:6504575


echo rs2bed.pl ASD.hg19.avinput | qsub -N 'ASD'
cd /gpfs/home/guosa/hpc/rheumatology/SSc/CCHCR1/1000G
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150
plink --vcf ~/hpc/db/hg19/1000Genome/chr6.recode.vcf --snps rs141380757,rs72856720,rs11540822,rs72856718,rs130068,rs144885162 --make-bed --out cchcr1
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink --snps rs141380757,rs72856720,rs11540822,rs72856718,rs130068,rs144885162 --make-bed --out cchcr1


plink --vcf ~/hpc/db/hg19/1000Genome/chr6.uni.vcf --snps rs141380757,rs72856720,rs11540822,rs72856718,rs130068,rs144885162 --make-bed --out cchcr1

cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF
for i in 13 19 
do
for j in `ls /gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phe*.Michigen.txt`
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R chr$i.update.vcf $j >> $i.job
done
done


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation
for i in {1..22}
do
echo \#PBS -N chr$i.Guo  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo zcat chr$i.dose.filter.9.vcf.gz \| awk \'{print \$1,\$2,\".\",\$4,\$5}\' OFS=\"\\t\" \> chr$i.vepinput  >> $i.job
qsub $i.job
done
cat *vepinput | grep -v '#' > pmrp.exom2.vepinput.hg19.vcf




d1<-read.table("ASD.hg19.txt")
d2<-read.table("RsMergeArch.bcp",sep="\t")
old<-paste("rs",d2[,1],sep="")
new<-paste("rs",d2[,2],sep="")
xx<-data.frame(d1[which(as.character(d1[,1])%in%old),])
yy<-data.frame(d2[na.omit(match(as.character(d1[,1]),old)),])
zz<-data.frame(xx,yy)
zz$V2=paste("rs",zz$V2,sep="")
out<-data.frame(old=zz[,1],new=zz[,3])
write.table(out,file="ASD.old2new.SNP.txt",sep="\t",quote=F,col.names=T,row.names=F)

d1<-read.table("Raw_ASD_List.txt",sep="\t")
d2<-read.table("ShichengList.txt",sep="\t")

match(d2[,6],d1[,1])

write.table(out,file="ASD.old2new.SNP.txt",sep="\t",quote=F,col.names=T,row.names=F)




table_annovar.pl ASD.hg19.avinput ~/hpc/tools/annovar/humandb/ --thread 24 -buildver hg19 --csvout -out ASD.hg19.avinput-4 -remove -protocol refGene -operation g -nastring . -otherinfo  -polish  -polish 

convert2annovar.pl -format rsid  ASD.hg19.avinput -dbsnpfile ~/hpc/tools/annovar/humandb/ -out ASD.hg19.avinput-4

table_annovar.pl ASD.hg19.avinput ~/hpc/tools/annovar/humandb/ --thread 24 -buildver hg19 --csvout -out ASD.hg19.avinput-3 -remove -protocol gnomad_genome -operation f -nastring . -otherinfo  -polish  -polish 


table_annovar.pl ASD.hg19.avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out ASD.hg19.avinput -remove -protocol Adipose_Subcutaneous.eqtl,Adipose_Visceral_Omentum.eqtl,Adrenal_Gland.eqtl,Artery_Aorta.eqtl,Artery_Coronary.eqtl,Artery_Tibial.eqtl,Brain_Amygdala.eqtl,Brain_Anterior_cingulate_cortex_BA24.eqtl,Brain_Caudate_basal_ganglia.eqtl,Brain_Cerebellar_Hemisphere.eqtl,Brain_Cerebellum.eqtl,Brain_Cortex.eqtl,Brain_Frontal_Cortex_BA9.eqtl,Brain_Hippocampus.eqtl,Brain_Hypothalamus.eqtl,Brain_Nucleus_accumbens_basal_ganglia.eqtl,Brain_Putamen_basal_ganglia.eqtl,Brain_Spinal_cord_cervical_c-1.eqtl,Brain_Substantia_nigra.eqtl,Breast_Mammary_Tissue.eqtl,Cells_EBV-transformed_lymphocytes.eqtl,Cells_Transformed_fibroblasts.eqtl,Colon_Sigmoid.eqtl,Colon_Transverse.eqtl,Esophagus_Gastroesophageal_Junction.eqtl,Esophagus_Mucosa.eqtl,Esophagus_Muscularis.eqtl,Heart_Atrial_Appendage.eqtl,Heart_Left_Ventricle.eqtl,Liver.eqtl,Lung.eqtl,Minor_Salivary_Gland.eqtl,Muscle_Skeletal.eqtl,Nerve_Tibial.eqtl,Ovary.eqtl,Pancreas.eqtl,Pituitary.eqtl,Prostate.eqtl,Skin_Not_Sun_Exposed_Suprapubic.eqtl,Skin_Sun_Exposed_Lower_leg.eqtl,Small_Intestine_Terminal_Ileum.eqtl,Spleen.eqtl,Stomach.eqtl,Testis.eqtl,Thyroid.eqtl,Uterus.eqtl,Vagina.eqtl,Whole_Blood.eqtl,bed -operation f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r --bedfile hg19_wgEncodeRegTfbsClusteredWithCellsV3.bed -argument '-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4' -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 






time(grep vep gnomad.genomes.r2.1.sites.chr21.vcf.bgz | wc -l)    3483001

wc -l gnomad.genomes.r2.1.sites.chr21.vcf.bgz             
wc -l gnomad.genomes.r2.1.sites.chr21.vcf.bgz.annovar.txt        3206099


#!/usr/bin/perl
use strict;
use Cwd;
chdir "/gpfs/home/guosa/hpc/db/Gnomad";
open (F, "<", "gnomad.genomes.r2.1.sites.chr5.vcf.bgz" or die $!;
while(<F>){
next /^#/;
if(/vep=T|(\w+)_variant/){
print "$1\n";
}
}

for i in {4..22} X Y
do
qsub gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.job
done

# build annotation file for vcf as the input of bcftools
cd /gpfs/home/guosa/hpc/db/Gnomad
for i in `ls *gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf $i >> $i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/db/Gnomad
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo perl recode.pl gnomad.genomes.r2.1.sites.chr$i.vcf.bgz \> gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.annovar.txt >> $i.job
qsub $i.job
done


# build annotation file for vcf as the input of bcftools
cd /gpfs/home/guosa/hpc/db/Gnomad
for i in `ls *gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo gunzip $i >> $i.job
qsub $i.job
done


wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chrX.vcf.bgz 
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chrY.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chrX.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chrY.vcf.bgz.tbi



perl split.pl > cpgSNP.hg19.bed.V12.avinput
table_annovar.pl hg19_cpgSNP_eQTL.V6.txt ~/hpc/tools/annovar/humandb/ --thread 24 -buildver hg19 --csvout -out cpgSNP.hg19.bed.V12.avinput -remove -protocol knownGene,gwasCatalog,tfbsConsSites -
operation gx,r,r -nastring . -otherinfo  -polish  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 

table_annovar.pl allSNP150.hg19.bed.avinput ~/hpc/tools/annovar/humandb/ --thread 24 -buildver hg19 --csvout -out allSNP150.hg19.bed.avinput -remove -protocol refGene,gnomad_genome,dbnsfp35c,gwasCatalog,tfbsConsSites -operation gx,r,r -nastring . -otherinfo  -polish  -polish 

avinput="test.avinput"
table_annovar.pl $avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out $avinput -remove -protocol refGene -operation g -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 

avinput="/gpfs/home/guosa/hpc/db/hg19/snp150annovar"
table_annovar.pl $avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out $avinput -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 

1:15764961-15773084
for i in {1..22}
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done

 
plink 1.9 and 2.0 do not currently support haplotype association, sorry.  You may want to look at BEAGLE 3.3 for this (http://faculty.washington.edu/browning/beagle/b3.html )

# add INFO annotation to LOSS ALLLELS
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF
for i in `ls *.update.vcf`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla 2LOF.R $i >> $i.job
qsub $i.job
done


perl -lane '{next if /^#/;print "@F[0]\t@F[1]\t"@F[1]+1"\t@F[3]\t@F[4]"}' chr22.chip.vcf > chr22.chip.map &
perl -lane '{next if /^#/;print "@F[0]\t@F[1]\t"@F[1]+1"\t@F[3]\t@F[4]"}' chr22.imputation.vcf > chr22.imputation.map &

awk '{print $1,$2,$2+1,$4,$5}' OFS="\t" chr22.chip.vcf > chr22.chip.map &
awk '{print $1,$2,$2+1,$4,$5}' OFS="\t" chr22.imputation.vcf > chr22.imputation.map&

plink --vcf chr22.imputation.vcf.gz --recode --make-bed  --out chr22.imputation
plink --vcf chr22.chip.vcf.gz --recode --make-bed  --out chr22.chip

Rheumatoid Arthritis                DRB1*0401                                              6.2             rs660895/ rs6910071/rs3817964                                                                                

rs27044,rs17482078,rs10050860,rs30187,rs2287987

#pipeline to merge two vcf files
plink --bfile chr22.chip --list-duplicate-vars 
awk '{print $4}' plink.dupvar | grep -v ID > plink.dupvar.id 
plink --bfile chr22.chip --exclude plink.dupvar.id --make-bed --out chr22.chip.rmdup
plink --bfile chr22.imputation --list-duplicate-vars 
awk '{print $4}' plink.dupvar | grep -v ID > plink.dupvar.id 
plink --bfile chr22.imputation --exclude plink.dupvar.id --make-bed --out chr22.imputation.rmdup
plink --bfile chr22.imputation.rmdup --bmerge chr22.chip.rmdup --make-bed --out merge
plink --bfile chr22.chip.rmdup --flip merge-merge.missnp --make-bed --out chr22.chip.rmdup.flip
plink --bfile chr22.imputation.rmdup --bmerge chr22.chip.rmdup.flip --make-bed --out merge
plink --bfile chr22.imputation.rmdup --exclude merge-merge.missnp --make-bed --out chr22.imputation.rmdup.rm3
plink --bfile chr22.chip.rmdup.flip --exclude merge-merge.missnp --make-bed --out chr22.chip.rmdup.flip.rm3
plink --bfile chr22.imputation.rmdup.rm3 --bmerge chr22.chip.rmdup.flip.rm3 --make-bed --out merge
plink --bfile merge  --genome --out merge.ibd



cd ~/hpc/db/hg19/1000Genome/

#HLA-DRB1*0401
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs660895 --snp rs6910071 --snp rs3817964 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out HLA-DRB1-0401 --hap-r2 --recode --geno-r2  --geno-chisq --freq  --counts --hapcount 
#HLA-B27
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs13202464 --snp rs4349859 --snp rs3819299 --snp rs116488202 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out HLA-B27 --hap-r2  --geno-r2  --geno-chisq --freq  --counts 
#ERAP1
vcftools --gzvcf ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snp rs27044 --snp rs17482078 --snp rs10050860 --snp rs30187 --snp rs2287987 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out ERAP1 --hap-r2  --geno-r2  --geno-chisq --freq  --counts 

vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --maf 0.05 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --hapcount Hap.bed


# /gpfs/home/guosa/hpc/project/pmrp/phase2/plink
cd /gpfs/home/guosa/hpc/project/pmrp/phase2/vcf
for i in {1..22}
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo plink --bfile ../S_Hebbring_Unr.Guo.Forward --recode vcf --chr $i --out ../vcf/Exom2_Forward.chr$i >> chr$i.job
echo plink --bfile ../S_Hebbring_Unr.Guo.Forward --make-bed --chr $i --out ../plink/Exom2_Forward.chr$i >> chr$i.job
qsub chr$i.job
done


cd /gpfs/home/guosa/hpc/db/hg19/fa
for i in {1..22}
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=16 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo perl cpgsnp.pl chr$i >> chr$i.job 
qsub chr$i.job
done

avinput="test.avinput"

table_annovar.pl $avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out $avinput -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 




grep rs35705950 ~/hpc/tools/annovar/humandb/*eqtl*

table_annovar.pl ~/hpc/db/hg19/cpgSNP.hg19.bed.avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out cpgSNP.hg19.bed.avinput -remove -protocol Adipose_Subcutaneous.eqtl,Adipose_Visceral_Omentum.eqtl,Adrenal_Gland.eqtl,Artery_Aorta.eqtl,Artery_Coronary.eqtl,Artery_Tibial.eqtl,Brain_Amygdala.eqtl,Brain_Anterior_cingulate_cortex_BA24.eqtl,Brain_Caudate_basal_ganglia.eqtl,Brain_Cerebellar_Hemisphere.eqtl,Brain_Cerebellum.eqtl,Brain_Cortex.eqtl,Brain_Frontal_Cortex_BA9.eqtl,Brain_Hippocampus.eqtl,Brain_Hypothalamus.eqtl,Brain_Nucleus_accumbens_basal_ganglia.eqtl,Brain_Putamen_basal_ganglia.eqtl,Brain_Spinal_cord_cervical_c-1.eqtl,Brain_Substantia_nigra.eqtl,Breast_Mammary_Tissue.eqtl,Cells_EBV-transformed_lymphocytes.eqtl,Cells_Transformed_fibroblasts.eqtl,Colon_Sigmoid.eqtl,Colon_Transverse.eqtl,Esophagus_Gastroesophageal_Junction.eqtl,Esophagus_Mucosa.eqtl,Esophagus_Muscularis.eqtl,Heart_Atrial_Appendage.eqtl,Heart_Left_Ventricle.eqtl,Liver.eqtl,Lung.eqtl,Minor_Salivary_Gland.eqtl,Muscle_Skeletal.eqtl,Nerve_Tibial.eqtl,Ovary.eqtl,Pancreas.eqtl,Pituitary.eqtl,Prostate.eqtl,Skin_Not_Sun_Exposed_Suprapubic.eqtl,Skin_Sun_Exposed_Lower_leg.eqtl,Small_Intestine_Terminal_Ileum.eqtl,Spleen.eqtl,Stomach.eqtl,Testis.eqtl,Thyroid.eqtl,Uterus.eqtl,Vagina.eqtl,Whole_Blood.eqtl,bed -operation f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r --bedfile hg19_wgEncodeRegTfbsClusteredWithCellsV3.bed -argument '-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4','-colsWanted 4' -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 


newdata<-data[,c(1:53,55)]
info=lapply(newdata$Otherinfo,function(x) strsplit(as.character(x),split="[\t]"))
status=unlist(lapply(info,function(x) unlist(x)[1]))
rs=unlist(lapply(info,function(x) unlist(x)[2]))
context=unlist(lapply(info,function(x) unlist(x)[3]))
newdata2<-data.frame(newdata[,1:53],status,rs,context)
write.table(newdata2,file="CpG-SNP_eQTL_Functional_analysis.txt",sep="\t",quote=F,col.names=NA,row.names=T)

data<-newdata2
any(data[])

as.numeric(data[,6])

table_annovar.pl ~/hpc/db/hg19/cpgSNP.hg19.bed.avinput ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out cpgSNP.hg19.bed.avinput -remove -protocol Adipose_Subcutaneous.eqtl,Adipose_Visceral_Omentum.eqtl,Adrenal_Gland.eqtl,Artery_Aorta.eqtl,Artery_Coronary.eqtl,Artery_Tibial.eqtl,Brain_Amygdala.eqtl,Brain_Anterior_cingulate_cortex_BA24.eqtl,Brain_Caudate_basal_ganglia.eqtl,Brain_Cerebellar_Hemisphere.eqtl,Brain_Cerebellum.eqtl,Brain_Cortex.eqtl,Brain_Frontal_Cortex_BA9.eqtl,Brain_Hippocampus.eqtl,Brain_Hypothalamus.eqtl,Brain_Nucleus_accumbens_basal_ganglia.eqtl,Brain_Putamen_basal_ganglia.eqtl,Brain_Spinal_cord_cervical_c-1.eqtl,Brain_Substantia_nigra.eqtl,Breast_Mammary_Tissue.eqtl,Cells_EBV-transformed_lymphocytes.eqtl,Cells_Transformed_fibroblasts.eqtl,Colon_Sigmoid.eqtl,Colon_Transverse.eqtl,Esophagus_Gastroesophageal_Junction.eqtl,Esophagus_Mucosa.eqtl,Esophagus_Muscularis.eqtl,Heart_Atrial_Appendage.eqtl,Heart_Left_Ventricle.eqtl,Liver.eqtl,Lung.eqtl,Minor_Salivary_Gland.eqtl,Muscle_Skeletal.eqtl,Nerve_Tibial.eqtl,Ovary.eqtl,Pancreas.eqtl,Pituitary.eqtl,Prostate.eqtl,Skin_Not_Sun_Exposed_Suprapubic.eqtl,Skin_Sun_Exposed_Lower_leg.eqtl,Small_Intestine_Terminal_Ileum.eqtl,Spleen.eqtl,Stomach.eqtl,Testis.eqtl,Thyroid.eqtl,Uterus.eqtl,Vagina.eqtl,Whole_Blood.eqtl -operation f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 


cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/

annotate_variation.pl -build hg19 -downdb gwasCatalog humandb/
annotate_variation.pl -downdb -build hg19 gtexEqtlTissueArteryAorta  humandb/
annotate_variation.pl -downdb -build hg19 gtexEqtlCluster  humandb/
table_annovar.pl hg19_cpgSNP_eQTL.V6.txt ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out hg19_cpgSNP_eQTL.V6.bed -remove -protocol gwasCatalog,tfbsConsSites -operation r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 


perl split.pl > cpgSNP.hg19.bed.V12.avinput
table_annovar.pl hg19_cpgSNP_eQTL.V6.txt ~/hpc/tools/annovar/humandb/ --thread 24 -buildver hg19 --csvout -out cpgSNP.hg19.bed.V12.avinput -remove -protocol knownGene,gwasCatalog,tfbsConsSites -operation gx,r,r -nastring . -otherinfo  -polish  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 
 



awk '{print $1,$2,$3,$4"_"$5"_"$6}' OFS="\t" wgEncodeRegTfbsClusteredWithCellsV3.bed > ~/hpc/tools/annovar/humandb/hg19_wgEncodeRegTfbsClusteredWithCellsV3.txt
perl -p -i -e 's/chr//g' ~/hpc/tools/annovar/humandb/hg19_wgEncodeRegTfbsClusteredWithCellsV3.bed
perl -p -i -e 's/,/_/g' ~/hpc/tools/annovar/humandb/hg19_wgEncodeRegTfbsClusteredWithCellsV3.bed
perl -p -i -e 's/\"//g' ~/hpc/db/hg19/cpgSNP.hg19.bed.avinput.hg19_multianno.csv
perl -p -i -e 's/\+//g' /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt
 
 

# build annotation file for vcf as the input of bcftools
cd /gpfs/home/guosa/hpc/project/pmrp/phase2/vcf
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo cd $(pwd) >> $i.job
echo bgzip -c Exom2_Forward.$i.vcf \> Exom2_Forward.$i.vcf.gz  >> $i.job
echo tabix -p vcf Exom2_Forward.$i.vcf.gz >> $i.job
qsub $i.job
done

# build plink for imputated dataset
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo plink --vcf $i.dose.filter.9.vcf.gz --make-bed --out ./Plink/Exom2.M.$i  >> $i.job
qsub $i.job
done

# plink merge
cd ~/hpc/project/pmrp/Exom2/IBD/

annotate_variation.pl -filter -dbtype generic -genericdbfile annovar_eqtl.hg19.bed -build hg19 -out ex1 test.anno ./
awk '{print $1,$2,$3,$4,$5}' OFS="\t" annovar_eqtl.hg19.sort.bed | head -n 200 | sort -u > test.anno

/gpfs/home/guosa/hpc/db/hg19/GTEx_Analysis_v7_eQTL
/gpfs/home/guosa/hpc/db/hg19/GTEx_Analysis_v7_eQTL



plink --bfile /gpfs/home/guosa/hpc/project/pmrp/phase2/plink/Exom2_Forward.chr21 --bmerge /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation/Plink/Exom2.M.chr21  --make-bed --out chr21
plink --bfile /gpfs/home/guosa/hpc/project/pmrp/phase2/plink/Exom2_Forward.chr21 --flip chr21-merge.missnp --make-bed --out Exom2_Forward.flip.chr21
plink --bfile /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation/Plink/Exom2.M.chr21 --flip chr21-merge.missnp --make-bed --out Exom2.M.chr21
plink --file ./Exom2_Forward.flip.chr21 --bmerge /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation/Plink/Exom2.M.chr21  --make-bed --out chr21


awk '{print $1,$4}' /gpfs/home/guosa/hpc/project/pmrp/phase2/plink/Exom2_Forward.chr21.bim > chr21.loci
bcftools view -T chr21.loci ..
view -T chr21.loci ../../Exom2/imputation/chr21.dose.vcf.gz

perl -p -i -e 's/|/_/g' chr13.9.hg19_multianno.csv

echo zcat chr13.dose.trim.9.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
 


java–Xmx 20g –jar beagle..jar 
cd /home/guosa/hpc/project/pmrp/Exom2/IBD/
bgzip -c /home/guosa/hpc/project/pmrp/phase2/phase/beagle/S_Hebbring_Unr.Guo.Forward.chr1.vcf > /home/guosa/hpc/project/pmrp/Exom2/IBD/F1.vcf.gz
tabix -p vcf F1.vcf.gz
1:752721
bcftools view -t 1:752721 F1.vcf.gz   -Ov -o X.vcf
bcftools view -t 1:752721 chr1.dose.vcf.gz  -Ov -o Y.vcf
bcftools isec F1.vcf.gz  chr1.dose.vcf.gz  -p bcfisec

cd ~/hpc/project/pmrp/Exom2/IBD
REF="~/hpc/tools/beagle/chr10.1kg.phase3.v5a.vcf.gz"
GTF1="~/hpc/project/pmrp/Exom2/imputation/chr10.dose.filter.9.vcf.gz"
GTF2="~/hpc/project/pmrp/phase2/vcf/Exom2_Forward.chr10.vcf.gz"

/gpfs/home/guosa/hpc/project/pmrp/phase2/vcf/Exom2_Forward.chr10.vcf.gz
/gpfs/home/guosa/hpc/project/pmrp/phase2/vcf/Exom2_Forward.chr10.vcf.gz

grep -v EUR ~/hpc/tools/beagle/integrated_call_samples_v3.20130502.ALL.panel| cut -f1 > non.eur.excl
java -jar ~/hpc/tools/beagle/conform-gt.24May16.cee.jar ref=/home/guosa/hpc/tools/beagle/chr10.1kg.phase3.v5a.vcf.gz gt=/home/guosa/hpc/project/pmrp/Exom2/imputation/chr10.dose.filter.9.vcf.gz chrom=10 out=chr10.M.vcf.gz excludesamples=non.eur.excl
java -jar ~/hpc/tools/beagle/conform-gt.24May16.cee.jar ref=/home/guosa/hpc/tools/beagle/chr10.1kg.phase3.v5a.vcf.gz gt=/home/guosa/hpc/project/pmrp/phase2/vcf/Exom2_Forward.chr10.vcf.gz chrom=10 out=chr10.B.vcf.gz excludesamples=non.eur.excl


plink --bfile /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Rel --snp rs1014988 --freq --out Rel
plink --bfile /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr --snp rs1014988 --freq --out Unr



# build annotation file for vcf as the input of bcftools
cd /gpfs/home/guosa/hpc/project/pmrp/phase2/phase/beagle
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bgzip -c S_Hebbring_Unr.Guo.Forward.$i.vcf \> $i.vcf.gz  >> $i.job
echo tabix -p $i.vcf.gz >> $i.job
qsub $i.job
done


10      1230968 10:1230968      C       T       .       .       TAG3=ADARB2     GT      1|0     1|1     1|1     0|1     0|1     1|1     1|1     0|0     0|1     0|1     1|1     1|1     0|0     0|0     0|1     0|1     1|1     1|1     1|0

ADARB2
rs2271275                 


tabix chr1.dose.vcf.gz 1:752721-752721 > X.vcf
tabix F1.vcf.gz 1:752721-752721 > Y.vcf
	
java–Xmx[GB]g –jar beagle.[version].jar [arguments]

data<-read.table("chr2.update.vcf",head=T,sep="\t",skip=16)

for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla sillyScript.R $i >> $i.job
qsub $i.job
done

Rscript --vanilla sillyScript.R iris.txt out.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("~/hpc/project/pmrp/Exom2/imputation")
file=paste(args[1],".dose.trim.9.vcf",sep="")
gdata<-read.table("chr22.update.vcf",head=T,sep="\t",skip=16,comment.char = "/",check.names = F)
save(gdata,file=paste(file,".RData",sep=""))

case<-
control<-
for(i in names(which(table(gdata$INFO)>=2))){
gdata[which(gdata$INFO %in% i),]
}


grep '^#' chr22.update.vcf > Exom2.vcf
for i in chr{1..22}
do
grep -v '^#' $i.update.vcf >> Exom2.vcf
done


cd ~/hpc/project/pmrp/Exom2/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view -i \'R2\>0.9\' $i.dose.vcf.gz -Oz -o $i.dose.filter.9.vcf.gz >> $i.job
echo tabix -p vcf $i.dose.trim.9.vcf.gz >> $i.job
echo zcat $i.dose.trim.9.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
echo table_annovar.pl ../annovar/$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 1 -buildver hg19 --csvout -out ../annovar/$i -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> $i.job
qsub $i.job
done

perl -p -i -e 's/|/_/g' chr13.9.hg19_multianno.csv

echo zcat chr13.dose.trim.9.vcf.gz \| awk \'{print \$1,\$2,\$3,\$4,\$5}\' OFS=\"\t\" \| grep -v '#' \> $i.vcf.avinput >> $i.job
 
 
# build annovar annotation for dbsnp150
table_annovar.pl hg19_avsnp150.txt ~/hpc/tools/annovar/humandb/ --thread 12 -buildver hg19 --csvout -out hg19_avsnp150 -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt 


cd ~/hpc/project/pmrp/Exom2/imputation
for i in chr{1..22} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf $i.dose.trim.9.vcf.gz >> $i.job
qsub $i.job
done

bcftools annotate -x FILTER,INFO,^FORMAT/GT,FORMAT/PL,FORMAT chr22.vcf.gz -Ov -o chr22.trim.vcf

# add INFO annotation to LOSS ALLLELS
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla readanno.R $i >> $i.job
qsub $i.job
done

open F,""
my $i=1;
while(<F>){
if($i==197836){
print ""
}
$i++
}


# 
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF/
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bcftools view ../imputation/$i.dose.trim.9.vcf.gz -R ../annovar/$i.anno.txt -Oz -o ../2LOF/$i.vcf.gz >> $i.job
echo bcftools annotate -x INFO,FILTER ../2LOF/$i.vcf.gz -Oz -o ../2LOF/$i.vcf.tmp.gz >> $i.job
echo bcftools annotate -a ../annovar/$i.anno.gz -h ../annovar/head.hdr -c CHROM,POS,REF,ALT,-,-,TAG3 $i.vcf.tmp.gz -Oz -o $i.update.vcf.gz >> $i.job
echo rm ../2LOF/$i.vcf.tmp.gz >> $i.job
qsub $i.job
done


# build annotation file for vcf as the input of bcftools
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo sort -k1,1 -k2,2n $i.anno.txt \> $i.sort.anno.txt >> $i.job
echo bgzip -c $i.sort.anno.txt \> $i.anno.gz  >> $i.job
echo tabix -s1 -b2 -e2 $i.anno.gz >> $i.job
qsub $i.job
done

# add INFO annotation to LOSS ALLLELS
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF/
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo bgzip -c ../imputation/$i.dose.trim.9.vcf  \> ../imputation/$i.dose.trim.9.vcf.gz >> $i.job
echo bcftools view ../imputation/$i.dose.trim.9.vcf.gz -R ../annovar/$i.anno.txt -Oz -o ../2LOF/$i.vcf.gz >> $i.job
echo bcftools annotate -a ../annovar/$i.anno.gz -h ../annovar/head.hdr -c CHROM,POS,REF,ALT,-,-,TAG3 $i.vcf.gz -Oz -o $i.update.vcf.gz >> $i.job
done
qsub $i.job



grep '^#' chr22.update.vcf > Exom2.vcf
for i in chr{1..22}
do
grep -v '^#' $i.update.vcf >> Exom2.vcf
done






vcftools --gzvcf chr22.dose.trim.9.vcf --bed ../annovar/chr22.
bcftools view chr22.dose.trim.9.vcf -R ../annovar/chr22.anno.snp
bcftools view chr22.dose.trim.9.vcf.gz -R ../annovar/chr22.anno.txt -Oz -o ../2LOF/chr22.vcf.gz

bcftools view chr22.dose.trim.9.vcf.gz -R ../annovar/chr22.anno.txt -Ov -o ../2LOF/chr22.vcf
bcftools annotate -x FILTER,INFO,^FORMAT/GT,FORMAT/PL,FORMAT chr22.vcf.gz -Ov -o chr22.trim.vcf

bcftools annotate -a src.bcf -c ID,QUAL,+TAG dst.bcf
bcftools annotate -a ../2LOF/chr22.vcf.gz -c ID,QUAL,+TAG dst.bcf


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/2LOF
bcftools annotate -a ../annovar/chr22.anno.gz -h ../annovar/head.hdr -c CHROM,POS,REF,ALT,TAG1,TAG2,TAG3 chr22.trim.vcf -Oz -o chr22.trim.update.vcf.gz
bcftools annotate -a src.bcf -c ID,QUAL,+TAG dst.bcf
bcftools annotate -a ../2LOF/chr22.vcf.gz -c ID,QUAL,+TAG dst.bcf
bcftools annotate -a ../annovar/chr22.anno.gz -c -,-,-,-,-,INFO,Tag chr22.trim.vcf
less -S ../annovar/chr21.anno.gz
less -S ../annovar/head.hdr
vim ../annovar/head.hdr


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("~/hpc/project/pmrp/Exom2/annovar")
file=paste(args[1],".9.hg19_multianno.csv",sep="")
data<-read.table(file,head=T,sep=",",check.names = F)

dataset<-subset(data,FATHMM_pred =="D"| MetaSVM_pred=="D"|MetaLR_pred =="D"| gwasCatalog !="." | targetScanS !="." |tfbsConsSites !="." | Eigen>0.2)
ID=paste(dataset$Chr,":",dataset$Start,sep="")
write.table(ID,file=paste(args[1],"anno.txt",sep="."),col.names=F,row.names=F)

data<-read.table("chr21.9.hg19_multianno.csv",head=T,sep=",",check.names = T)
data<-data.frame(data)



dataset<-subset(data,FATHMM_pred =="D"| MetaSVM_pred=="D"|MetaLR_pred =="D" | M.CAP_pred=="D"|gwasCatalog !="." | targetScanS !="." |tfbsConsSites !=".")
ID=paste(dataset$Chr,":",dataset$Start,sep="")
write.table(ID,file=paste(args[1],"anno.txt",sep="."),col.names=F,row.names=F)

GscPred<-apply(data[,grep("pred",colnames(data))],1,function(x) sum(x %in% "D")>=1)

xx<-data[,grep("pred",colnames(data))]

xx[GscPred,]

data<-read.table(annovar_csv,head=T,sep=",",check.names = T)
RiskAAC<-function(data){
risk1<-grep("D",data$SIFT_pred)
risk2<-grep("D|P",data$Polyphen2_HDIV_pred)
risk3<-grep("D|P",data$Polyphen2_HVAR_pred)
risk4<-grep("D",data$LRT_pred)
risk5<-grep("D",data$MutationTaster_pred)
risk6<-grep("H|M",data$MutationAssessor_pred)
risk7<-grep("D",data$FATHMM_pred)
risk8<-grep("D",data$PROVEAN_pred)
risk9<-grep("D",data$MetaSVM_pred)
risk10<-grep("D",data$MetaLR_pred)
risk11<-grep("D",data$fathmm.MKL_coding_pred)
risk12<-grep("D",data$M.CAP_pred)
risk13<-grep("Name",data$gwasCatalog)
risk14<-grep("Name",data$tfbsConsSites)
risk15<-grep("Name",data$targetScanS)
risk16<-grep("Name",data$wgRna)
rlt=unique(c(risk1,risk2,risk3,risk4,risk5,risk6,risk7,risk8,risk9,risk10,risk11,risk12,risk13,risk14,risk15,risk16))
return(rlt)
}


cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -webfrom annovar ljb23_fathmm humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar ljb23_metasvm humandb -buildver hg19  
annotate_variation.pl -downdb -webfrom annovar ljb23_metalr humandb -buildver hg19  
annotate_variation.pl -downdb -webfrom annovar eigen humandb -buildver hg19  

# PBS code with R2 filter
# Memory Size should be careful.
cd ~/hpc/project/pmrp/Exom2/imputation
for i in 7 20 22
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=16 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.dose.vcf.9.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../chr$i.9 -remove -protocol refGene,dbnsfp33a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done


ln -s /mnt/bigdata/Genetic/Projects/shg047 hpc
zcat chr22.dose.filter.vcf.gz | awk -F'[R2=\s]' '$5>0.1{print}' 
zcat chr22.dose.filter.vcf.gz | awk -F'[=\s]' '{print $4}'  | head 
/home/yez/Schrodi_2ALOF/script
data<-read.table("chr22.dose.trim.9.vcf",head=T,sep="\t",skip=16)

# PBS code with R2 filter
cd ~/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo bcftools annotate -x 'FORMAT' chr$i.dose.filter.9.vcf.gz -Oz -o chr$i.dose.trim.9.vcf.gz >> chr$i.job
qsub chr$i.job
done



cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/annovar
for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> $i.job
echo Rscript --vanilla readanno.R $i >> $i.job
qsub $i.job
done

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("~/hpc/project/pmrp/Exom2/annovar")
file=paste(args[1],".9.hg19_multianno.csv",sep="")
data<-read.table(file,head=T,sep=",",check.names = F)
dataset<-subset(data,FATHMM_pred =="D"| MetaSVM_pred=="D"  |MetaLR_pred =="D" | gwasCatalog !="." | targetScanS !="." |tfbsConsSites !="." | Eigen>0.2)
ID=paste(dataset$Chr,":",dataset$Start,sep="")
write.table(ID,file=paste(args[1],"anno.txt",sep="."),col.names=F,row.names=F)


data<-read.table("chr22.dose.trim.9.vcf",comment.char="",head=T,skip=16)
annotate_variation.pl -webfrom annovar -downdb avdblist
## Step 1:  Download Annotation Files
cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -build hg19 gtexEqtlTissueArteryAorta  humandb/
annotate_variation.pl -downdb -build hg19 gtexEqtlCluster  humandb/

-colsWanted all
-colsWanted 4,5,6
cd ~/hpc/project/pmrp/Exom2/annovar
table_annovar.pl chr22.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out chr22 -remove -protocol refGene,gwasCatalog,gtexEqtlCluster -operation gx,r,r -nastring . 


annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp33a  humandb/
annotate_variation.pl -downdb -webfrom annovar -build hg19 ljb23_all  humandb/
annotate_variation.pl -downdb -webfrom annovar ljb23_fathmm humandb -buildver hg19  &
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar eigen humandb/
annotate_variation.pl -downdb -webfrom annovar -downdb avdblist humandb/ -build hg19 
annotate_variation.pl -downdb -webfrom annovar -build hg19 seq humandb/hg19_seq/
annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19
retrieve_seq_from_fasta.pl -format genericGene -seqdir ./humandb/hg19_seq/ -outfile ./humandb/hg19_wgEncodeGencodeBasicV19Mrna.fa ./humandb/hg19_wgEncodeGencodeBasicV19.txt 
annotate_variation.pl -downdb -webfrom annovar ljb23_metalr humandb/ -build hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_gerp++ humandb/ -build hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_all humandb/ -build hg19 &
annotate_variation.pl -downdb -webfrom annovar ljb23_metasvm humandb/ -build hg19 &
annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19
annotate_variation.pl -buildver hg19 -downdb wgEncodeGencodeBasicV19 humandb/
annotate_variation.pl -buildver hg19 -downdb dgvMerged humandb/
annotate_variation.pl -buildver hg19 -downdb tfbsConsSites humandb/
annotate_variation.pl -buildver hg19 -downdb wgRna humandb/
annotate_variation.pl -buildver hg19 -downdb targetScanS humandb/
annotate_variation.pl -buildver hg19 -downdb gwasCatalog humandb/
annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/ &
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar regsnpintron humandb/ &

# just for allele frequency
cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -webfrom annovar exac03 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2_all humandb -buildver hg19 &

cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -webfrom annovar gnomad_exome humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar gnomad_genome humandb -buildver hg19 &

# whole-exome data
annotate_variation.pl -downdb -webfrom annovar 1000g2015aug humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar kaviar_20150923 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar hrcr1 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cg69 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar gnomad_genome humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar dbnsfp30a humandb -buildver hg19 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg19 &
annotate_variation.pl -downdb esp6500siv2 humandb -buildver hg19 &
#  whole-genome data
cd ~/hpc/tools/annovar
annotate_variation.pl -downdb -webfrom annovar gerp++ humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cadd humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cadd13 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar fathmm humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar eigen humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar gwava humandb -buildver hg19  &
# for CNV
annotate_variation.pl -downdb -webfrom annovar dbscsnv11 humandb -buildver hg19  
annotate_variation.pl -downdb -webfrom annovar spidex humandb -buildver hg19  
# disease-specific variants
annotate_variation.pl -downdb -webfrom annovar clinvar_20160302 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar cosmic70 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar icgc21 humandb -buildver hg19  &
annotate_variation.pl -downdb -webfrom annovar nci60 humandb -buildver hg19  &
# miSNP
# CpG-SNP
# eQTL

## Step 2:  Creat annovar input from VCF files
echo convert2annovar.pl -format vcf4 chr$i.dose.filter.vcf.gz  \> ../annovar/chr$i.dose.vcf.avinput >> chr$i.job

## Step 3:  Annotation input vcf 
cd ~/hpc/project/pmrp/Exom2/annovar
table_annovar.pl chr22.dose.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out chr22 -remove -protocol refGene,cytoBand,exac03,ljb23_all,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . 

awk -F"," '{print $6}' chr22.ts.hg19_multianno.csv | sort -u

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV28lift37.txt.gz
 
 table_annovar.pl ../annovar/chr22.dose.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr22 -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfb
 
 annotate_variation.pl -geneanno -dbtype refGene -buildver hg19 ../annovar/chr22.dose.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/
 
table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt

 
 
:miR

grep 17072086 S_Hebbring_Unr.Guo.Forward.chr22.vcf
grep 17073134 S_Hebbring_Unr.Guo.Forward.chr22.vcf
grep 17073275 S_Hebbring_Unr.Guo.Forward.chr22.vcf


exonic;splicing ncRNA_exonic;splicing, UTR5, UTR3, ncRNA_exonic, exonic, 



# PBS code with R2 filter
cd ~/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo bcftools view -i \'R2\>0.9\' chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.9.vcf.gz >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq chr$i.dose.filter.9.vcf.gz  \> ../annovar/chr$i.dose.vcf.9.avinput >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.dose.vcf.9.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr$i.9 -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done


convert2annovar.pl -format vcf4 example/ex2.vcf -outfile ex2.avinput -allsample -withfreq -include



# Annotation and Compound heterozygotes scanning
cd ~/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo -withzyg chr$i.dose.filter.9.vcf.gz  \> ../annovar/chr$i.dose.vcf.9.avinput >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.dose.vcf.9.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr$i.9 -remove -protocol refGene,ljb23_fathmm,ljb23_metasvm,ljb23_metalr,eigen,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done




 
#PBS -N chr21
#PBS -l nodes=1:ppn=1
cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation
bcftools view -i 'R2>0.9' chr21.dose.vcf.gz -Oz -o chr21.dose.filter.9.vcf.gz
convert2annovar.pl -format vcf4 -allsample -withfreq chr21.dose.filter.9.vcf.gz > ../annovar/chr21.dose.vcf.9.avinput
table_annovar.pl ../annovar/chr21.dose.vcf.9.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr21 -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt


cd ~/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.dose.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr$i -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
done

table_annovar.pl ../annovar/test.vcf /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr21 -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,r,f,f,r,r,r,r -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt


# PBS code without R2 filter
cd ~/hpc/project/pmrp/Exom2/imputation
for i in 7 20 22
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=16 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq chr$i.dose.vcf.gz  \> ../annovar/chr$i.dose.vcf.avinput >> chr$i.job
echo annotate_variation.pl -out ../annovar/chr$i.annovar.txt -build hg19 ../annovar/chr$i.dose.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.dose.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ../annovar/chr$i -remove -protocol refGene,cytoBand,exac03,dbnsfp30a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation g,r,f,f,r,r,r,r -nastring . >> chr$i.job
qsub chr$i.job
done


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation

1544_14014
11863_2780
2653_1120
10992_11507
11321_7649
 


/gpfs/home/guosa/hpc/db/ExAC

cp /mnt/bigdata/Center/CHG/valenzur/GWAS/ExomeII/rayner_goncalo_without_update_script/imputation_results/* ./


cd /mnt/bigdata/Center/CHG/valenzur/GWAS/ExomeII/rayner_goncalo_without_update_script/imputation_results/

/home/guosa/hpc/project/pmrp/Exom2

/home/guosa/hpc/tools/annovar

echo \#PBS -N chr$i  > phen$j.chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> phen$j.chr$i.job
echo cd $(pwd) >> phen$j.chr$i.job
echo plink --file exomechip_SNV_PASS_BEAGLE_chr$i\_phased_sel2 --allow-no-sex --pheno two_alof_all_combed_v4.phen --mpheno $j --assoc --out /gpfs/home/guosa/hpc/project/CIBM/plink/phen$j.chr$i >> phen$j.chr$i.job
qsub  phen$j.chr$i.job


cd /gpfs/home/guosa/hpc/project/pmrp/Exom2/imputation
for i in {1..22} 
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo bcftools view -i \'R2\>0.9\' chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.vcf.gz>> chr$i.job
echo convert2annovar.pl -format vcf4 chr$i.dose.filter.vcf.gz  \> ../annovar/chr$i.dose.vcf.avinput >> chr$i.job
echo annotate_variation.pl -out ../annovar/chr$i.annovar.txt -build hg19 ../annovar/chr$i.dose.vcf.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ >> chr$i.job
qsub chr$i.job
done


open my F, '<:gzip' shift @ARGV;
while( <F> ) {
next if /#/;
my @geno;
my @line=split/\s+/;
$geno[0]=$line[3];
$geno[1]=$line[4];
my ($A1,$A2)=split/\|/,$line[9];
print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$geno[$A1]\t$line[$A2]\n";
}	


tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 16:56995835-57017756 > genotypes.vcf
vcftools --vcf genotypes.vcf --plink-tped --out plinkformat
plink --tped plinkformat.tped --tfam plinkformat.tfam --make-bed --out delete
plink --bfile delete --ld rs11485349 rs1823549
plink --bfile delete --r2 --ld-snp-list list.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ld_results
plink --bfile delete --r2 --ld-snp-list list.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --recode HV --out haploview

AHRR: rs169000804, rs11960760

tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 5:304292-438405 > AHRR.vcf
plink --vcf AHRR.vcf --maf 0.01 --make-bed --out AHRR
plink --bfile AHRR --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out AHRR
plink --bfile AHRR --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --recode HV --out AHRR

AHRR: rs169000804, rs11960760

tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 5:138373-543772 > AHRR.vcf
plink --vcf AHRR.vcf --maf 0.1 --make-bed --out AHRR
plink --bfile AHRR --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out AHRR
plink --bfile AHRR --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --recode HV --out AHRR


rs11949577; rs11960760
rs13166205; rs727000694


	
# VCF to PGEN
for i in {1..22} X Y
do
plink2 --vcf chr$i.recode.vcf.gz --make-pgen --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --out ./pgen/chr$i
done

chr1     102906142       102906481   rs11485349.rs1823549.CHB.CHS.r.log
plink2 --pfile chr1 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --ld rs11485349 rs1823549 --
plink2 --pfile chr1 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt --ld rs11485349 rs1823549 --export haps --out rs11485349
plink2 --pfile chr1 --keep ~/hpc/db/hg19/1000Genome/CHB_CHS_221.txt  --keep-nosex --snps rs11485349,rs1823549 --export haps --out rs11485349


/gpfs/home/guosa/hpc/db/hg19/1000Genome
/gpfs/home/guosa/hpc/db/hg19/1000Genome/plink

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

# 1000Genome Sample
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150
awk '{cmd="plink --bfile ~/hpc/db/hg19/1000Genome/plink/"$1 " --ld "$4" "$5 " --out './ld/'"$4"."$5".1000G.r | qsub -N "$4"."$5 " -e ./temp/ -o ./temp/";system(cmd)}'  WGBS.bed.cSNP150.bed.sort.bed.pair.bed

vcftools --gzvcf ~/hpc/db/hg19/1000Genome/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --snp WGBS.SNP.txt

grep R-sq *log | awk -F'[=\sD]' '$5>0.1{print}' 
grep R-sq *log | awk -F'[=\sD]' '{print $7}' 

hpc/hemochromatosis/plink
two_alof_all_combed_v4.phen
plink --file FGF6-C11 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out  FGF6-C11
--allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov --linear --out test

cd /gpfs/home/guosa/hpc/hemochromatosis/plink
for i in {1..22}
do
for j in {1..19}
do
echo \#PBS -N chr$i  > phen$j.chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> phen$j.chr$i.job
echo cd $(pwd) >> phen$j.chr$i.job
echo plink --file exomechip_SNV_PASS_BEAGLE_chr$i\_phased_sel2 --allow-no-sex --pheno two_alof_all_combed_v4.phen --mpheno $j --assoc --out /gpfs/home/guosa/hpc/project/CIBM/plink/phen$j.chr$i >> phen$j.chr$i.job
qsub  phen$j.chr$i.job
done
done

cd /gpfs/home/guosa/hpc/project/CIBM/plink/
for j in {1..19}
do
cat phen$j.chr*assoc >> phen$j.pvalue.txt
done
awk '{print "chr"$1,$3-1,$3,$2,$4,$5,$6,$7,$10,$9}' OFS="\t" > phen$j.pvalue.bed
done


#!/bin/bash
if [ ! -f beagle.28Sep18.793.jar ]; then
  echo
  echo "Downloading beagle.28Sep18.793.jar"
  wget http://faculty.washington.edu/browning/beagle/beagle.28Sep18.793.jar
fi
if [ ! -f bref3.28Sep18.793.jar ]; then
  echo
  echo "Downloading bref3.28Sep18.793.jar"
  wget http://faculty.washington.edu/browning/beagle/bref3.28Sep18.793.jar
fi
echo

if [ ! -f test.28Sep18.793.vcf.gz ]; then
    echo
    echo "*** Downloading some 1000 Genomes Project data to file: test.28Sep18.793.vcf.gz ***"
    wget http://faculty.washington.edu/browning/beagle/test.28Sep18.793.vcf.gz
fi

if [ ! -f conform-gt.24May16.cee.jar ]; then
    echo
    echo "*** Downloading some 1000 Genomes Project data to file: test.28Sep18.793.vcf.gz ***"
    wget http://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar
fi

wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/20140625_related_individuals.txt
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_male_samples_v3.20130502.ALL.panel
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel

wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/filterlines.jar
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/gtstats.jar
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/remove.ids.jar
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/simplify-vcf.jar

grep EUR integrated_call_male_samples_v3.20130502.ALL.panel | cut -f1 > EUR.excl

java -jar conform-gt.24May16.cee.jar ref=/home/local/MFLDCLIN/guosa/hpc/tools/beagle/chr1.1kg.phase3.v5a.vcf.gz gt=chr1.vcf.gz chrom=1 out=mod.chr1 


cd /gpfs/home/guosa/hpc/project/pmrp/phase2/phase
for i in {1..24}
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=32 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo plink --bfile S_Hebbring_Unr.Guo --recode vcf --chr $i --snps-only just-acgt --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i >> chr$i.job
qsub chr$i.job
done

for i in {1..24}
do
plink --bfile S_Hebbring_Unr.Guo --recode vcf --chr $i --snps-only just-acgt --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i
done


cd /gpfs/home/guosa/hpc/tools/beagle/
for i in {1..24}
do
plink --vcf chr$i.1kg.phase3.v5a.vcf.gz --recode vcf --keep EUR.incl --make --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i
done


cd /gpfs/home/guosa/hpc/project/pmrp/phase2/phase/beagle
for i in {1..24}
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=32 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo java -jar conform-gt.24May16.cee.jar ref=/gpfs/home/guosa/hpc/tools/beagle/chr$i.1kg.phase3.v5a.vcf.gz gt=S_Hebbring_Unr.Guo.Forward.chr$i.vcf chrom=$i out=Mod.S_Hebbring_Unr.Guo.chr$i >> chr$i.job
qsub chr$i.job
done


echo
echo "*** Creating test files: ref.28Sep18.793.vcf.gz target.28Sep18.793.vcf.gz ***"
echo
zcat test.28Sep18.793.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > ref.28Sep18.793.vcf.gz
zcat test.28Sep18.793.vcf.gz | cut -f1-9,191-200 | gzip > target.28Sep18.793.vcf.gz

echo
echo "*** Running test analysis with \"gt=\" argument ***"
echo
java -jar beagle.28Sep18.793.jar gt=test.28Sep18.793.vcf.gz out=out.gt

echo
echo "*** Running test analysis with \"ref=\" and \"gt=\" arguments ***"
echo
java -jar beagle.28Sep18.793.jar ref=ref.28Sep18.793.vcf.gz gt=target.28Sep18.793.vcf.gz out=out.ref

echo
echo "*** Making \"bref3\" file ***"
echo
java -jar bref3.28Sep18.793.jar ref.28Sep18.793.vcf.gz > ref.28Sep18.793.bref3

echo
echo "*** Running test analysis with \"bref3\" file ***"
echo
java -jar beagle.28Sep18.793.jar ref=ref.28Sep18.793.bref3 gt=target.28Sep18.793.vcf.gz out=out.bref3



awk '{print "chr"$1,$3-1,$3,$2}' OFS="\t" ASA-MD.map > ASA-MD.hg19.txt
awk '{print "chr"$2,$3-1,$3}' OFS="\t" ASA-CHIA.sites.txt > ASA-CHIA.hg19.txt
bedtools intersect -wao -a ASA-CHIA.hg19.txt -b ~/hpc/db/hg19/snp150.hg19.txt > ASA-CHIA.hg19.snp150.txt


plink2 --bfile PMRP.PhaseII.Steven.Guo.RA --pca approx  --maf 0.05 --memory 40000 --threads 32 --out PMRP.PhaseII.Steven.Guo.RA.pca
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip.fam --extract rs7550295.txt --recode --tab --out test


/home/guosa/hpc/project/SGLT1

cd /home/guosa/hpc/rheumatology/RA/miRNASNP/
for i in chr{1..22} chrX chrY
do
plink --bfile ~/hpc/db/hg19/1000Genome/plink/$i --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --extract Target_92-SNP.txt --freq --out $i
done



for i in chr{1..22} chrX chrY
do
cat $i.frq >> FRQ.txt
done

plink --bfile ~/hpc/db/hg19/1000Genome/chr5 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --r2 --ld-window 800 --ld-window-r2 0 --ld-snps rs10050860,rs2287987 --out Test

plink --bfile ~/hpc/db/hg19/1000Genome/chr5 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --extract two.test --recode fastphase --out two 

plink --bfile ~/hpc/db/hg19/1000Genome/plink/chr1 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --extract PADI3.SNP.txt --freq --out PADI3


chr1:17,387,360-17,731,458

UK biobank
ALL of US  storage 

# plan A:  
cd /home/guosa/hpc/rheumatology/RA/miRNASNP/planB
bedtools intersect -wo -a ~/hpc/db/hg38/miRNA.hg38.bed -b  ~/hpc/db/hg38/allSNP150.GRCH38.bed > miRNA.hg38.allSNP150.bed
awk '{print $8}' miRNA.hg38.allSNP150.bed | sort -u > miRNA.hg38.allSNP150.list.txt

for i in chr{1..22} chrX chrY
do
plink --bfile ~/hpc/db/hg19/1000Genome/plink/$i --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --extract miRNA.hg38.allSNP150.list.txt --maf 0.01 --freq --out $i
done

rm Freq.txt
cat *frq >> Freq.txt
awk '$5>0.01' Freq.txt | grep -v CHR > miRNA.hg38.allSNP150.maf.list.txt
data<-read.table("miRNA.hg38.allSNP150.maf.list.txt")
r1<-which(nchar(as.character(data$V4))>1)
r2<-which(nchar(as.character(data$V3))>1)
newdata=data[-unique(c(r1,r2)),]
write.table(newdata,file="miRNA.hg38.allSNP150.maf.list.uni.snp.txt",sep="\t",quote=F,row.names=F,col.names=F)


# Go to CHG1
cd /home/guosa/hpc/rheumatology/RA/miRNASNP/planC
bedtools intersect -wo -a ~/hpc/db/hg38/miRNA.mature.hg38.bed -b ~/hpc/db/hg38/allSNP150.GRCH38.bed > miRNA.hg38.allSNP150.bed
awk '{print $8}' miRNA.hg38.allSNP150.bed | sort -u > miRNA.hg38.allSNP150.list.txt
for i in chr{1..22} chrX chrY
do
plink --bfile ~/hpc/db/hg19/1000Genome/plink/$i --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --extract miRNA.hg38.allSNP150.list.txt --maf 0.01 --freq --out $i
done
rm Freq.txt
cat *frq >> Freq.txt
awk '$5>0.01' Freq.txt | grep -v CHR > miRNA.hg38.allSNP150.maf.list.txt

data<-read.table("miRNA.hg38.allSNP150.maf.list.txt")
r1<-which(nchar(as.character(data$V4))>1)
r2<-which(nchar(as.character(data$V3))>1)
r3<-which(data$V5<0.05)
newdata=data[-unique(c(r1,r2,r3)),]
write.table(newdata,file="miRNA.hg38.allSNP150.maf.list.uni.snp.txt",sep="\t",quote=F,row.names=F,col.names=F)

# SNP to bed with UCSC table
bedtools window -w 500000 -a miRNA.hg38.allSNP150.maf.list.uni.snp.bed -b ../AutoImmue.GWAS.SNP.hg38.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.hg38.allSNP150.maf.list.uni.snp.bed -b ../RA.GWAS.SNP.hg38.bed | awk '{print $4}' | sort -u | wc -l


bedtools window -v -w 500000 -a miRNA.hg38.allSNP150.maf.list.uni.snp.bed -b ../AutoImmue.GWAS.SNP.hg38.bed > xxx


data<-read.table("miRNA.hg38.allSNP150.maf.list.txt")
r1<-which(nchar(as.character(data$V4))>1)
r2<-which(nchar(as.character(data$V3))>1)
r3<-which(data$V5<0.05)
newdata=data[-unique(c(r1,r2,r3)),]
write.table(newdata,file="miRNA.hg38.allSNP150.maf0.05.N75.list.uni.snp.txt",sep="\t",quote=F,row.names=F,col.names=F)


plink --bfile ~/hpc/db/hg19/1000Genome/chr5 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --r2 --extrac --ld-snps rs10050860,rs17482078,rs2287987,rs2549803,rs27044,rs30187 --out ERAP1
plink --bfile ~/hpc/db/hg19/1000Genome/chr6 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --r2 --extrac --ld-snps rs116488202,rs13202464,rs3819299,rs4349859 --out HLA-B27
plink --bfile ~/hpc/db/hg19/1000Genome/chr1 --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --r2 --extrac --ld-snps rs1004819,rs10489629,rs10889677,rs11209026,rs11465804,rs1343151 --out IL23R

plink --bfile ~/hpc/db/hg19/1000Genome/chr6 --r2 --ld-snps rs116488202,rs13202464,rs3819299,rs4349859 --out Test






browser position chr6:32457705-32995816
track type=bedGraph name=CpG-SNP description=CpG_SNP Marshfield Clinic visibility=full color=0,0,255
chr10   14743   14744   2
chr10   15390   15391   2
chr10   21624   21625   1
chr10   26658   26659   5
chr10   27821   27822   5
chr10   28940   28941   5

	
cd /home/guosa/hpc/db/hg38/fa/chroms
for i in {1..22} X Y
do
bedtools window -w 5000 -a chr$i.CpGSnp.bed -b chr$i.CpGSnp.bed > chr$i.count
perl cpgsnp2bedgraph.pl chr$i.count 
bedtools sort chr$i.count.bedgraph > chr$i.count.sort.bedgraph
echo $i
done


cd /home/guosa/hpc/db/hg38/fa/chroms
for i in {1..22} X Y
do
bedtools sort -i chr$i.count.bedgraph > chr$i.count.sort.bedgraph
echo $i
done


cd ~/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd $(pwd) >> chr$i.job
echo perl fa2mask.pl chr$i >> chr$i.job
qsub chr$i.job
done


cd ~/hpc/db/hg38/fa/chroms
for i in {14..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=32 >> $i.job
echo cd $(pwd) >> chr$i.job
echo perl cpgsnp.pl chr$i >> chr$i.job
qsub chr$i.job
done

cd /home/guosa/hpc/db/hg38/fa/chroms
for i in {1..13} X Y
do
echo chr$i
perl episnp.pl chr$i
done

cd /home/guosa/hpc/db/hg38/fa/chroms
for i in {14..22}
do
echo chr$i
perl episnp.pl chr$i
done

Sequence alignment of human FGF subfamilies. The sequences were aligned with ClustalX and Dendroscope [41,42].


# vcf2plink
cd ~/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo plink --vcf ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --make-bed --out ./plink/chr$i >> chr$i.job
qsub chr$i.job
done

bedtools intersect -wo -a AutoImmue.GWAS.SNP.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > AutoImmue.GWAS.SNP.hg38.commonSNP.bed

bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $8}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $16}' | sort -u | wc -l

bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $8}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $16}' | sort -u | wc -l

bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $1,$2,$3,$4,$8,$16}' OFS="\t" | sort -u > miRNA.64.AutoImmue.GWAS.SNP.hg38.commonSNP.bed



# CHB+CHS 221 Sample
cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
awk '{cmd="plink --bfile ~/hpc/db/hg19/1000Genome/plink/"$1 " --keep '~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt' --ld "$5" "$6 " --out './ld/'"$5"."$6".CHB.CHS.r | qsub -N "$5"."$6 " -e ./temp/ -o ./temp/";system(cmd)}'  miRNA.64.AutoImmue.GWAS.SNP.hg38.commonSNP.bed

cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
awk '{cmd="plink --bfile ~/hpc/db/hg19/1000Genome/plink/"$1 " --ld "$5" "$6 " --out './ld/'"$5"."$6".CHB.CHS.r | qsub -N "$5"."$6 " -e ./temp/ -o ./temp/";system(cmd)}'  miRNA.64.AutoImmue.GWAS.SNP.hg38.commonSNP.bed


plink --bfile ~/hpc/db/hg19/1000Genome/chrX --keep ~/hpc/rheumatology/RA/miRNASNP/CHB_CHS_221.txt --ld rs782045310 rs782045310 --out ./ld/rs782045310.rs782045310.CHB.CHS.r 


cd /gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP
R
data<-read.table("immune_disease.txt",head=T,sep="\t")
data<-read.table("RA.GWAS.SNP.txt",head=F,sep="\t")
grep("chr",data$V4)
snp150<-read.table("~/hpc/db/hg19/snp150.hg19.txt")

bedtools window -w 500000 -a ~/hpc/db/hg38/miRNA.hg38.bed -b ~/hpc/rheumatology/RA/miRNASNP/RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $4}' | sort -u | wc -l

bedtools intersect -wo -a miRNA.hg38.commonSNP.bed -b miRNA.hg38.commonSNP.bed > AutoImmue.GWAS.SNP.hg38.commonSNP.bed
  --ld rs3745199 rs2233152

grep rs3745199 ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
grep rs2233152 ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
grep rs3745199 ~/hpc/db/hg19/1000Genome/chr19.bim
grep rs2233152 ~/hpc/db/hg19/1000Genome/chr19.bim
 


awk '{print $4}' AutoImmue.GWAS.SNP.hg38.commonSNP.bed
awk '{print $8}' AutoImmue.GWAS.SNP.hg38.commonSNP.bed | sort -u | wc -l


perl -lane 'print $1 if /(miR-\d+-*\d*\w*)/ig' RA-miRNA.txt | sort -u > RA_108_miRNA.txt
tr '[:upper:]' '[:lower:]' < RA_108_miRNA.txt > RA_miRNA.txt
sort -u RA_miRNA.txt > RA_miRNA.uni.txt
 
perl -lane 'print $1 if /(miR-\d+)/ig' RA-miRNA.txt | sort -u > RA_108_miRNA.txt
tr '[:upper:]' '[:lower:]' < RA_108_miRNA.txt > RA_miRNA.txt
sort -u RA_miRNA.txt > RA_miRNA.uni.txt

 
/home/guosa/hpc/db/hg19/hg19.win2k.bed
/home/guosa/hpc/epimarker/bedgraph
/home/guosa/hpc/db/hg38/window200

wget http://epigenomesportal.ca/tracks/DEEP/hg19/62968.DEEP.01_HepG2_LiHG_Ct1.WGB-Seq.methylation_profile.bigWig ./
wget http://epigenomesportal.ca/tracks/DEEP/hg19/63053.DEEP.41_Hf01_LiHe_Ct.WGB-Seq.methylation_profile.bigWig ./
wget http://epigenomesportal.ca/tracks/DEEP/hg19/63078.DEEP.41_Hf02_LiHe_Ct.WGB-Seq.methylation_profile.bigWig ./
wget http://epigenomesportal.ca/tracks/DEEP/hg19/63099.DEEP.41_Hf03_LiHe_Ct.WGB-Seq.methylation_profile.bigWig ./

for i in `ls *bigWig`
do
cp $i $i.bw
done

cd /gpfs/home/guosa/run/bedmethyl
for i in `ls *bw`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $i.job
echo bigWigAverageOverBed $i /gpfs/home/guosa/hpc/nash/methcov/4521NASH.candidate.hg38.bed /gpfs/home/guosa/hpc/nash/methcov/PBMC/$i.tab >> $i.job
echo $i.job
qsub $i.job
done



AHRR: chr5:317989-325025
cd /home/guosa/hpc/methylation/methtype
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr5 --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --freq --out /home/guosa/hpc/methylation/methtype/AHRR
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr5 --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --freq --out /home/guosa/hpc/methylation/methtype/AHRR

cd /home/guosa/hpc/db/hg19/plan2/commonCpGSNP
rm cpgSNP.hg19.uni.bed
for i in `ls *cpgsnp.bed`
do
perl unitrim.pl $i >> cpgSNP.hg19.uni.bed
done
awk '{print $4}' cpgSNP.hg19.uni.bed > cpgSNP.hg19.list.txt 

cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr$i --extract /home/guosa/hpc/db/hg19/plan2/commonCpGSNP/cpgSNP.hg19.list.txt --make-bed --out chr$i.cpgSNP
done

cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
perl assignRefence.pl chr$i.cpgSNP.bim > chr$i.cpgSNP.bim.tr
mv chr$i.cpgSNP.bim.tr chr$i.cpgSNP.bim
awk '{print $2"\tC"}' chr$i.cpgSNP.bim > chr$i.C.ref
plink --bfile chr$i.cpgSNP --reference-allele cpgSNP.ref.txt --make-bed --out chr$i.cpgSNP.C
done

plink --bfile /home/guosa/hpc/cpgSNP/chr5.cpgSNP.C --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --reference-allele cpgSNP.ref.txt  --freq --out /home/guosa/hpc/cpgSNP/AHRR
cd /home/guosa/hpc/cpgSNP
perl swith2csv.pl fastphase_hapguess_switch.out > AHRR.hap.txt

plink --bfile FGF6 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out data2
plink --bfile FGF6 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq case-control --allow-no-sex --out FGF6-C2
data<-read.table("AHRR.txt",row.names=1)
colSums(input)[colSums(input)>50]
awk '{print $1}' two_alof_all_combed_v4.phen > two.phen.txt
awk '{print $1}' FGF6.fam > one.phen.txt
diff one.phen.txt two.phen.txt

MC.154314@1075641920 MC.154314@1075641920 0 0 0 -9
MC.154316@1075680154 MC.154316@1075680154 0 0 0 -9

data<-read.table("two_alof_all_combed_v4.phen",head=T)
output<-data.frame(FID=data$FID,FID=data$FID,0,0,0,PheTyp7_Iron_C1=data$PheTyp7_Iron_C1)
write.table(output,file="FGF6-C1.phen",sep="\t",quote=F,col.names=F,row.names=F)
data<-read.table("FGF6.dip.txt",head=F,row.names=1)
saminfo<-read.table("two_alof_all_combed_v4.phen",head=T)

levels(data[,1])[4]=0
levels(data[,2])[3]=0
levels(data[,1])[2]="CTCGGCCCAC"
levels(data[,1])[7]="M1"
levels(data[,2])[6]="M1"
levels(data[,1])[5]="M4"
levels(data[,2])[4]="M4"

levels(data[,1]) %in% levels(data[,2])
levels(data[,1])[2]<-"CTCGGCCCAC"
table(data[,1],data[,2])


write.table(table(data[,1],data[,2]),file="FGF6.PMRP.txt",sep="\t",col.names=NA,row.names=T,quote=F)

CON1<-which(saminfo$PheTyp7_Iron_C1==1)
CAS1=which(saminfo$PheTyp7_Iron_C1==2)
write.table(table(data[CAS1,1],data[CAS1,2]),file="FGF6.PMRP.C1.Case.tsv",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(table(data[CON1,1],data[CON1,2]),file="FGF6.PMRP.C1.CONTROL.tsv",sep="\t",col.names=NA,row.names=T,quote=F)

CON2=which(saminfo$PheTyp7_Iron_C2==1)
CAS2=which(saminfo$PheTyp7_Iron_C2==2)
write.table(table(data[CAS2,1],data[CAS2,2]),file="FGF6.PMRP.C2.Case.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(table(data[CON2,1],data[CON2,2]),file="FGF6.PMRP.C2.CONTROL.txt",sep="\t",col.names=NA,row.names=T,quote=F)

tab1<-table(data[CAS2,1],data[CAS2,2])
tab2<-table(data[CON2,1],data[CON2,2])

table(C3)
table(C4)
C4[C4$V2=="M1",]


/home/guosa/hpc/wgbs/GSE17972
/home/guosa/hpc/methylation/methtype

for i in `ls *.cov.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd $(pwd) >> ${i}.job
echo gunzip $i >> ${i}.job
qsub ${i}.job
done


mkdir ../mf
for i in `ls *.cov.gz`
do
j=${i/_L001_R1_001_00_bismark_bt2_pe.bam/}
BismarkRefereDb="~/hpc/db/hg19/bismrk/"
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=32 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd $(pwd) >> ${i}.job
echo samtools sort -n $i $i.nsort >> ${i}.job
echo $j
done


wigToBigWig GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.wig ~/hpc/db/hg19/hg19.chrom.sizes GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.bw 
bigWigAverageOverBed GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.bw ~/hpc/db/hg19/BUR.GRCH37.bed GSM916049_BI.Adult_Liver.Bisulfite-Seq.3.bw.tab

# cov2bedgraph
# cov2bw
for i in `ls *.cov`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo awk \'{print \$1,\$2-1,\$3,\$4}\' $i OFS=\"\\t\" \>$i.bedgraph >>$i.job
echo wigToBigWig $i.bedgraph ~/hpc/db/hg19/hg19.chrom.sizes $i.bedgraph.bw >> $i.job
echo $i.job
done

# bw2tab
for i in `ls GSM*bw`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo awk \'{print \$1,\$2-1,\$3,\$4}\' $i OFS=\"\\t\" \>$i.bedgraph >>$i.job
echo bigWigAverageOverBed $i ~/hpc/db/hg19/BUR.GRCH37.bed $i.tab >> $i.job
echo $i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/nash/methcov
for i in `ls *bw`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo bigWigAverageOverBed $i ~/hpc/db/hg19/LINE1.hg19.bed ./LINE-1/$i.tab >> $i.job
echo $i.job
qsub $i.job
done


for i in `ls *bigWig.bw`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo bigWigAverageOverBed $i ~/hpc/db/hg19/BUR.GRCH37.bed $i.tab >> $i.job
echo $i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/nash/methcov
for i in `ls *bw`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/methcov >> $i.job
echo bigWigAverageOverBed $i ~/hpc/db/hg19/BUR.GRCH37.bed $i.tab >> $i.job
echo $i.job
qsub $i.job
done




for i in `ls *txt`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=32 >> $i.job
echo cd /gpfs/home/guosa/hpc/wgbs/GSE17972 >> $i.job
echo ./qmap2bedgraph $i \>$i.bedgraph >> $i.job
echo $i.job
qsub $i.job
done

for i in `ls *txt`
do
./qmap2bedgraph $i > $i.bedgraph
done


Unbiased RNAi screen for hepcidin regulators links hepcidin suppression to proliferative Ras/RAF and nutrient-dependent mTOR signaling.

awk '{print $1,$2,$3,$1":"$2"-"$3}' ~/hpc/db/hg19/BUR.GRCH37.bed > ~/hpc/db/hg19/BUR.GRCH37..bed
mv ~/hpc/db/hg19/BUR.GRCH37..bed ~/hpc/db/hg19/BUR.GRCH37.bed 

If the above is true, it makes perfect sense. Let me explain, so mutant FGF6 have decreased expression of HAMP which in turn leads to decreased hepcidin. Hepcidin controls the ubiquitination of ferroportin (iron exporter which is very important in the transport of Fe2+ out of a cell). So, if you have lower HAMP--> lower hepcidin--> decreased breakdown of ferroportin--> increased export of ferrous iron--> decreased intracellular Fe2+. 

it make sense for small intestine cell. However, we are using liver cell, we 



cd /gpfs/home/guosa/hpc/hemochromatosis/plink/fgf6

/gpfs/home/guosa/hpc/tools/phase.2.1.1.linux

--recode fastphase



/gpfs/home/guosa/hpc/db/hg19/allsnp150.hg19
grep rs148135963 ~/hpc/db/hg19/allsnp150.hg19

chr12   4553360 4553361 rs148135963     +       G       C/G


plink --file exomechip_SNV_PASS_BEAGLE_chr12_phased_sel2 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq case-control --out ./fgf/FGF6

plink --file exomechip_SNV_PASS_BEAGLE_chr12_phased_sel2 --chr chr12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --out FGF6


12      rs71583765      0       4543445 T       C
12      rs7961645       0       4543487 A       T
12      rs372503617     0       4543494 A       C
12      rs140796118     0       4543523 A       G
12      rs140216440     0       4553300 A       G
12      rs145168026     0       4553304 T       C
12      rs148135963     0       4553361 C       G
12      rs148657794     0       4554550 T       C
12      rs11613495      0       4554630 G       A
12      rs139098855     0       4554651 A       C

plink --bfile FGF6 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out data2
plink --bfile FGF6 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq case-control --allow-no-sex --out FGF6-C2


awk '{print $1}' two_alof_all_combed_v4.phen > two.phen.txt
awk '{print $1}' FGF6.fam > one.phen.txt
diff one.phen.txt two.phen.txt

MC.154314@1075641920 MC.154314@1075641920 0 0 0 -9
MC.154316@1075680154 MC.154316@1075680154 0 0 0 -9
MC.154321@1075678065 MC.154321@1075678065 0 0 0 -9
MC.154332@1075679003 MC.154332@1075679003 0 0 0 -9
MC.154345@1075680506 MC.154345@1075680506 0 0 0 -9
MC.154346@1076253986 MC.154346@1076253986 0 0 0 -9

data<-read.table("two_alof_all_combed_v4.phen",head=T)
output<-data.frame(FID=data$FID,FID=data$FID,0,0,0,PheTyp7_Iron_C1=data$PheTyp7_Iron_C1)
write.table(output,file="FGF6-C1.phen",sep="\t",quote=F,col.names=F,row.names=F)


data<-read.table("FGF6.dip.txt",head=F,row.names=1)
saminfo<-read.table("two_alof_all_combed_v4.phen",head=T)

levels(data[,1])[4]=0
levels(data[,2])[3]=0
levels(data[,1])[2]="CTCGGCCCAC"
levels(data[,1])[7]="M1"
levels(data[,2])[6]="M1"
levels(data[,1])[5]="M4"
levels(data[,2])[4]="M4"

levels(data[,1]) %in% levels(data[,2])
levels(data[,1])[2]<-"CTCGGCCCAC"
table(data[,1],data[,2])


write.table(table(data[,1],data[,2]),file="FGF6.PMRP.txt",sep="\t",col.names=NA,row.names=T,quote=F)

CON1<-which(saminfo$PheTyp7_Iron_C1==1)
CAS1=which(saminfo$PheTyp7_Iron_C1==2)
write.table(table(data[CAS1,1],data[CAS1,2]),file="FGF6.PMRP.C1.Case.tsv",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(table(data[CON1,1],data[CON1,2]),file="FGF6.PMRP.C1.CONTROL.tsv",sep="\t",col.names=NA,row.names=T,quote=F)

CON2=which(saminfo$PheTyp7_Iron_C2==1)
CAS2=which(saminfo$PheTyp7_Iron_C2==2)
write.table(table(data[CAS2,1],data[CAS2,2]),file="FGF6.PMRP.C2.Case.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(table(data[CON2,1],data[CON2,2]),file="FGF6.PMRP.C2.CONTROL.txt",sep="\t",col.names=NA,row.names=T,quote=F)

tab1<-table(data[CAS2,1],data[CAS2,2])
tab2<-table(data[CON2,1],data[CON2,2])

table(C3)
table(C4)
C4[C4$V2=="M1",]




data<-read.table("two_alof_all_combed_v4.phen",head=T)
output<-data.frame(FID=data$FID,FID=data$FID,0,0,0,PheTyp7_Iron_C1=data$PheTyp7_Iron_C2)
write.table(output,file="FGF6-C2.phen",sep="\t",quote=F,col.names=F,row.names=F)

plink --bfile FGF6-C11 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out  FGF6-C11
plink --bfile FGF6-C11 --chr 12 --from-bp 4536137 --to-bp 4561951 --freq case-control --recode fastphase --allow-no-sex --out  FGF6-C11
plink --file FGF6-C11 --hap myfile.hlist --hap-phase

plink --bfile FGF6-C2 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq counts --allow-no-sex --out data1
plink --bfile FGF6-C2 --chr 12 --from-bp 4536137 --to-bp 4561951 --make-bed --freq case-control --allow-no-sex --out FGF6-C21
plink --bfile FGF6-C21 --chr 12 --from-bp 4536137 --to-bp 4561951 --freq case-control --recode fastphase --allow-no-sex --out  FGF6-C21




CHR             SNP   A1   A2        MAF_A        MAF_U  NCHROBS_A  NCHROBS_U
  12       exm975531    T    C            0      0.00323         38      13930
  12       exm975536    A    T            0    0.0007897         38      13930
  12   variant.22576    A    C            0    7.179e-05         38      13930
  12       exm975538    0    G            0            0         38      13930
  12       exm975541    0    G            0            0         38      13930
  12       exm975542    0    C            0            0         38      13930
  12       exm975549    C    G            0    0.0002872         38      13930
  12       exm975571    T    C            0    0.0006461         38      13930
  12       exm975575    G    A       0.3947        0.148         38      13930
  12       exm975576    0    C            0            0         38      13930
cannot open file '/home/guosa/hpc/Schrodi_2ALOF/result/PheTyp7_Iron_C1_res_summary_v4.csv': No such file or directory
time (tabix -R snp_30000.bed ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz > /dev/null)
( tabix -R snp_30000.bed  > /dev/null; )  26,14s user 0,68s system 82% cpu 32,342 total
$ time (bcftools view -R snp_30000.bed ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bcf.gz > /dev/null)                                                 
( bcftools view -R snp_30000.bed ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bcf.gz  > /dev/null; )  19,30s user 0,68s system 75% cpu 26,340 total
$ time (bcftools view -Ou -R snp_30000.bed ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bcf.gz > /dev/null)
( bcftools view -Ou -R snp_30000.bed ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bcf.gz > /dev/null; )  16,97s user 0,21s system 99% cpu 17,190 total

cd /gpfs/home/guosa/nc/plink
awk '{print $1,2*$2-$7,2*$2-$6,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.mirror.hg19.bed 
bedtools intersect -wao -a DMER.mirror.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.mirror.hg19.SNP150.bed
awk '{print $1,$6,$7,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.hg19.bed 
bedtools intersect -wao -a DMER.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.hg19.SNP150.bed
mkdir temp
mkdir LD

for i in CEU CHB YRI
do

 perl summ.pl CEU
 perl summ.pl CHB
 perl summ.pl YRI

 
# Whole Sample 
cd /gpfs/home/guosa/nc/plink
awk '{print $1,2*$2-$7,2*$2-$6,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.mirror.hg19.bed 
bedtools intersect -wao -a DMER.mirror.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.mirror.hg19.SNP150.bed
awk '{print $1,$6,$7,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.hg19.bed 
bedtools intersect -wao -a DMER.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.hg19.SNP150.bed
mkdir temp
mkdir LD
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

 
# CHB Sample
cd /gpfs/home/guosa/nc/plink
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CHB.r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CHB.r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

# CEU Sample
cd /gpfs/home/guosa/nc/plink
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CEU.r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CEU.r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

# YRI Sample
cd /gpfs/home/guosa/nc/plink
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".YRI.r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".YRI.r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed





o conf mbuildpl_arg "--install_base /gpfs/home/guosa/hpc/tools/perl"
o conf makepl_arg "PREFIX=/gpfs/home/guosa/hpc/tools/perl"
install List::Util

o conf mbuildpl_arg "--install_base /gpfs/home/guosa/hpc/tools/perl"
o conf makepl_arg "REFIX=/gpfs/home/guosa/hpc/tools/perl"
install List::Util
 
cd /gpfs/home/guosa/nc/plink
awk '{print $1,2*$2-$7,2*$2-$6,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.mirror.hg19.bed 
bedtools intersect -wao -a DMER.mirror.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.mirror.hg19.SNP150.bed
awk '{print $1,$6,$7,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.hg19.bed 
bedtools intersect -wao -a DMER.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.hg19.SNP150.bed
mkdir temp
mkdir LD
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed


awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './LD/'"$4"."$8".r2 | qsub -N "$4"."$8;print(cmd)}' DMER.mirror.hg19.SNP150.bed
join -t $'\t' -1 1 -2 2 <(sort -t $'\t' -k1,1 input.txt) <(sort -t $'\t' -k2,2 ref2.txt) | uniq | awk -F '\t' '{line=sprintf("%s\t%s\t%s\t%s\t%s",$1,$2,$3,$4,$5);if($7>=$2 && $7<=$3) {a[line]+=int($6);} else {a[line]+=0;}} END {for(line in a) printf("%s\t%d\n",line,a[line]);}'
plink --bfile ~/hpc/db/hg19/1000Genome/chrX --ld rs12863738 rs5931090 --out ./LD/rs12863738.rs5931090.r2 | qsub -N rs12863738.rs5931090 -e ./temp/
# vcf2gz
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo bgzip -c chr$i.recode.vcf  \> chr$i.recode.vcf.gz >> chr$i.job
echo tabix -p vcf chr$i.recode.vcf.gz >> chr$i.job
qsub chr$i.job
done

# vcf2plink
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo plink --vcf chr$i.recode.vcf --make-bed --out chr$i >> chr$i.job
qsub chr$i.job
done


# vcf2gz
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo bgzip -c ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf  \> ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> chr$i.job
echo tabix -p ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  >> chr$i.job
qsub chr$i.job
done


# vcf2gz
cd /gpfs/home/guosa/nc/Shuffle
for i in `ls *.ld`
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo perl ldtrim.pl $i >> chr$i.job
qsub chr$i.job
done

for i in `ls *.ld`
do
echo \#PBS -N chr$i  > $i.job
echo cd $(pwd) >> $i.job
echo perl ../ldc.pl $i >> $i.job
qsub $i.job
done

for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --maf 0.01 --remove-indels --recode --out chr$i >> chr$i.job
qsub chr$i.job
done

plink --vcf ~/hpc/db/hg19/1000Genome/chr6.uni.vcf --extract chr6.rs6457620.txt --r inter-chr dprime --out rs6457620

for i in `ls *bed.sort.bed`
do
bedtools closest -a GWAS-378.GRCH37.sort.bed -b $i > GWAS-378.$i.distance.bed
done


scp *vcf.gz nu_guos@128.105.244.191:/mnt/gluster/nu_guos/1000Genome/

cd /gpfs/home/guosa/hpc/db/hg19/CpGs
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo awk \'{print \$1\"\\t\"\$2-1\"\\t\"\$2}\' chr$i.CpG.positions.txt \> /gpfs/home/guosa/hpc/db/hg19/CpGs/bed/chr$i.CpG.positions.bed >> chr$i.job
echo bedtools closest -a ./bed/chr$i.CpG.positions.bed -b /gpfs/home/guosa/hpc/db/hg19/refGene.hg19.bed \> ./bed/chr$i.CpG.positions.refGene.hg19.bed >> chr$i.job
qsub chr$i.job
done

cd /gpfs/home/guosa/hpc/db/hg19/CpGs/bed
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$8}\' chr$i.CpG.positions.refGene.hg19.bed \> /gpfs/home/guosa/hpc/db/hg19/CpGs/bed/chr$i.CpG.Gene.hg19.bed >> chr$i.job
echo gzip chr$i.CpG.Gene.hg19.bed  >> chr$i.job
qsub chr$i.job
done

for i in `ls *CpG.Gene.hg19.bed`
do
gzip $i
done



bedtools shuffle -i RA-OA.DMER.GRCH37.bed -g ~/hpc/db/hg19/hg19.chrom.sizes | bedtools sort -i -


comm listA listB -1 -2 -3

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/GWAS-reported-gene.txt
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/RA_FLS_DEG_11vs11_Nature_Commnication.txt


a=/gpfs/home/guosa/hpc/db/hg19/refGene.hg19.bed
b=/gpfs/home/guosa/nc/GWAS-RA-378.GRCH37.bed
bedtools window -a $b -b $a -w 500000 | awk '{print $9}' | sort -u  >  GWAS-378SNP-50K-Gene.txt
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/RA_FLS_DEG_11vs11_Nature_Commnication.txt
awk '{print $1}' RA_FLS_DEG_11vs11_Nature_Commnication.txt | sort -u > DMER.GENE.List.txt
comm GWAS-378SNP-50K-Gene.txt  DMER.GENE.List.txt -1 -2

ACOXL
ACYP1
ARHGEF3-AS1
BCL2L11
BIRC3
CDH6
CFTR
CREB5
FRAS1
GPR157
LINC00565
LPIN2
MINK1
NUDT6
RTTN
SLC2A5
TYR


## background distribution for shuffling of DMERs and LD

for i in $(seq 0 0.05 1)
do 
echo \#PBS -N $i > $i.job
echo \#PBS -o ./temp/ >>$i.job
echo \#PBS -e ./temp/ >> $i.job
echo cd $(pwd) >> $i.job
echo perl ../ALL2000/checklinkage.pl $i \| awk \'{print \$1}\' \| sort -u \|wc -l \> "tres.R2."$i".txt" >> $i.job
qsub $i.job
done

for i in $(seq 0 0.05 1)
do 
echo \#PBS -N $i > $i.job
echo \#PBS -o ./temp/ >>$i.job
echo \#PBS -e ./temp/ >> $i.job
echo cd $(pwd) >> $i.job
echo perl ../ALL2000/checklinkageDP.pl $i \| awk \'{print \$1}\' \| sort -u \|wc -l \> "tres.DP."$i".txt" >> $i.job
qsub $i.job
done


for i in $(seq 0 0.05 1)
do 
perl checklinkage.pl $i | awk '{print $1}' | sort -u |wc -l 
done


time plink --vcf ~/hpc/db/hg19/1000Genome/chr1.uni.vcf --extract chr1.rs12131057.txt --r inter-chr dprime --out ../test
time plink --vcf ~/hpc/db/hg19/1000Genome/chr1.uni.vcf --extract chr1.rs12131057.txt --out ../test


--r inter-chr dprime



# Calculate LD 
cd /gpfs/home/guosa/nc/Snpset
for k in ALL2000 CEU CHB YRI 
do
mkdir $k
for i in DMER
do 
for j in `ls *.txt.uni`
do
pop=~/hpc/db/hg19/1000Genome/$k.txt
vcf=/gpfs/home/guosa/nc/chr$j.dmer.recode.vcf
echo \#PBS -N $j.$k  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo plink --vcf $vcf --keep $pop --extract $j --r2 inter-chr dprime --out ./$k/$i.$j.ld  >> $i.$j.job
qsub $i.$j.job
done
done
done


cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
vcftools --gzvcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --snps /gpfs/home/guosa/nc/GWAS-RA-378.SNPlist.txt  --recode --recode-INFO-all  --out GWAS-RA-378.Genome1000


# bedtools sort
cd /gpfs/home/guosa/nc/gwassnpsampling
for i in `ls *.txt.bed`
do
echo \#PBS -N chr$i > chr$i.job
echo cd $(pwd) >> chr$i.job
echo bedtools sort -i $i \> $i.sort >> chr$i.job
echo bedtools closest -a $i.sort -b  ~/hpc/db/hg19/RA-OA.DMER.GRCH37.bed \> $i.DMER.dis.bed >>chr$i.job
qsub chr$i.job
done

cd /gpfs/home/guosa/nc/gwassnpsampling
for i in `ls *.txt.bed`
do
echo \#PBS -N chr$i > chr$i.job
echo cd $(pwd) >> chr$i.job
echo bedtools sort -i $i \> $i.sort >> chr$i.job
echo bedtools closest -a $i.sort -b  ~/hpc/db/hg19/RA-OA.DMER.GRCH37.bed \> $i.DMER.dis.bed >>chr$i.job
qsub chr$i.job
done



# vcf2gz
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo bgzip -c chr$i.uni.vcf  \> chr$i.vcf.gz >> chr$i.job
echo tabix -p vcf chr$i.vcf.gz >> chr$i.job
qsub chr$i.job
done


cd /gpfs/home/guosa/nc
mkdir ./intersect_between_dmer_and_gwas
cp RA-OA.DMER.GRCH37.bed ./intersect_between_dmer_and_gwas/
cp GWAS-RA-378.GRCH37.bed ./intersect_between_dmer_and_gwas/
cd ./intersect_between_dmer_and_gwas/
mkdir ./Shuffle/
mkdir ./temp/
for i in `ls *RA-OA.DMER.GRCH37.bed`
do
for j in {1..10000}
do
echo \#PBS -N $i.$j  > $i.$j.job
echo \#PBS -o ./temp/ >>$i.$j.job
echo \#PBS -e ./temp/>> $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo bedtools shuffle -i $i -g ~/hpc/db/hg19/hg19.chrom.sizes \> ./Shuffle/$i.$j.shuffle >> $i.$j.job
echo bedtools sort -i ./Shuffle/$i.$j.shuffle  \> ./Shuffle/$i.$j.shuffle.sort  >> $i.$j.job 
echo bedtools closest -a ./Shuffle/$i.$j.shuffle.sort -b GWAS-RA-378.GRCH37.bed \> ./Shuffle/$i.$j.GWAS.Cloest >> $i.$j.job
qsub $i.$j.job
done
done



## Statistic 
bedtools closest -a GWAS-RA-378.GRCH37.bed -b > RA-OA.DMER-GWAS.N378.bed
## R 
file="RA-OA.DMER-GWAS.N378.bed"
data<-read.table(file)
xx1<-abs(as.numeric(as.character(data[,2]))-as.numeric(as.character(data[,7])))
xx2<-abs(as.numeric(as.character(data[,3]))-as.numeric(as.character(data[,7])))
xx2[which(xx1<xx2)]=xx1[which(xx1<xx2)]
len<-c(mean(xx2,na.rm=T))
sum<-c(sum(xx2<50000,na.rm=T)/length(xx2))
# R
setwd("./")
file=list.files(pattern="*GWAS.Cloest")
Length<-c()
Sum<
for(i in 1:length(file)){
  marker<-unlist(strsplit(file[i],"[.]"))[1]
  data<-read.table(file[i])
  xx1<-abs(as.numeric(as.character(data[,2]))-as.numeric(as.character(data[,7])))
  xx2<-abs(as.numeric(as.character(data[,3]))-as.numeric(as.character(data[,7])))
  xx2[which(xx1<xx2)]=xx1[which(xx1<xx2)]
  Length<-c(Length,mean(xx2,na.rm=T))
  Sum<-c(Sum,sum(xx2<50000,na.rm=T)/length(xx2))
}
LengthS1000000<-sample(Length,100000,replace=T)



bedtools intersect -wa  -a GWAS-SNP-359.GRCH37.bed -b ~/hpc/db/hg19/RA-OA.DMER.GRCH37.bed  | awk '{print $4}' | sort -u | wc -l
bedtools sort -i GWAS-RA-378.GRCH37.bed > GWAS-RA-378.GRCH37.bed.sort

bedtools sort -i ~/hpc/db/hg38/RA-OA.DMER.GRCH38.bed > ~/hpc/db/hg38/RA-OA.DMER.GRCH38.bed.sort

bedtools closest -a GWAS-RA-378.GRCH38.bed -b ~/hpc/db/hg38/RA-OA.DMER.GRCH38.bed.sort > RA-OA.DMER.RA-GWAS.distance.bed
bedtools closest -b GWAS-RA-378.GRCH38.bed -a ~/hpc/db/hg38/RA-OA.DMER.GRCH38.bed.sort > RA-GWAS.RA-OA.DMER.distance.bed

bedtools intersect -wao -a ~/hpc/db/hg38/RA-OA.DMER.GRCH38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed | awk '($9 !="."){print $1"\t"$2"\t"$3"\t"$9}' > RA-OA.DMER.GRCH38.SNP150.bed
bedtools closest -a  RA-OA.DMER.GRCH38.SNP150.bed -b GWAS-RA-378.GRCH38.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > RA-OA.DMER.GRCH37.SNP150.RA-GWAS.bed
awk '{print $4}'  RA-GWAS-RA-OA.DMER.GRCH37.SNP150.bed | sort -u | wc -l
awk '{print $5}'  RA-GWAS-RA-OA.DMER.GRCH37.SNP150.bed | sort -u | wc -l
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' RA-GWAS-RA-OA.DMER.GRCH37.SNP150.bed > RA-GWAS-RA-OA.DMER.GRCH37.SNP150.ReV.bed
awk '{print $4}'  RA-GWAS-RA-OA.DMER.GRCH37.SNP150.ReV.bed | sort -u | wc -l
awk '{print $5}'  RA-GWAS-RA-OA.DMER.GRCH37.SNP150.ReV.bed | sort -u | wc -l
cat RA-GWAS-RA-OA.DMER.GRCH37.SNP150.ReV.bed RA-OA.DMER.GRCH37.SNP150.RA-GWAS.bed > RA-OA-DMER-GWAS.GRCH37.bed
bedtools sort -i RA-OA-DMER-GWAS.GRCH37.bed > RA-OA-DMER-GWAS.GRCH37.sort.bed
sort -u RA-OA-DMER-GWAS.GRCH37.sort.bed > RA-OA-DMER-GWAS.GRCH37.sort.uni.bed
awk '{print $4}'  RA-OA-DMER-GWAS.GRCH37.sort.uni.bed | sort -u | wc -l
awk '{print $5}'  RA-OA-DMER-GWAS.GRCH37.sort.uni.bed | sort -u | wc -l
awk '{print $4"\n"$5}' RA-OA-DMER-GWAS.GRCH37.sort.uni.bed | sort -u > DMER.GWAS.SNP.FullList.txt

# prepare DMER-GWAS-1000G-Dataset
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --vcf chr$i.uni.vcf --recode --recode-INFO-all --snps /gpfs/home/guosa/nc/DMER.GWAS.SNP.FullList.txt --out /gpfs/home/guosa/nc/newVcf/chr$i.dmer >> chr$i.job
qsub chr$i.job
done

# creat SNP-set for fast LD-calculation
awk '{print $4"\n"$5 > "./snpset/"$1"."$4".snp.list"}' RA-GWAS-RA-OA.DMER.GRCH37.SNP150.bed
cd /gpfs/home/guosa/nc/snpset
for i in `ls *.list`
do
sort -u $i > $i.uni
done

# Calculate LD 
cd /gpfs/home/guosa/nc/Snpset
for k in ALL2000 CEU CHB YRI 
do
mkdir $k
for i in DMER
do 
for j in `ls *.txt.uni`
do
pop=~/hpc/db/hg19/1000Genome/$k.txt
vcf=/gpfs/home/guosa/nc/chr$j.dmer.recode.vcf
echo \#PBS -N $j.$k  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo plink --vcf $vcf --keep $pop --extract $j --r2 inter-chr dprime --out ./$k/$i.$j.ld  >> $i.$j.job
qsub $i.$j.job
done
done
done

# Calculate LD 






awk '{print $4}' RA-OA.DMER.GRCH37.SNP150.RA-GWAS.bed | sort -u | wc -l
bedtools closest -a GWAS-RA-378.GRCH38.bed -b RA-OA.DMER.GRCH38.SNP150.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > RA-OA.DMER.GRCH37.SNP150.RA-GWAS.bed
bedtools sort -i RA-OA.DMER.GRCH37.SNP150.RA-GWAS.bed > RA-OA.DMER.GRCH37.SNP150.RA-GWAS.sort.bed
sort -u RA-OA.DMER.GRCH37.SNP150.RA-GWAS.sort.bed > RA-OA.DMER.GRCH37.SNP150.RA-GWAS.sort.uni.bed
awk '{print $4}'  RA-OA.DMER.GRCH37.SNP150.RA-GWAS.sort.uni.bed | sort -u | wc -l
awk '{print $5}'  RA-OA.DMER.GRCH37.SNP150.RA-GWAS.sort.uni.bed | sort -u | wc -l

awk '{print $4"\n"$5 > "./snpset/$1.$3.txt}' RA-GWAS-RA-OA.DMER.GRCH37.SNP150.bed
awk '{print $1,$2 > ""$1""$4".txt"}' RA-GWAS-RA-OA.DMER.GRCH37.SNP150.bed

rm 
cd /gpfs/home/guosa/nc/snpset
perl ldcalc.pl

cd /gpfs/home/guosa/nc/
for k in CEU CHB YRI 
do
mkdir $k
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X
do
pop=~/hpc/db/hg19/1000Genome/$k.txt
vcf=/gpfs/home/guosa/nc/chr$j.dmer.recode.vcf
echo \#PBS -N $i.$j.$k  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo plink --vcf $vcf --keep $pop --r inter-chr dprime --out ./$k/$i.chr$j.ld  >> $i.$j.job
qsub $i.$j.job
done
done
done

for i in $(seq 0 0.05 1)
do 
echo \#PBS -N $i > $i.job
echo \#PBS -o ./temp/ >>$i.job
echo \#PBS -e ./temp/ >> $i.job
echo cd $(pwd) >> $i.job
echo perl linked.pl $i \| wc -l \> "tres."$i".txt" >> $i.job
qsub $i.job
done

for j in {1..1}
do
for i in `ls *bed.sort.bed `
do
echo \#PBS -N $i.$j  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo \#PBS -o ./temp/ >>$i.$j.job
echo \#PBS -e ./temp/>> $i.$j.job
echo bedtools shuffle -i $i -excl wgEncodeDukeMapabilityRegionsExcludable.bed -g ~/hpc/db/hg19/hg19.chrom.sizes \> ./Shuffle/$i.$j.shuffle >> $i.$j.job
echo bedtools sort -i ./Shuffle/$i.$j.shuffle  \> ./Shuffle/$i.$j.shuffle.sort  >> $i.$j.job 
echo bedtools closest -a ./Shuffle/$i.$j.shuffle.sort -b RA_GWAS_475_Catalog_GRCH37.bed \> ./Shuffle/$i.$j.GWAS.Cloest >> $i.$j.job
qsub $i.$j.job
done
done

bedtools shuffle -i ~/hpc/db/hg38/RA-OA.DMER.GRCH38.bed.sort -excl wgEncodeDukeMapabilityRegionsExcludable.bed -g ~/hpc/db/hg38/hg38.chrom.sizes > RA-OA.DMER.GRCH38.Shuffle.bed

cd /gpfs/home/guosa/nc/SNP
for i in `ls *snp.list.uni` 
do
set -f  
chr=(${i//./ })

awk -F. '{print $1}' $i
done


${i}

# creat SNP-set
awk '{print $4"\n"$5 > "./SNP/"$1"."$5".snp.list"}' dmer.uni.bed
for i in `ls *.list`
do
sort -u $i > $i.uni
done


hr16   68378232        68378422        rs16957831      rs16973500      NA
chr16   69253591        69253780        rs8052096       rs16973500      NA
chr16   69968969        69969249        rs73580166      rs16973500      NA
chr16   69968969        69969249        rs115740373     rs16973500      NA
chr16   69968969        69969249        rs77253126      rs16973500      NA
chr16   69968969        69969249        rs114735867     rs16973500      NA
chr16   77911799        77912060        rs11150022      rs78507369      NA
chr16   78329523        78329717        rs140326318     rs78507369      NA
chr16   78477479        78478001        rs12445747      rs78507369      NA
chr16   78477479        78478001        rs7199023       rs78507369      NA
chr16   78477479        78478001        rs74030257      rs78507369      NA
chr16   78477479        78478001        rs72799912      rs78507369      NA




bedtools closest -i /gpfs/home/guosa/hpc/db/hg19/OA.DMER.GRCH37.bed -b RA_GWAS_475_Catalog_GRCH37.bed



cd /gpfs/home/guosa/nc/
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in CHB YRI 
do
perl ldcal.pbs.pl $i $j
done
done


cd /gpfs/home/guosa/nc/
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in CEU CHB YRI 
do
echo \#PBS -N $j.$i >$i.$j.job
echo \#PBS nodes=1:ppn=32 >>$i.$j.job
echo \#PBS -o ./temp/ >>$i.$j.job
echo \#PBS -e ./temp/>> $i.$j.job
echo  cd $(pwd) >>$i.$j.job
echo perl ldsum.pl $i $j >>$i.$j.job
qsub $i.$j.job
done
done



use Cwd;
use strict;
my $dir=getcwd;
chdir $dir;
my $marker=shift @ARGV;
my $pop=shift @ARGV;
open F,"/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/$marker.bed.sort.bed.SNP150.pair.bed" || die "cannot open $marker\n";
while(<F>){
chomp;
my ($chr,$start,$end,$rs1,$rs2)=split/\s+/;
open OUT,">$marker.$pop.$rs1.$rs2.job";
print OUT "#PBS -N $marker.$pop.$rs1.$rs2\n";
print OUT "#PBS -o ./temp/\n";
print OUT "#PBS -e ./temp/\n";
print OUT "cd $dir\n";
print OUT "plink --vcf /gpfs/home/guosa/nc/$chr.dmer.recode.vcf --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/$pop.txt --ld $rs1 $rs2 --out ./$pop/$marker.$pop.$rs1.$rs2\n";
print OUT "rm ./$pop/*nosex\n";
close OUT;
system("qsub $marker.$pop.$rs1.$rs2.job");
system("rm $marker.$pop.$rs1.$rs2.job");
}
close F;




# CHS ACB ASW BEB CDX CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI 
cd /gpfs/home/guosa/nc/
for k in CEU CHB YRI 
do
mkdir $k
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
for j in {1..22} X
do
pop=~/hpc/db/hg19/1000Genome/$k.txt
vcf=/gpfs/home/guosa/nc/chr$j.dmer.recode.vcf
echo \#PBS -N $i.$j.$k  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo plink --vcf $vcf --keep $pop --r inter-chr dprime --out ./$k/$i.chr$j.ld  >> $i.$j.job
qsub $i.$j.job
done
done
done

plink --vcf /gpfs/home/guosa/nc/chrX.dmer.recode.vcf --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt --r2 inter-chr --ld-window-r2 0 --out ./WGBS.chrX.ld

plink --vcf /gpfs/home/guosa/nc/chr15.dmer.recode.vcf --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt --ld rs139622314 rs12719740

 
H3K4me1.chr15.ld.ld
 
 
chr15   102035946       102042500       rs11247303      rs12719740
chr15   102035946       102042500       rs7179692       rs12719740
chr15   102035946       102042500       rs7181883       rs12719740
chr15   102035946       102042500       rs2898886       rs12719740
chr15   102035946       102042500       rs895395        rs12719740
chr15   102035946       102042500       rs735504        rs12719740
chr15   102035946       102042500       rs746428        rs12719740
chr15   102035946       102042500       rs735505        rs12719740
chr15   102035946       102042500       rs4965891       rs12719740
chr15   102035946       102042500       rs1992730       rs12719740
chr15   102035946       102042500       rs1992729       rs12719740
chr15   102035946       102042500       rs5815019       rs12719740
chr15   102035946       102042500       rs74347782      rs12719740
chr15   102250827       102251385       rs139622314     rs12719740
chr15   102250827       102251385       rs76748706      rs12719740


 
 
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do
for j in CEU CHB YRI
do
perl ldcalc.pl $i $j
done
done



for i in {12..22} X Y
do
for j in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do
qsub $j.$i.job
doneqsta
done





#CEU
cd ~/hpc/rheumatology/RA/NatureCommunication
for k in CEU CHB CHS ACB ASW BEB CDX CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI 
do
rm -rf $k
done


cd /gpfs/home/guosa/nc
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo bgzip -c chr$i.uni.vcf \> chr$i.uni.vcf.gz >> chr$i.job
echo tabix -p vcf chr$i.uni.vcf.gz >> chr$i.job
qsub chr$i.job
done



cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
awk '{print $4}' /gpfs/home/guosa/hpc/db/hg19/commonsnp150.hg19.bed > commonsnp150.rs.list.tmp
sort -u   commonsnp150.rs.list.tmp >  commonsnp150.rs.list
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --recode --recode-INFO-all --snps ~/hpc/rheumatology/RA/NatureCommunication/dmer.snp.list --out /gpfs/home/guosa/nc/chr$i.common150 >> chr$i.job
#qsub chr$i.job
done


cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
awk '{print $4}' /gpfs/home/guosa/hpc/db/hg19/commonsnp150.hg19.bed > commonsnp150.rs.list.tmp
sort -u   commonsnp150.rs.list.tmp >  commonsnp150.rs.list
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --recode --recode-INFO-all --snps dmer.snp.list --out chr$i.dmer >> chr$i.job
qsub chr$i.job
done

scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/wangmh/data/* ./

sftp guosa_ftp@sftp.mfldclin.edu

sshpass -p 'H6s=!KC=7Sy^'  autoimmune.txt scp guosa_ftp@sftp.mfldclin.edu:/test

sftp guosa_ftp@sftp.mfldclin.edu

cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcf2dedup ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \> chr$i.uni.vcf >> chr$i.job
qsub chr$i.job
done


wget -e robots=off -nH -nd  -r -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17972/suppl/

awk '{print $1"\t"$2}' 1000GenomeSampleInfo.txt >>$3.txt
awk '{print $1,$2 > ""$3".txt"}' 1000GenomeSampleInfo.txt

awk '{print $1,$2 > ""$3".txt"}' 1000GenomeSampleInfo.txt
awk '{ split($2, a, "_"); print $1"\t"a[2]"\t"$3 >> a[1]".txt"; }' foo.txt


############################################
cd /gpfs/home/guosa/hpc/wgbs
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo gunzip GSE17972_HUMtg5lib.qmap.chr$i.txt.gz >> chr$i.job
qsub chr$i.job
done
############################################

rs75289190
rs2843401
rs61768776
rs72632736

1       13413938        rs200269882     C       T       100     PASS    AC=93;
1       13634735        rs200269882     C       T       100     PASS  
grep rs200269882 chr1.uni.vcf


chr1    1735021 1735199 rs75289190      rs2843401
chr1    4104130 4104228 rs61768776      rs72632736

plink --vcf /gpfs/home/guosa/hpc/db/hg19/1000Genome/chr1.uni.vcf --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt --r2 --ld-snp-list TEST.pairSNP --out ./CEU/TEST.chr1.ld


###########################
grep rs147777396 WGBS.pairSNP
grep rs147777396 WGBS.bed.sort.bed.SNP150.pair.bed
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do
for j in CEU JPT CHINA
do
perl ldcalc.pl $i $j
done
done
###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo gunzip ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz >> chr$i.job
qsub chr$i.job
done
###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo bcftools norm -d both --threads=4 ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O z -o chr$i.vcf.gz >> chr$i.job
qsub chr$i.job
done
############################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo ./vcf2dedup chr$i.common150.recode.vcf \> chr$i.uni.vcf >> chr$i.job
qsub chr$i.job
done

for i in 617{27..50}.bright
do
qdel $i
done
############################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo bcftools view --threads=4 -k --types=snps chr$i.common150.recode.vcf -O z -o chr$i.uni.vcf.gz >> chr$i.job
#qsub chr$i.job
done
###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
awk '{print $4}' /gpfs/home/guosa/hpc/db/hg19/commonsnp150.hg19.bed > commonsnp150.rs.list.tmp
sort -u   commonsnp150.rs.list.tmp >  commonsnp150.rs.list
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --recode --recode-INFO-all --snps commonsnp150.rs.list --out chr$i.common150 >> chr$i.job
#qsub chr$i.job
done
###########################
for i in `ls *bed.sort.bed `
do
echo \#PBS -N $i  > $i.job
echo cd $(pwd) >>$i.job
echo bedtools intersect -wao -a $i -b  ~/hpc/db/hg19/commonsnp150.hg19.bed \> $i.SNP150 >> $i.job
qsub $i.job
done
###########################
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication
for i in `ls *SNP150`
do
bedtools closest -a $i -b RA_GWAS_475_Catalog_GRCH37.bed | awk '($10!="." && $11!="."){print $1"\t"$2"\t"$3"\t"$10"\t"$18}' > $i.pair.bed  &
echo $i
done
###########################
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
awk '{print $4"\n"$5}' $i.bed.sort.bed.SNP150.pair.bed | sort -u > $i.pairSNP &
done

##########################
for i in `ls *bed.sort.bed`
do
wc -l $i
done

for i in `ls *bed.sort.bed`
do
bedtools window -b $i -a RA_GWAS_475_Catalog_GRCH37.bed -w 500000 | awk '{print $4}' |grep rs | sort -u | wc -l 
done


for i in $(seq 0 5000 500000)
do
bedtools window -b DMER.sort.bed -a RA_GWAS_475_Catalog_GRCH37.bed -w $i | awk '{print $4}' |grep rs | sort -u | wc -l 
done
for i in $(seq 0 5000 500000)
do
echo $i
done

data<-read.table("LenvsPerc.txt",sep="\t")
pdf("Supp.Fig-4.pdf")
plot(data[,2]~data[,1],cex=2,pch=16,main="Percentage of sauration vs distance",xlab="Distance",ylab="percentage",lwd=2)
dev.off()

###########################
for j in {1..1000}
do
for i in `ls *bed.sort.bed `
do
echo \#PBS -N $i.$j  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo bedtools shuffle -i $i -excl wgEncodeDukeMapabilityRegionsExcludable.bed -g ~/hpc/db/hg19/hg19.chrom.sizes \> ./Shuffle/$i.$j.shuffle >> $i.$j.job
echo bedtools sort -i ./Shuffle/$i.$j.shuffle  \> ./Shuffle/$i.$j.shuffle.sort  >> $i.$j.job 
echo bedtools closest -a ./Shuffle/$i.$j.shuffle.sort -b RA_GWAS_475_Catalog_GRCH37.bed \> ./Shuffle/$i.$j.GWAS.Cloest >> $i.$j.job
qsub $i.$j.job
done
done
###########################



###########################
#CEU
cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication

for k in CEU CHB CHS ACB ASW BEB CDX CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI 
do
mkdir $k
for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS 
do 
for j in {1..22} X Y
do
echo \#PBS -N $i.$k  > $i.$j.job
echo cd $(pwd) >> $i.$j.job
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
pop=~/hpc/db/hg19/1000Genome/$k.txt
vcf=~/hpc/db/hg19/1000Genome/chr$j.uni.vcf
echo plink --vcf $vcf --keep $pop --r2 --ld-snp-list $i.pairSNP --out ./$k/$i.chr$j.ld >> $i.$j.job
qsub $i.$j.job
done
done
done

TEST.pairSNP

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
###########################


rs10523519

plink --vcf /gpfs/home/guosa/hpc/db/hg19/1000Genome/chr11.uni.vcf --snp rs10523519 --recode tab --out rs10523519
rs10523519.log

bcftools norm -d both --threads=2 ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O z -o chr$i.vcf.gz >> chr$i.job

grep rs10627346 /gpfs/home/guosa/hpc/db/hg19/1000Genome/chr20.common150.recode.vcf
###########################

for i in `ls *.bed `
do
bedtools sort -i $i > $i.sort.bed
echo $i
done

for i in `ls *sort.bed `
do
bedtools closest -a $i -b RA_GWAS_475_Catalog_GRCH37.bed > $i.GWAS.Cloest &
done



for i in H3K27AC H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3 OpenChrom WGBS
do 
dmer_snp_db=/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt
vcf=~/hpc/db/hg19/1000Genome/
plink --vcf $vcf --keep $ceu --extract $dmer_snp_db --make-bed --out ./plink/DMER.GRCH37
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz  -O chrX.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi -O chrX.vcf.gz.tbi
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi
mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz  ALL.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi ALL.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
mv ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz ALL.chrY.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
mv ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi ALL.chrY.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz ALL.chrMT.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 
mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi ALL.chrMT.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz

cd /mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication

for i in `ls *.vcf.gz`
do
chr=echo ${string//[[:digit:]]/}
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/db/1000Genome >> ${i}.job
echo ulimit -n 3000 >> ${i}.job
echo plink --vcf $i --snps-only --out $i >>  ${i}.job
echo ${i}.job
qsub $i.job
done


cd ~/hpc/rheumatology/RA/NatureCommunication/snp150/
dmer_snp_db=~/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=~/hpc/db/hg19/1000Genome/CEU.txt
dir=~/hpc/db/hg19/1000Genome


for i in chr{1..22} chrX chrY
do
plink --vcf $dir/ --keep $ceu --ld-snp-list --make-bed --out ./plink/DMER.GRCH37
done


use strict;
open F,$shift @ARGV;
while(<F>){
chomp;
my @line=split/\s+/;
print "$line[0]\t$line[1]\t$line[2]\n";
}

for i in `ls *bed`
do
perl trim.pl $i > $i.trim
done


for i in `ls *bed`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication >> ${i}.job
echo bedtools intersect -wao -a $i -b ~/hpc/db/hg19/commonsnp150.hg19.bed  \> $i.cSNP150.bed >>  ${i}.job
echo ${i}.job
qsub $i.job
done





for i in `ls *.bed`
do
awk '{print $1,$2,$3}' $i | sort -u | wc -l 
done

for i in `ls *cSNP150.bed`
do
awk '$7!="."{print $1,$2,$3}' $i | sort -u | wc -l 
done


for i in `ls *cSNP150.bed`
do
bedtools sort -i $i > $i.sort.bed 
done

for i in `ls *.sort.bed `
do
bedtools closest -a $i -b ../RA_GWAS_475_Catalog_GRCH37.bed  | awk '($10!="." && $11!="."){print $1"\t"$2"\t"$3"\t"$10"\t"$18}' > $i.pair.bed &
echo $i
done

plink --bfile mydata --ld rs2840528 rs7545940

rm ./plink/*
for i in chr{1..22} chrX chrY
do
for j in `ls *.pair.bed`
do
grep  -w $i $j | awk '{print $4"\t"$5}'>> ./plink/$j.$i.plink.bed
done
done


rm ./plink/*
for i in chr{1..22} chrX chrY
do
for j in `ls *.pair.bed`
do
grep  -w $i $j >> ./plink/$j.$i.plink.bed
done
done


cd /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink/

for i in chr{1..22} chrX  chrY
do
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink/$i --extract $dmer_snp_db --make-bed --snps-only --out $i.dmer
done


for i in `ls *.cSNP150.bed`
do 
awk '{print $10}' $i | grep rs >> SNP.db
done

85,056,574 

dmer_snp_db=/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/snp150/SNP.uni.db
ceu=/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt
for i in chr{1..22} chrX  chrY
do
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink/$i --keep $ceu --extract $dmer_snp_db --make-bed --out $i.dmer
done

for f in ../ALL.chr*.vcf.gz
do
f2=${f##*/}       
fout=${f%.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz} 
echo \#PBS -N $fout  > $fout.job
echo \#PBS -l nodes=1:ppn=16 >> $fout.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $fout.job
echo \#PBS -m abe  >> $fout.job
echo \#PBS -q longq  >> $fout.job
echo plink --vcf $f --make-bed --snps-only --keep $ceu --extract $dmer_snp_db --out $fout  >> $fout.job
done

 
 
plink --bfile chr1 --list-duplicate-vars ids-only suppress-first --out chr1.uni

plink --bfile chr1 --list-duplicate-vars --out chr1.uni
plink --bfile chr1 --exclude chr1.uni.dupvar --out chr1.uni


grep rs10625753 chr1.uni.dupvar

plink --bfile chr1 --r2 --ld-snp-list ~hpc/rheumatology/RA/NatureCommunication/snp150/plink/H3K27AC.bed.cSNP150.bed.sort.bed.pair.bed.chr1.plink.bed/


rs60778745      rs11121380
plink --bfile chr1 --r2 --ld-snp-list ~hpc/rheumatology/RA/NatureCommunication/snp150/plink/H3K27AC.bed.cSNP150.bed.sort.bed.pair.bed.chr1.plink.bed/

H3K27AC.bed.cSNP150.bed.sort.bed.pair.bed.chr1.plink.bed 
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink
plink --bfile chr1 --ld-snp-list /mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication/snp150/plink/H3K27AC.bed.cSNP150.bed.sort.bed.pair.bed.chr1.plink.bed
cd /mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication/snp150/plink
grep rs10625753 /mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication/snp150/plink/H3K27AC.bed.cSNP150.bed.sort.bed.pair.bed.chr1.plink.bed
bedtools intersect -wao -a WGBS.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed  > WGBS.bed.cSNP150.bed
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
CollapsABEL: Generalized CDH (GCDH) Analysis
wget -e robots=off -nH -nd  -r -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29127/suppl/
wget -e robots=off -nH -r -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83147/suppl/
GSE83nnn
cd /gpfs/home/guosa/hpc/rheumatology/RA/lncRNA
bedtools intersect -wao -a candidate-87.bed -b /gpfs/home/guosa/hpc/db/miRNASNP/Gwas.Catalog.immune.disease.GRCH38.bed > candidate.Guo.bed
cd /gpfs/home/guosa/hpc/rheumatology/RA/lncRNA
bedtools intersect -wao -a candidate-112.bed -b /gpfs/home/guosa/hpc/db/miRNASNP/Gwas.Catalog.immune.disease.GRCH38.bed > candidate.Guo.112.GWAS.bed
bedtools intersect -wao -a RA-OA.DMER.GRCH38.bed -b /gpfs/home/guosa/hpc/db/miRNASNP/Gwas.Catalog.immune.disease.GRCH38.bed > RA-OA.DMER.GWAS.GRCH38.bed
59/103
less RA-OA.DMER.GWAS.GRCH38.bed | grep rheuma | awk '{print $9}' | sort -u  | wc -l
scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/wangmh/data/* ./

scp autoimmune.txt guosa_ftp@sftp.mfldclin.edu:/test


git clone https://github.com/dnanexus/htslib.git
git clone https://github.com/dnanexus/samtools.git
make -C samtools
tar -czvf samtools.tar.gz samtools

cp Pool_A-plasma_S1_L001_R1_001.nsort.bam ../Pool_A-plasma_S1_L001_R1_001.nsort.bam
cp Pool_A-tissue_S4_L002_R1_001.nsort.bam ../Pool_A-tissue_S4_L002_R1_001.nsort.bam
cp Pool_B-tissue_S5_L002_R1_001.nsort.bam ../Pool_B-tissue_S5_L002_R1_001.nsort.bam
cp Pool_C-tissue_S6_L002_R1_001.nsort.bam ../Pool_C-tissue_S6_L002_R1_001.nsort.bam


GSE29127
wget -e robots=off -nH -nd  -r -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29127/suppl/

wget -e robots=off -nH -nd  -r -nd FTP_URL

FTP_URL is as the following:

ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1001/soft/GDS1001.soft.gz  
ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1001/soft/GDS1001_full.soft.gz  
ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL10/soft/GPL10_family.soft.gz   
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE1/soft/GSE1_family.soft.gz   
ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL10/miniml/GPL10_family.xml.tgz   
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE1/miniml/GSE1_family.xml.tgz   
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE1/matrix/GSE1_series_matrix.txt.gz   
ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1073/suppl/   
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/suppl/GSE1000_RAW.tar   
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1nnn/GSM1137/suppl/GSM1137.CEL.gz  
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1nnn/GSM29nn/suppl/

wget https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29127&targ=self&view=brief&form=text
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE29127&format=file
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/suppl/GSE1000_RAW.tar  

mkdir ../methyfreq
option1=$(echo --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end)
option2=$(echo --bedGraph --ignore 1 --buffer_size 4G --gzip --comprehensive)

for i in `ls Pool_*bam`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/pool/bam >> ${i}.job
echo bismark_methylation_extractor ${option1} ${option2} --output ../methyfreq  ./$i >> ${i}.job
echo ${i}.job
done

mkdir ../methyfreq
option1=$(echo --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end)
option2=$(echo --bedGraph --ignore 1 --buffer_size 4G --gzip --comprehensive)
for i in `ls *bam`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/ >> ${i}.job
echo bismark_methylation_extractor ${option1} ${option2} --output ../methyfreq  ./$i >> ${i}.job
echo ${i}.job
done




qsub Pool_A-plasma_S1_L001_R1_001.bam.job
qsub Pool_B-plasma_S2_L001_R1_001.bam.job
qsub Pool_C-plasma_S3_L001_R1_001.bam.job

qsub Pool_A-tissue_S4_L002_R1_001.nsort.bam.job
qsub Pool_B-tissue_S5_L002_R1_001.nsort.bam.job
qsub Pool_C-tissue_S6_L002_R1_001.nsort.bam.job


Pool_A-plasma_S1_L001_R1_001.nsort.bam.job
mv Pool_B-plasma_S2_L001_R1_001re.nsort.bam.job Pool_B-plasma_S2_L001_R1_001.bam.job 
mv Pool_C-plasma_S3_L001_R1_001re.nsort.bam.job Pool_C-plasma_S3_L001_R1_001.bam.job
mv Pool_C-plasma_S3_L001_R1_001re.nsort.bam Pool_C-plasma_S3_L001_R1_001.bam
mv Pool_B-plasma_S2_L001_R1_001re.nsort.bam Pool_B-plasma_S2_L001_R1_001.bam
mv Pool_A-plasma_S1_L001_R1_001.nsort.bam Pool_A-plasma_S1_L001_R1_001.bam
mv Pool_A-tissue_S4_L002_R1_001.nsort.bam Pool_A-tissue_S4_L002_R1_001.bam
mv Pool_B-tissue_S5_L002_R1_001.nsort.bam Pool_B-tissue_S5_L002_R1_001.bam
mv Pool_C-tissue_S6_L002_R1_001.nsort.bam Pool_C-tissue_S6_L002_R1_001.bam


mkdir ../mf
for i in `ls CpG_context_*txt`
do
j=${i/_L001_R1_001_00_bismark_bt2_pe.bam/}
BismarkRefereDb="~/hpc/db/hg19/bismrk/"
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=32 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd $(pwd) >> ${i}.job
echo samtools sort -n $i $i.nsort >> ${i}.job
echo $j
done

for i in `ls *.bam`
do
j=${i/.bam/}
echo \#PBS -N $i.sort  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/pool/bam >> ${i}.job
echo cd $(pwd) >> ${i}.job
echo samtools sort -n -@ 6 $i $j.nsort >> ${i}.job
echo $i.job
qsub $i.job
done


for i in `ls *.bam`
do
j=${i/.bam/}
echo \#PBS -N $i.sort  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/ >> ${i}.job
echo cd $(pwd) >> ${i}.job
echo samtools sort -n -@ 6 $i $j.nsort >> ${i}.job
echo $i.job
done



#qsub $i.job
echo coverage2cytosine --merge_CpG ../methyfreq/ --genome_folder $BismarkRefereDb -o ../methyfreq/$j >> ${i}.job





bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq /home/shg047/luministuscdata/Pool_B-tissue_S5_L002_R1_001.bam

bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ./Pool_B-tissue_S5_L002_R1_001.bam

samtools sort -n





Genome-wide DNA methylation mapping in breast cancer cells (HCC1954) and normal breast cells (HMEC)
Design: MethylC-Seq on breast cancer HCC1954 and normal breast HMEC. 100 cycles of sequencing on Illumina platform.

scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata3/*_pe.bam ./
scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata/*.bam* ./

for i in `ls *bigWig`
do
bigWigAverageOverBed $i ~/hpc/db/hg19/BUR.GRCH37.bed $i.tab 
done


grep chr20:59005276-59007276 ~/hpc/db/hg19/BUR.GRCH37.bed
grep chr20:59005276-59007276 ~/hpc/db/hg19/BUR.GRCH37.bed



for i in chr{1..22} chrX chrY chrM
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
# echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
# echo \#PBS -m abe  >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedgraph >> $i.job
echo Rscript --vanilla dmr.R $i.tab.matrix.rlt >> $i.job
qsub $i.job
done


/gpfs/home/guosa/hpc/epi-allele/study/sortbam

/gpfs/home/guosa/hpc/epi-allele/study/sortbam

scp BUR.GRCH37.bed nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/db/hg19/
scp BUR.GRCH38.bed nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/db/hg38/


/home/nu_guos/nash/BUR.GRCH37.bed


CpG.GRCH37.positions.txt

/home/nu_guos/db/hg19/cpgs/CpG.GRCH37.positions.txt

-rwxrwxr-x   1 nu_guos nu_guos 9.2G Aug 14 21:26 A02_S3_L001_R1_001_00_bismark_bt2_pe.bam.sorted.bam
-rw-------   1 nu_guos nu_guos 1.4G Aug 15 00:50 A02_S3_L001_R1_001_00_bismark_bt2_pe.bam.sorted.bam.gz



for i in chr{1..22} chrX chrY
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedgraph >> $i.job
echo Rscript --vanilla BUR.R $i.tab.matrix.rlt >> $i.job
qsub $i.job
done


condor_q -hold
ssh-keygen -t rsa
ssh shg047@23.99.137.107 'mkdir -p .ssh'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'cat >> .ssh/authorized_keys'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'chmod 700 ~/.ssh'

cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'cat >> .ssh/authorized_keys2'
cat /home/nu_guos/.ssh/id_rsa.pub | ssh shg047@23.99.137.107 'chmod 700 ~/.ssh/authorized_keys2'
ssh shg047@23.99.137.107
alias chtc="ssh nu_guos@submit-3.chtc.wisc.edu"


./sshpass -p "luminist%%01" scp -r shg047@23.99.137.107:/home/shg047/luministuscdata3/A02_S3_L001_R1_001_00_bismark_bt2_pe.bam.sorted.bam /mnt/gluster/nu_guos
./sshpass -p "luminist%%01" scp -r shg047@23.99.137.107:/home/shg047/luministuscdata3/A02_S3_L001_R1_001_00_bismark_bt2_pe.bam.sorted.bam.bai /mnt/gluster/nu_guos
perl ./bam2methhap.pl mybed.bed A02_S3_L001_R1_001_00_bismark_bt2_pe.bam.sorted.bam bismark hg19.chrom.sizes chr6.CpG.positions.txt


tar -czvf command.tar.gz /home/nu_guos/tools/sshpass-1.05/sshpass  bam2methhap.sh  bam2methhap.pl
/home/yez/Schrodi_2ALOF/rpgmcd 

shg047@23.99.137.107:/home/shg047/haplotype

scp *.txt nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/db/hg19/cpgs/
/squid/nu_guos



\\mcrfnas2\bigdata\Genetic\Projects\shg047\pmrp\phase2\RA\C2
\mcrfnas2\bigdata\Genetic\Projects\shg047\pmrp\phase2\RA\C2
/mnt/bigdata/Genetic/Projects/shg047/pmrp/phase2/RA/C2
/home/guosa/hpc/pmrp/phase2/RA/C2
/home/local/MFLDCLIN/guosa/hpc/pmrp/phase2/RA/C2cd 

pwd | sed 's/\//\\/g' | sed 's/mnt/\\mcrfnas2/g' 


plink --bfile PMRP.PhaseII.Steven.Guo.RA.C1 --recode vcf --out PMRP.PhaseII.Steven.Guo.RA.C1
plink --bfile PMRP.PhaseII.Steven.Guo.RA.C2 --recode vcf --out PMRP.PhaseII.Steven.Guo.RA.C2
plink2 --bfile PMRP.PhaseII.Steven.Guo.RA --pca approx  --maf 0.05 --memory 40000 --threads 32 --out PMRP.PhaseII.Steven.Guo.RA.pca

setwd("/mnt/bigdata/Genetic/Projects/shg047/pmrp/phase2/RA")
eigenvec<-read.table("PMRP.PhaseII.Steven.Guo.RA.pca.eigenvec",head=F)
saminfo<-read.table("~/hpc/pmrp/phase2/S_Hebbring_Release_Sample_Sheet.txt",head=T,sep="\t")
sam<-saminfo[match(as.character(eigenvec$V1),saminfo$Sample_Name),]

pdf("phase2.pca-population.pdf")
Legends<-unique(data.frame(Population=sam$Population,Col=as.numeric(sam$Population)))
plot(eigenvec$V4~eigenvec$V3,cex=0.55,col=as.numeric(sam$Population),pch=as.numeric(sam$Population),xlab="PC1",ylab="PC2")
legend("topleft",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()

threshold=min(eigenvec$V3[which(sam$Population=="African")])
exclude=eigenvec[which(eigenvec$V3>threshold),1]
write.table(data.frame(exclude,exclude),file="PCA.exclude.individal.thres0.0123604.txt",sep="\t",col.names=F,row.names=F,quote=F)

pdf("phase2.pca-gender.pdf")
Legends<-unique(data.frame(Gender=sam$Gender,Col=as.numeric(sam$Gender)))
plot(eigenvec$PC2~eigenvec$PC1,cex=0.55,col=as.numeric(sam$Gender),pch=as.numeric(sam$Gender))
legend("topleft",legend=Legends$Gender,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()

plink --bfile PMRP.PhaseII.Steven.Guo.RA --remove PCA.exclude.individal.thres0.0123604.txt --make-bed --out PMRP.PhaseII.Steven.Guo.RA.CEU

78,75

grep 8536628-1-0224009354 *phen   F       
grep 6584333-1-0238041894 *phen   F
grep 1176608-1-0238062177 *phen   F
grep 4860568-1-0238041718 *phen   F
grep 7016015-1-0224008233 *phen   M
grep 3358239-1-0238094627 *phen   M
grep 6552643-1-0238041362 *phen   F

8536628-1-0224009354    8536628-1-0224009354
6584333-1-0238041894    6584333-1-0238041894
1176608-1-0238062177    1176608-1-0238062177
4860568-1-0238041718    4860568-1-0238041718
7016015-1-0224008233    7016015-1-0224008233
3358239-1-0238094627    3358239-1-0238094627
6552643-1-0238041362    6552643-1-0238041362

cp S_Hebbring_Unr.Guo.Forward.fam ./RA/PMRP.PhaseII.Steven.Guo.RA.fam
cp S_Hebbring_Unr.Guo.Forward.bim ./RA/PMRP.PhaseII.Steven.Guo.RA.bim
cp S_Hebbring_Unr.Guo.Forward.bed ./RA/PMRP.PhaseII.Steven.Guo.RA.bed

plink --bfile PMRP.PhaseII.Steven.Guo.RA --keep RA.SamList.txt --remove PCA.exclude.individal.thres0.0123604.txt --make-bed --out PMRP.PhaseII.Steven.Guo.RA.CEU
perl famupdate.pl
mv PMRP.PhaseII.Steven.Guo.RA.C.fam PMRP.PhaseII.Steven.Guo.RA.fam
mkdir C1
cp PMRP.PhaseII.Steven.Guo.RA.CEU.C1.fam ./C1/PMRP.PhaseII.Steven.Guo.RA.CEU.C1.fam 
cp PMRP.PhaseII.Steven.Guo.RA.CEU.bed ./C1/PMRP.PhaseII.Steven.Guo.RA.CEU.C1.bed 
cp PMRP.PhaseII.Steven.Guo.RA.CEU.bim ./C1/PMRP.PhaseII.Steven.Guo.RA.CEU.C1.bim 
plink --bfile ./C1/PMRP.PhaseII.Steven.Guo.RA.CEU.C1 --recode vcf --out ./C1/PMRP.PhaseII.Steven.Guo.RA.CEU.C1
mkdir C2  
cp PMRP.PhaseII.Steven.Guo.RA.CEU.C2.fam ./C2/PMRP.PhaseII.Steven.Guo.RA.CEU.C2.fam 
cp PMRP.PhaseII.Steven.Guo.RA.CEU.bed ./C2/PMRP.PhaseII.Steven.Guo.RA.CEU.C2.bed 
cp PMRP.PhaseII.Steven.Guo.RA.CEU.bim ./C2/PMRP.PhaseII.Steven.Guo.RA.CEU.C2.bim
plink --bfile ./C2/PMRP.PhaseII.Steven.Guo.RA.CEU.C2 --recode vcf --out ./C2/PMRP.PhaseII.Steven.Guo.RA.CEU.C2

cd /home/guosa/hpc/pmrp/phase2
cp S_Hebbring_Unr.Guo.Forward.fam ./PA/PMRP.PhaseII.Steven.Guo.PA.fam
cp S_Hebbring_Unr.Guo.Forward.bim ./PA/PMRP.PhaseII.Steven.Guo.PA.bim
cp S_Hebbring_Unr.Guo.Forward.bed ./PA/PMRP.PhaseII.Steven.Guo.PA.bed
awk '{print $1,$1}' ./PA/PA.phen > ./PA/PA.SamList.txt
cd ./PA
plink --bfile PMRP.PhaseII.Steven.Guo.PA --keep PA.SamList.txt --remove PCA.exclude.individal.thres0.0123604.txt --make-bed --out PMRP.PhaseII.Steven.Guo.PA.CEU
perl famupdate.pl
mv PMRP.PhaseII.Steven.Guo.RA.C.fam PMRP.PhaseII.Steven.Guo.RA.fam
plink --bfile PMRP.PhaseII.Steven.Guo.PA.CEU --recode vcf --out PMRP.PhaseII.Steven.Guo.PA.CEU



SampleID        GENDER  PheTyp3_PA_C1   Study_Age_90
8072113-1-0238039048    M       1       62
7137284-1-0238040982    F       -9      89
4810747-1-0224008214    F       -9      62
3986861-1-0238042063    M       1       51
9435642-1-0238096441    M       2       61
6676256-1-0238040704    F       -9      90



/home/guosa/hpc/epimarker/pdf
for i in `ls *pdf`
do
time(convert $i $i.tiff)
done 


 
plink --bfile ../S_Hebbring_Unr.Guo.Forward --make-bed --keep RA.SamList.txt --out PMRP.PhaseII.Steven.Guo.RA

/gpfs/home/guosa/hpc/db/hg38

CACTCCCAACCCCTTT

AAAGGGGGTTGGGSGTG

>X58075.1 Human LINE-1 transposon (L1Hs) DNA
GGGGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGACGCAGAAG
ACGGGTGATTTCTGCATTTCCATCTGAGGTACCGGGTTCATCTCACTAGGGAGTGCCAGACAGTGGGCGC
AGGCCACTGTGTGCGCGCACCGTGCGCGAGCCGAAGCAGGGCGAGGCATTGCCTCACCTGGGAAGCGCAA
GGGGTCAGGGAGTTCCCTTTCCGAGTCAAAGAAAGGGGTGACGGACGCACCTGGAAAATCGGGTCACTCC
CACCCGAATATTGCGCTTTTCAGACCGGCTTAAGAAACGGCGCACCACGAGACTATATCCCACACCTGGC
TCAGAGGGTCCTACGCCCACGGAATCTCGCTGATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCG
GCAACGAGGCTGGGGGAGGGGCGCCCGCCATTGCCCAGGCTTGCTTAGGTAAACAAAGCAGCCGGGAAGC
TCGAACTGGGTGGAGCCCACCACAGCTCAAGGAGGCCTACCTGCCTCTGTAGGCTCCACCTCTGGGGGCA
GGGCACAGACAAACAAAAAGACAGCAGTAACCTCTGCAGACTTAAGTGTCCCTGTCTGACAGCTTTGAAG
AGAGCAGTGGTTCTCCCAGCACGCAGCTGGAGATCTGAGAACGGGCAGACTGCCTCCTCAAGTGGGTCCC
TGACCCCTGACCCCCGAGCAGCCTAACTGGGAGGCACCCCCCAGCAGGGCACACTGACACCTCACACGGC
AGGGTATTCCAACAGACCTGCAGCTGAGGGTCCTGTCTGTTAGAAGGAAAACTAACAACCAGAAAGGACA
TCTACACGAAAACCCATCTGTACATCACCATCATCAAAGACCAAAAGTAGATAAAACCACAAAGATGGGG
AAAAAACAGAACAGAAAAACTGGAAACTCTAAAACGCAGAGGCCCTCTCCTCCTCCAAAGGAACGCAGTT
CCTCACCAGCAACAGAACAAAGCTGGATGGAGAATGATTTTGACGAGCTGAGAGAAGAAGGCTTCAGACG
ATCAAATTACTCTGAGCTACAGGAGGACATTCAAACCAAAGGCAAAGAAGTTGAAAACTTTGAAAAAAAT
TTAGAAGAATGTATAACTAGAATAACC

/gpfs/home/guosa/hpc/db/miRNASNP/Gwas.Catalog.immune.disease.GRCH38.bed 

bedtools intersect -wao -a miRNA.SNP150.hg38.bed.bed -b Gwas.Catalog.immune.disease.GRCH38.bed | awk '(match($10,"chr"))'


How to calculate cost for WGBS:
http://epicore.med.cornell.edu/pricelist.php
Price per sample for WGBS library prep (internal core client): $250.00
Price per lane for PE150 sequencing (internal core client): $2650.00
Cost of library prep for 8 samples: $2000.00
Cost of sequencing 2 PE150 lanes: $5300.00
Total cost: $7300.00
Cost per sample: $912.50



for i in chr{1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  > $i.job >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedgraph >> $i.job
qsub $i.job
done

commonSNP150.hg38

awk '(match($4,/rs/)){print $4}' commonSNP150.hg38


awk '$3< -0.4 & M1<0.3 & M2>0.6' chr*.tab.matrix.rltSig.pvalue.rlt | sort -u | awk -F':|-|\t' '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > ./pdf/region.bed
cd /home/guosa/hpc/epimarker/bedgraph/pdf
bedtools sort -i region.bed > region.bed.sort
bedtools merge -i region.bed.sort > region.bed
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' region.bed > region.cor.bed 
awk '{print $3"\t"$5-5000"\t"$6+1000"\t"$13}' ../../../db/hg38/refGene | sort -u | awk '$2>0' > RefGene.hg38.5k.bed
bedtools intersect -wao -a region.cor.bed -b RefGene.hg38.5k.bed | awk '{print $8"\t"$4}' | sort -u 
bedtools intersect -wao -a region.cor.bed -b RefGene.hg38.5k.bed | awk '{print $8"\t"$4}' | sort -u | awk '$1 != "."' > download.list.txt
awk '{print $1}' download.list.txt | wc -l
cd ../



| awk '{print $1}' | sort -u | wc -l

bedtools intersect -wao -a region.bed -b RefGene.hg38.5k.bed | awk '{print $8"\t"$4}' | sort -u | awk '$1 != "."' | awk '{print $2}' | sort -u | wc -l


cd /gpfs/home/guosa/run/bedmethyl
for i in chr{7..22} chrX chrY chr6 chr5 chr4 chr3 chr2 chr1 
do
for j in `ls *bw`
do
echo \#PBS -N $j.$i  > $j.$i.job
echo \#PBS -l nodes=1:ppn=1 >> $j.$i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $j.$i.job
echo bigWigAverageOverBed $j /gpfs/home/guosa/hpc/db/hg38/window200/hg38.$i.win2K.bed $j.$i.tab >> $j.$i.job
qsub $j.$i.job
done
done


for i in chr{10..22} chrX chrY chrM
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org
echo \#PBS -m abe  > $i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $i.job
echo perl ~/hpc/bin/tab2matrix.pl $i \> $i.tab.matrix.rlt >> $i.job
echo Rscript --vanilla dmr.R $i.tab.matrix.rlt >> $i.job
qsub $i.job
done


for i in chr{1..22} chrX chrY chrM
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org
echo \#PBS -m abe  > $i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $i.job
echo Rscript --vanilla ../dmr.R $i.tab.matrix.rlt >> $i.job
qsub $i.job
done




for i in chr{10..12} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
echo cd /home/guosa/hpc/epimarker/bedgraph >> $i.job
echo perl ~/hpc/bin/tab2matrix.pl $i \> $i.tab.matrix.rlt >> $i.job
echo Rscript --vanilla dmr.R $i.tab.matrix.rlt >> $i.job
sh $i.job &
done


cd /home/guosa/hpc/db/hg19
mkdir window2000/
for i in {1..22} X Y M
do
perl ~/hpc/bin/cutchrosome.pl chr$i 500 >  ./window2000/hg38.chr$i.win2K.bed
done



for i in `ls *.fa`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/db/hg38/fa/chroms >> $i.job
echo perl cgpositionFinder.pl $i >> $i.job
qsub $i.job
done

/gpfs/home/guosa/hpc/db/hg38/fa/chroms/chr1.CpG.positions.txt
/gpfs/home/guosa/hpc/db/hg38/blueprint

ENCFF023MDK   100 scale


scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata3/*sorted.bam* ./
scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata/*.bam* ./
Pool_C-tissue_S6_L002_R1_001.bam
Pool_A-tissue_S4_L002_R1_001.bam
Pool_B-tissue_S5_L002_R1_001.bam
Pool_A-plasma_S1_L001_R1_001.bam
Pool_C-plasma_S3_L001_R1_001re.bam
Pool_B-plasma_S2_L001_R1_001re.bam



for i in `ls *.bedgraph`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedmethyl >> $i.job
echo wigToBigWig $i ~/hpc/db/hg38/hg38.chrom.sizes $i.bw >> $i.job
qsub $i.job
done

/gpfs/home/guosa/hpc/db/hg38/hg38.chrom.sizes

sh ENCFF103DNU.bed.bedgraph.job


/gpfs/home/guosa/run/bedmethyl

for i in `ls *.bed`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $i.job
echo perl bedMethyl2bedgraph.pl $i >> $i.job
echo wigToBigWig $i.bedgraph ~/hpc/db/hg38/hg38.chrom.sizes $i.bedgraph.bw >> $i.job
qsub $i.job
done

cd /gpfs/home/guosa/run/bedmethyl
for i in chr{7..22} chrX chrY chr6 chr5 chr4 chr3 chr2 chr1 
do
for j in `ls *bw`
do
echo \#PBS -N $j.$i  > $j.$i.job
echo \#PBS -l nodes=1:ppn=1 >> $j.$i.job
echo cd /gpfs/home/guosa/run/bedmethyl >> $j.$i.job
echo bigWigAverageOverBed $j /gpfs/home/guosa/hpc/db/hg38/window200/hg38.$i.win2K.bed $j.$i.tab >> $j.$i.job
qsub $j.$i.job
done
done

use strict;
use Cwd;
chdir getcwd;


for i in {1..22} X Y M
do
perl chrosomeCut.pl chr$i 2000 >>  hg38.cut2K.bed
done 

for i in {1..22} X Y M
do
perl chrosomeCut.pl chr$i 500 >  ./window200/hg38.chr$i.win2K.bed
done

/gpfs/home/guosa/hpc/db/hg38/hg38.win2K.bed



my $CHR=shift @ARGV;
my $STEP=shift @ARGV;
my $chrLen="/gpfs/home/guosa/hpc/db/hg38/hg38.chrom.sizes";
open F,$chrLen || die "cannot open $chrLen\n";
my %len;
while(<F>){
chomp;
my($chr,$len)=split/\s+/;
$len{$chr}=$len;
}
my $start=1;
my $end=$len{$CHR};
while($start<($len{$CHR}-$STEP)){
$end=$start+2000;
my $id="$CHR:$start-$end";
print "$CHR\t$start\t$end\t$id\n";
$start=$start+$STEP+1;
}


 
/gpfs/home/guosa/hpc/db/hg38/hg38.cut2K.bed


TDIR="/gpfs/home/guosa/hpc/epimarker/GbedGraph";
for i in `ls *.bigWig *.bw`
do
echo \#PBS -N $i  > $TDIR/$i.job
echo \#PBS -l nodes=1:ppn=1 >> $TDIR/$i.job
echo bigWigAverageOverBed /gpfs/home/guosa/hpc/epimarker/wig/$i /gpfs/home/guosa/hpc/db/hg38/window2000cut/hg38.cut2K.bed $TDIR/$i.tab >> $TDIR/$i.job
done

for(i in )




http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/bone_marrow/S01GRF/Acute_Lymphocytic_Leukemia/Bisulfite-Seq/CNAG/S01GRFA1.CPG_methylation_calls.bs_call.GRCh38.20160531.bw
http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/bone_marrow/15548/Multiple_Myeloma/Bisulfite-Seq/CNAG/S00XDKU1.CPG_methylation_calls.bs_call.GRCh38.20160531.bw

wget -e robots=off -nH -r -nd --reject="index.html*" --cut-dirs=10 --reject-regex 'GRCh37|RNA-Seq|ChIP-Seq' --accept-regex 'Bisulfite-Seq'  -l 6 -A "*CPG_methylation_calls.bs_call.GRCh38.20160531.bw" --no-parent http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/

wget http://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Adult_Liver/BI.Adult_Liver.Bisulfite-Seq.3.wig.gz

cd /gpfs/home/guosa/run/roadmap
8260aae218fcfbe0d8075c383dd86fb5  BI.Adult_Liver.Bisulfite-Seq.3.wig.gz
8260aae218fcfbe0d8075c383dd86fb5  BI.Adult_Liver.Bisulfite-Seq.3.wig.gz


mv /gpfs/home/guosa/hpc/db/hg38/roadmap/marker/include_sample.txt /gpfs/home/guosa/hpc/db/epimarker

include_sample.txt

/gpfs/home/guosa/hpc/db/hg38
/gpfs/home/guosa/hpc/db/epimarker/QC/CpGI


wget -e robots=off -l 6 --no-parent -R "index.htm*" -A "*_call.GRCh38.20160531*.bw" -nH -r -nd http://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38

scl enable devtoolset-2 bash

use strict;
use strict;
use Cwd;
use POSIX;
my $dir = getcwd;

my @file=glob("*.bigWig.bw");
foreach my $file(@file){
    my ($fn)=split/.bw/,$file;
	print "$fn\n";
}



use strict;
use Cwd;
use POSIX;
my $dir = getcwd;

my @file=glob("*.bigWig.1");
foreach my $file(@file){
    my ($fn)=split/.1/,$file;
	system("mv $file $fn")
}

for i in `ls *bw`
do
md5sum $i >>md5sum.roadmap.txt
done


wigToBigWig Bcell_PB_I44_fracmeth.wig /gpfs/home/guosa/hpc/db/hg38/hg38.chrom.sizes Bcell_PB_I44_fracmeth.bw





cd /gpfs/home/guosa/run
for i in `ls *wig.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run/roadmap >> $i.job
echo gunzip $i >> $i.job
echo wigToBigWig $i.w /gpfs/home/guosa/hpc/db/hg38/hg38.chrom.sizes $i.bw  >> $i.job
echo rm $i >> $i.job
qsub $i.job
done

wget http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O /gpfs/home/guosa/hpc/db/hg38/hg38.chrom.sizes


bigWigToBedGraph 65158.CEEHRC.CEMT0026.gDNA.WGB-Seq.methylation_profile.bigWig 65158.CEEHRC.CEMT0026.gDNA.WGB-Seq.methylation_profile.bedgraph
bigWigToBedGraph 63942.KNIH.KNIH_001_genomic_DNA.WGB-Seq.methylation_profile.bigWig  63942.KNIH.KNIH_001_genomic_DNA.WGB-Seq.methylation_profile.bedgraph


bigWigToBedGraph 24895.Blueprint.ERS433789.WGB-Seq.signal.bigWig 24895.Blueprint.ERS433789.WGB-Seq.signal.bedgraph



bigWigAverageOverBed  24895.Blueprint.ERS433789.WGB-Seq.signal.bigWig CpGI.test.bed 24895.Blueprint.ERS433789.WGB-Seq.signal.CpGI
bigWigAverageOverBed  24896.Blueprint.ERS433789.WGB-Seq.signal.bigWig CpGI.test.bed 24896.Blueprint.ERS433789.WGB-Seq.signal.CpGI


bigWigAverageOverBed  24833.Blueprint.ERS222241.WGB-Seq.signal.bigWig CpGI.test.bed 24833.Blueprint.ERS222241.WGB-Seq.signal.CpGI
bigWigAverageOverBed  24834.Blueprint.ERS222241.WGB-Seq.signal.bigWig CpGI.test.bed 24834.Blueprint.ERS222241.WGB-Seq.signal.CpGI

paste 24833.Blueprint.ERS222241.WGB-Seq.signal.CpGI 24834.Blueprint.ERS222241.WGB-Seq.signal.CpGI > test2.txt


chr7	128939699	128941569	chr7:128939699-128941569

for i in `ls *.bw`
do
bigWigAverageOverBed $i IRF5.bed $i.IRF5
done

/gpfs/home/guosa/hpc/epimarker/bedmethyl/*bed


sam<-read.table("/gpfs/home/guosa/run/cell.index.txt",head=F,sep="\t")
data=read.table(args[1],head=T,sep="\t",row.names=1,check.names=F)
filename=unlist(strsplit(colnames(data),split=".IRF5"))
newsam=sam[match(filename,sam[,1]),]
blood=which(newsam[,4]=="Blood")
solid=which(newsam[,4]=="Solid")
rlt<-TtestPValue(data,blood,solid)
output=paste(args[1],".pvalue.rlt",sep="")
colnames(rlt)=c("Pvalue","delta","M1","M2","SD1","SD2")


wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O ../hpc/db/hg19/hg19.chrom.sizes

wget http://epigenomesportal.ca/tracks/KNIH/hg19/63942.KNIH.KNIH_001_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63946.KNIH.KNIH_002_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63950.KNIH.KNIH_006_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63954.KNIH.KNIH_003_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63958.KNIH.KNIH_004_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63962.KNIH.KNIH_005_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63966.KNIH.KNIH_007_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63970.KNIH.KNIH_008_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63974.KNIH.KNIH_009_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63978.KNIH.KNIH_010_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/63982.KNIH.KNIH_011_genomic_DNA.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64277.KNIH.CKD23_C_Mesan_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64282.KNIH.CKD24_C_Podo_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64287.KNIH.CKD25_C_Podo_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64292.KNIH.CKD27_C_Mesan_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64297.KNIH.DB31_N_Alpha_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64302.KNIH.IPS01_N_Fibroblast_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64307.KNIH.IPS02_N_NPC_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64312.KNIH.IPS03_N_ENeuron_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64317.KNIH.IPS04_X_Fibroblast_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64322.KNIH.IPS05_X_NPC_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64327.KNIH.IPS06_X_ENeuron_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64332.KNIH.OB56_N_PreA_WGBS.WGB-Seq.methylation_profile.bigWig
wget http://epigenomesportal.ca/tracks/KNIH/hg19/64337.KNIH.OB57_D_PreA_WGBS.WGB-Seq.methylation_profile.bigWig

cd /gpfs/home/guosa/run
for i in `ls *DEEP*.bigWig *CEEHRC*.bigWig`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run >> $i.job
echo bigWigToWig $i $i.wig  >> $i.job
echo awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4/100}\' $i.wig \> $i.w >> $i.job
echo wigToBigWig $i.w /gpfs/home/guosa/hpc/db/hg19/hg19.chrom.sizes $i.bw  >> $i.job
echo rm $i >> $i.job
qsub $i.job
done

for i in `ls *KNIH*.bigWig`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run >> $i.job
echo bigWigToWig $i $i.wig  >> $i.job
echo awk \'{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4/10}\' $i.wig \> $i.w >> $i.job
echo wigToBigWig $i.w /gpfs/home/guosa/hpc/db/hg19/hg19.chrom.sizes $i.bw  >> $i.job
echo rm $i >> $i.job
qsub $i.job
done

for i in chr{1..22} chrX chrY
do
for j in `ls *.bigWig`
do
echo \#PBS -N $j.$i  > ./tab/$j.$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./tab/$j.$i.job
echo cd /gpfs/home/guosa/run/ >> ./tab/$j.$i.job
echo bigWigAverageOverBed $j ./tab/$i.window.bed ./tab/$j.$i.tab >> ./tab/$j.$i.job
done
done

for i in `ls *.chr22.job`
do
qsub $i
done

mkdir IRF5
for i in `ls *.bigWig`
do
bigWigAverageOverBed $i IRF5.bed ./IRF5/$i.IRF5
done

bigWigToBedGraph 65158.CEEHRC.CEMT0026.gDNA.WGB-Seq.methylation_profile.bigWig 65158.CEEHRC.CEMT0026.gDNA.WGB-Seq.methylation_profile.bedgraph
bigWigToBedGraph 63982.KNIH.KNIH_011_genomic_DNA.WGB-Seq.methylation_profile.bigWig 63982.KNIH.KNIH_011_genomic_DNA.WGB-Seq.methylation_profile.bedgraph
bigWigToBedGraph 24881.Blueprint.ERS337125.WGB-Seq.signal.bigWig 24881.Blueprint.ERS337125.WGB-Seq.signal.bedgraph


bigWigCorrelate 24871.Blueprint.ERS337105.WGB-Seq.signal.bigWig 24872.Blueprint.ERS337105.WGB-Seq.signal.bigWig
# R2=0.93
bigWigCorrelate  24885.Blueprint.ERS433781.WGB-Seq.signal.bigWig 24886.Blueprint.ERS433781.WGB-Seq.signal.bigWig
# R2=0.90

bedtools intersect -wao -a test.bed -b chr7.tab.txt.pvalue.rlt.ref.sort > chr7.tab.txt.pvalue.rlt.ref.sort.bed

DEEP and KNIH should /100

63297.DEEP.51_Hf01_BlCM_Ct.WGB-Seq.methylation_profile.bigWig.chr7.tab

bigWigToBedGraph -chrom=chr7 -start=128939735 -end=128941534  ENCFF087KYO.bigWig  ENCFF087KYO.bigWig.bedgraph
bigWigToBedGraph -chrom=chr7 -start=128939735 -end=128941534  C005VG51.CPG_methylation_calls.bs_call.GRCh38.20160531.bw  C005VG51.CPG_methylation_calls.bs_call.GRCh38.20160531.bw.bedgraph




data.frame(data[match("chr7:128579146-128581146",rownames(data)),], newsam[,4])
sam<-read.table("/gpfs/home/guosa/run/cell.index.txt",head=F,sep="\t")
filename=unlist(strsplit(colnames(data),split=paste(".",chr,".tab",sep="")))
newsam=sam[match(filename,sam[,1]),]
rlt<-data.frame(IRF5=t(data[match("chr7:128579146-128581146",rownames(data)),]),celltype=newsam[,1],celltype2=newsam[,2],Samtype=newsam[,4])
write.table(rlt,file="IRF5.txt",sep="\t",quote=F,col.names=NA,row.names=T)

for i in {1..22} X Y
do
for j in `ls *.$i.job`
do
echo $j
done
done

unlist(strsplit(colnames(data),split=paste(".",chr,".tab",sep="")))


for i in chr{1..22} chrX chrY
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run/tab >> $i.job
echo perl bigWigAverageOverBed2Matrix.pl $i  >> $i.job
qsub $i.job
echo $i.job
done


for i in chr{1..22} chrX chrY
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/run/tab >> $i.job
echo Rscript --vanilla dmr.R $i.tab.txt >> $i.job
qsub $i.job
echo $i.job
done

cp /gpfs/home/guosa/run/tab

for i in chr{1..22} chrX chrY
do
perl windowcut.pl $i 500 > $i.window.bed &
done

qsub 63950.KNIH.KNIH_006_genomic_DNA.WGB-Seq.methylation_profile.bigWig.2.job

 scp *.pl nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/bigwig

I just come back from a seminar in our department about HLA-DRB*58:01 and severe cutaneous adverse reactions (SCAR) caused by Allopurinol. I found the first paper was published in PNAS and we can talk about some Pharmacogenetic 


roadmap /home/guosa/hpc/db/hg19/wgbs

awk '{print $1"\t"$2"\t"$3"\t"$4}' GWAS-SNP-List_hg19_full.txt | sort -u > GWAS-SNP-List_hg19_full.txt.uni
bedtools window -w 500 -c  -a GWAS-SNP-List_hg19_full.txt.uni -b GWAS-SNP-List_hg19_full.txt.uni > GWAS-SNP-List_hg19_full.txt.uni.bedgraph

for i in `ls *bigWig`
do
bigWigAverageOverBed $i test.bed $i.IRF5.tab
echo $i
done



wget -nd -r -c --tries=20 ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/Bisulfite-Seq/ -o log

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes


bigWigAverageOverBed 65787.CEEHRC.CEMT0087.gDNA.WGB-Seq.methylation_profile.bigWig test.bed 65787.CEEHRC.CEMT0087.gDNA.WGB-Seq.methylation_profile.bigWig.IRF5.tab


Omics association study mediated by epigenetic emissary CpG-SNP 


bigWig



for i in chr{1..22} chrX chrY
do
cp $i.fa commonCpGSNP &
done


for i in chr{1..22} chrX chrY
do
plink --bfile $i --recode --tab --out $i &
done

cd /gpfs/home/guosa/hpc/db/1000Genome


/home/local/MFLDCLIN/guosa/hpc/db/hg19/GTEx_Analysis_v7_eQTL

vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &



cat *cpgsnp.bed >> hg19.CpG-SNP.bed
awk '{print $1,"\t",$2,"\t",$3,"\n"}' hg19.CpG-SNP.bed > Hg19.CpG-SNP.bed 
mv Hg19.CpG-SNP.bed.2 Hg19.CpG-SNP.bed

bedtools window -w 500 -c  -a Hg19.CpG-SNP.bed -b Hg19.CpG-SNP.bed > Hg19.CpG-SNP.bedgraph


cd /home/guosa/hpc/db/hg19/plan2/7mer
GTACGCA.positions.bed
GATCGCA.positions.bed
/home/guosa/hpc/db/hg19/wgbs

for i in A T C G
do
perl ~/hpc/bin/GenomeReferenceNucleotideSummary.pl ~/hpc/db/hg38/hg38.fa $i >> SingleNucleotideSummary.txt
done

for i in A T C G
do
for j in A T C G
do
perl ~/hpc/bin/GenomeReferenceNucleotideSummary.pl ~/hpc/db/hg38/hg38.fa $i$j >> diNucleotideSummary.txt
done
done 



/home/guosa/hpc/db/hg19/plan3

perl ~/hpc/bin/GenomeReferenceNucleotideSummary.pl ~/hpc/db/hg19/fa/chr1.fa A
perl ~/hpc/bin/GenomeReferenceNucleotideSummary.pl ~/hpc/db/hg38/hg38.fa A



wget ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
wget ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz  
wget ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi
wget ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz
wget ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz.tbi


mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz  ALL.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi ALL.chrX.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
mv ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz ALL.chrY.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
mv ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz.tbi ALL.chrY.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz
ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz


BIRCDEV13-LC


cg00000108      37459205        37459207        37459206        37459256                        rs9857774       1       LC      .;.     48      1       0       XY_NO   A_NO    37458757        37458758        449     C3orf35 CCDS46792

cg00000108      cg00000108      12709357        ATACAATAAAACAAACCTAAAATAATCCTAACTCCRCTATCATCCTAACC                      II                      TCCATTTTGAAGGAAAAAAATGAAGGCTCTGAAAGTGTAAATCGCTTACTGAAGGGCACA[CG]GCCAGGATGACAGCGGAGCCAGGATCACC


cd /home/local/MFLDCLIN/guosa/hpc/db/hg19/plan2
for i in `ls *bed.sort`
do 
awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $i > $i.bed
done

ssh nu_guos@transfer.chtc.wisc.edu 

wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi
  

for i in chr{1..22} chrX chrY
do
perl dinucleotideFinder.pl ../fa/$i.fa CG
done

for i in chr{1..22} chrX chrY
do
perl ../dinucleotideFinder.pl $i CG  &
done
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.G.positions.txt --counts --out chr22.G &

for i in chr{1..22} chrX chrY
do
perl ../dinucleotideFinder.pl chr22 $i  &
done


for i in chr{1..22} chrX chrY
do
for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.$i.*.genotypes.vcf.gz --recode --positions chr$i.$j.vcf.C.positions.txt --counts --out chr22.$j.1 &
done
done


#!/usr/bin/sh

chr="chr1"
for i in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
perl ../dinucleotideFinder.pl $chr $i  &
done

for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
wc -l $chr.$j.positions.bed
done


for i in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions $chr.$i.vcf.C.positions.txt --counts --out $chr.$i.1 &
done

for i in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions $chr.$i.vcf.G.positions.txt --counts --out $chr.$i.2 &
done

for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
wc -l $chr.$j.1.frq.count
done

for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
wc -l $chr.$j.2.frq.count
done

for chr in chr{1..22} chrX chrY
do
for i in GTACGCA
do
perl 7merMotifFinder.pl $chr $i  &
done
done
cat *GTACGCA.positions.bed >> GTACGCA.positions.txt
bedtools sort -i GTACGCA.positions.txt > GTACGCA.positions.bed

for chr in chr{1..22} chrX chrY
do
for i in GATCGCA
do
perl 7merMotifFinder.pl $chr $i  &
done
done
cat *GATCGCA.positions.bed >> GATCGCA.positions.txt
bedtools sort -i GATCGCA.positions.txt > GATCGCA.positions.bed

/home/guosa/hpc/db/hg19/plan2/7mer







for i in chr{} chrX chrY
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.G.positions.txt --counts --out chr22.G &


/home/guosa/hpc/db/1000Genome/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
vcftools --gzvcf /home/guosa/hpc/db/1000Genome --positions chr22.AT.positions.txt --recode  --counts --out chr22.AT.vcf 

vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions test.pos --recode --counts --out chr22

vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.G.positions.txt --counts --out chr22.G &


cd /gpfs/home/guosa/hpc/db/1000Genome
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --recode --snp CpGSNP.list.hg19.txt --derived --out CpG-SNP-MAF


cd /gpfs/home/guosa/hpc/db/1000Genome
vcftools --gzvcf /home/local/MFLDCLIN/guosa/hpc/db/1000Genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --recode --snps CpGSNP.list.hg19.txt --derived --out CpG-SNP-MAF
vcftools --gzvcf /home/local/MFLDCLIN/guosa/hpc/db/1000Genome/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --recode --snp rs782048968 --derived --out CpG-SNP-MAF




for i in `ls *bed.sort.bed.uni.bed`
do 
sort -u $i >> allCpG-SNP.hg19.bed
done

bedtools sort -i allCpG-SNP.hg19.bed > allCpG-SNP.hg19.bed.sort
bedtools cluster -d 1000 -i  allCpG-SNP.hg19.bed.sort >  allCpG-SNP.hg19.G1000.bed.sort.bed
bedtools cluster -d 500 -i  allCpG-SNP.hg19.bed.sort >  allCpG-SNP.hg19.G500.bed.sort.bed
bedtools cluster -d 250 -i  allCpG-SNP.hg19.bed.sort >  allCpG-SNP.hg19.G250.bed.sort.bed &

Rscript --vanilla summary.R allCpG-SNP.hg19.G1000.bed.sort.bed.summary
Rscript --vanilla summary.R allCpG-SNP.hg19.G500.bed.sort.bed.summary
Rscript --vanilla summary.R allCpG-SNP.hg19.G250.bed.sort.bed.summary


* Web server with PHP 7.0.0 or HHVM 3.18.5 or higher.
* A SQL server, the following types are supported
** MySQL 5.5.8 or higher
** PostgreSQL 9.2 or higher
** SQLite 3.3.7 or higher
** Oracle 9.0.1 or higher
** Microsoft SQL Server 2005 (9.00.1399)


http://www.10.103.135.149/index.php

sudo mount -t cifs //luministuscdata.file.core.windows.net/luministphase1data ./luministuscdata -o vers=3.0,username=luministuscdata,password=02grOPN59qcoAQ2EVaSad1z/28YkDe7j1SA6woz36VnbIdXkfhn8tf40JB+WPRvrUnNGWo7SQRTLaJRGANjH1Q==,dir_mode=0777,file_mode=0777,sec=ntlmssp

sudo mount -t cifs //luministuscdata.file.core.windows.net/sequencingdatarun2 ./luministuscdata2 -o vers=3.0,username=luministuscdata,password=02grOPN59qcoAQ2EVaSad1z/28YkDe7j1SA6woz36VnbIdXkfhn8tf40JB+WPRvrUnNGWo7SQRTLaJRGANjH1Q==,dir_mode=0777,file_mode=0777,sec=ntlmssp

sudo mount -t cifs //luministuscdata.file.core.windows.net/sequencingdatarun3 ./luministuscdata3 -o vers=3.0,username=luministuscdata,password=02grOPN59qcoAQ2EVaSad1z/28YkDe7j1SA6woz36VnbIdXkfhn8tf40JB+WPRvrUnNGWo7SQRTLaJRGANjH1Q==,dir_mode=0777,file_mode=0777,sec=ntlmssp


scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata3/*sorted.bam* ./
scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata/*.bam* ./

Pool_C-tissue_S6_L002_R1_001.bam
Pool_A-tissue_S4_L002_R1_001.bam
Pool_B-tissue_S5_L002_R1_001.bam
Pool_A-plasma_S1_L001_R1_001.bam
Pool_C-plasma_S3_L001_R1_001re.bam
Pool_B-plasma_S2_L001_R1_001re.bam


cutadapt -o ./test/Pool_B-tissue_S5_L002_R1_001.fastq.gz -p ./test/Pool_B-tissue_S5_L002_R2_001.fastq.gz  Pool_B-tissue_S5_L002_R1_001.fastq.gz Pool_B-tissue_S5_L002_R2_001.fastq.gz

cutadapt -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC Pool_B-tissue_S5_L002_R1_001.fastq.gz

cutadapt -o ./test/Pool_B-tissue_S5_L002_R1_001.fastq.gz -p ./test/Pool_B-tissue_S5_L002_R2_001.fastq.gz  Pool_B-tissue_S5_L002_R1_001.fastq.gz Pool_B-tissue_S5_L002_R2_001.fastq.gz

perl smartbismark.pl input2.txt submit
qsub Pool_B-tissue_S5_L002_R1_001.fastq.gz.job

 
for i in `ls Pool*job`
do
qsub $i
done

scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/data2/input2.txt ./


-rwxrwxrwx 1 root root  50G Jun  6 20:22 Pool_C-plasma_S3_L001_R2_001.fastq.gz
-rwxrwxrwx 1 root root  53G Jun  6 20:22 Pool_B-plasma_S2_L001_R2_001.fastq.gz
-rwxrwxrwx 1 root root  46G Jun  6 20:22 Pool_B-plasma_S2_L001_R1_001.fastq.gz
-rwxrwxrwx 1 root root  44G Jun  6 20:22 Pool_A-plasma_S1_L001_R2_001.fastq.gz
-rwxrwxrwx 1 root root  44G Jun  6 20:22 Pool_C-plasma_S3_L001_R1_001.fastq.gz
-rwxrwxrwx 1 root root  39G Jun  6 21:29 Pool_A-plasma_S1_L001_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    40G Jun 15 18:44 Pool_A-tissue_S4_L002_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    44G Jun 15 18:44 Pool_A-tissue_S4_L002_R2_001.fastq.gz
-rwxrwxrwx  1 root   root    42G Jun 15 18:44 Pool_B-tissue_S5_L002_R2_001.fastq.gz
-rwxrwxrwx  1 root   root    38G Jun 15 18:44 Pool_B-tissue_S5_L002_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    48G Jun 15 18:44 Pool_C-tissue_S6_L002_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    53G Jun 16 02:51 Pool_C-tissue_S6_L002_R2_001.fastq.gz

hpc/tools/FastQC

ssh 'guosa@mfldclin.org'@10.103.160.225
ssh -L 22:23.99.137.107:22 -p 22 nu_guos:submit-3.chtc.wisc.edu

ssh nu_guos@submit-3.chtc.wisc.edu 'shg047@23.99.137.107 "cp /home/shg047/data1/fastq/input.txt"' > input.txt
scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu 1 nc %h %p' user@remote2:path/to/file .

scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/data1/fastq/* ./
scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/data2/*gz ./
scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/bin/* ./

scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/bin/* ./

fastqc_v0.11.7.zip

 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/bowtie2-2.3.4.1-linux-x86_64.zip
 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/*.zip
 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/*.bz2
 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/*.gz


scp .ssh/id_rsa.pub -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:~/.ssh/authorized_keys
scp .ssh/id_rsa.pub -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:~/.ssh/authorized_keys2

ssh-keygen -t rsa
ssh nu_guos@submit-1.chtc.wisc.edu 'mkdir -p .ssh'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-1.chtc.wisc.edu 'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-1.chtc.wisc.edu 'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-1.chtc.wisc.edu 'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-1.chtc.wisc.edu 'chmod 640 ~/.ssh/authorized_keys2'
ssh nu_guos@submit-1.chtc.wisc.edu
alias chtc="ssh nu_guos@submit-3.chtc.wisc.edu"


###########################
cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo vcftools --gzvcf chr$i.vcf.gz --exclude mydup.txt --out chr$i.unique >> chr$i.job
qsub chr$i.job
done
#########################

ssh-keygen -t rsa
ssh shg047@23.99.137.107 mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 640 ~/.ssh/authorized_keys2'
ssh shg047@23.99.137.107 
alias azure="shg047@23.99.137.107"

MiSeq, ~25 million reads,  1.5G, 88Gene, 250K, 0.5G, 3 Sample, 65 hours,750$/3=250$/sample
HiSeq 2500, ~400 million reads, 6 G, 88Gene, 250K, 0.5G, 12 Sample, 24 hours,300$/12=25$/sample
HiSeq 4000, ~ 5 billion reads,  1.5T, 250K, 0.5G, 3000 Sample, 84 hours,20000$/3000=12$/sample
HiSeq X Ten, ~ 6 billion reads, 1.8T, 250K, 0.5G, 3600 Sample, 72 hours,12000/3600=3$/sample

https://designstudio-array.illumina.com/#/upload-targets
plink --bfile plink --extract AllCpGSNP150.hg19.RS_SNP.uni.txt --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input

for i in `ls *cpgsnp.bed`
do
bedtools sort -i $i > $i.sort
done

for i in `ls *cpgsnp.bed.sort`
do
awk '{print $1"\t",$}'
done


awk '{print $1"\t",$2,"\t",$3,"\t",$7}' Hg19.CpG-SNP.bed > Hg19.CpG-SNP.bedgraph

for i in chr{1..22} chrX chrY
do
cat $i.hg19_cpgsnp.bed.sort >> AllCpGSNP150.hg19.bed
done


rm AllCpGSNP150.hg19.RS_SNP.txt
rm AllCpGSNP150.hg19.RS_SNP.uni.txt
for i in chr{1..22} chrX chrY
do
awk '{print $7}' $i.hg19_cpgsnp.bed.sort >> AllCpGSNP150.hg19.RS_SNP.txt
done
sort -u AllCpGSNP150.hg19.RS_SNP.txt > AllCpGSNP150.hg19.RS_SNP.uni.txt


/home/local/MFLDCLIN/guosa/hpc/pmrp/phase2/RA


a7:b2:da:ad:ca:09:25:e7:80:8a:59:46:90:31:b1:b6 MFLDCLIN\guosa@birc7-lc

ssh-keygen -t rsa
ssh nu_guos@submit-3.chtc.wisc.edu mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'chmod 640 ~/.ssh/authorized_keys2'
ssh nu_guos@submit-3.chtc.wisc.edu
alias chtc="ssh nu_guos@submit-3.chtc.wisc.edu"


ssh-keygen -t rsa
ssh nu_guos@128.105.244.191 'mkdir -p .ssh'
cat .ssh/id_rsa.pub | ssh nu_guos@128.105.244.191 'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh nu_guos@128.105.244.191 'chmod 700 ~/.ssh'
#cat .ssh/id_rsa.pub | ssh nu_guos@128.105.244.191 'chmod 640 ~/.ssh/authorized_keys2'
#cat .ssh/id_rsa.pub | ssh nu_guos@128.105.244.191 'cat >> .ssh/authorized_keys2'
ssh nu_guos@128.105.244.191
alias chtc="ssh nu_guos@submit-1.chtc.wisc.edu"





ssh-keygen -t rsa
ssh shg047@23.99.137.107 mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 640 ~/.ssh/authorized_keys2'
ssh shg047@23.99.137.107 
alias azure="shg047@23.99.137.107"

bgzip	-c	S_Hebbring_Unr.Guo.vcf	>	S_Hebbring_Unr.Guo.vcf.gz
tabix	-p	vcf	S_Hebbring_Unr.Guo.vcf.gz

vcftools	--gzvcf	S_Hebbring_Unr.Guo.vcf.gz	--remove-indels	--recode	--recode-INFO-all	--out	S_Hebbring_Unr.Guo.reindel
bgzip	-c	S_Hebbring_Unr.Guo.reindel.vcf	>	S_Hebbring_Unr.Guo.reindel.vcf.gz
tabix	-p	vcf	S_Hebbring_Unr.Guo.reindel.vcf.gz


plink	--bfile	..	S_Hebbring_Unr.Guo	--keep	parmloss.txt	--chr	23	--allow-no-sex	--recode	--tab	--transpose	--out	6055529-1-0224008846

setwd("C:\\Users\\guosa\\Downloads")
data<-read.table("1176608-1-0238062177.tped",head=F)
het<-apply(data,1,function(x)	sum(!	as.character(x[5])==as.character(x[6])))
plot(het~data$V4,col="red",cex=2,xlab="Chromosome	X",ylab="Heterozygote")



/gpfs/home/guosa/hpc/hemochromatosis/haplotype

/home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype/exomechip_SNV_PASS_BEAGLE_chr6_phased_sel2.map


sel.map


methylation related CpG-SNPs in PMRP Project

for i in `ls *map.map`
do
awk '{print "chr"$2,$5,$7,$6}' $i > $i.newbed
done


HLA-B chr6:31323299-31324,734
MICA chr6:31367561-31383090
NOTCH4 chr6:32163725-32165371


NOTCH4 chr6:31323299-32165371
28493789-33425320

use lib qw(~/hpc/tools/Cwd-Ext-1.06/modulos/share/perl5); # You may need to change this path

/gpfs/home/guosa/hpc/tools/Statistics-R-0.02/modulos/share/perl5

Statistics/R.pm

perl Makefile.PL PREFIX=./modulos
make
make install


for i in chr{1..22} X Y
do
awk '{$1==chr$i}' ../snp150.hg19.txt >>chr$i.vcf.bed
done


for i in chr{1..22} chrX chrY chrM
do
awk -v chr="$i" '$1==chr' ../snp150.ls .txt >> $i.vcf.bed
echo $i
done

for i in chrM
do
awk -v chr="$i" '$1==chr' ../snp150.hg19.txt >> $i.vcf.bed
echo $i
done



'A/G' => 'R',
'C/T' => 'Y',
'A/C' => 'M',
'G/T' => 'K',
'C/G' => 'S',
'A/T' => 'W',
'A/C/T' => 'H',
'C/G/T' => 'B',
'A/C/G' => 'V',
'A/G/T' => 'D',
'A/C/G/T' => 'N',

chr17   80679961        80679962
chr17   80679983        80679984
chr17   80679984        80679985
chr17   80679989        80679990
chr17   80679990        80679991
chr17   80680036        80680037

chr19   93503   93504   rs2353749       -       A       G/T


chr19   1939700 1939701
grep 11138 ../snp150.hg19.txt | grep chr18
grep 93503 ../snp150.hg19.txt | grep chr19
method 1: chr19   668570/4804043
chr19   119331  119332




# http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/index.html
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/HumanOmni2.5-4v1_B-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/HumanOmni25-8v1-2_A1-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/humanomniexpress-24v1-0_a-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/HumanOmniZhongHua-8v1-2_A-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/Immuno_BeadChip_11419691_B-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/InfiniumCoreExome-24v1-1_A-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/InfiniumExome-24v1-1_A1-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/InfiniumPsychArray-24v1-1_A2-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/Consortium-OncoArray_15047405_A-b37.Source.strand.zip

for i in `ls *.strand`
do
awk '{print "chr"$2,"\t",$3-1,"\t",$3,"\t","chr"$2":"$3}' $i > $i.bed
done

for i in `ls *.strand.bed`
do
bedtools sort -i $i > $i.sort
bedtools intersect -wa -a $i.sort -b ../allSNP150.hg19.cpg-allele-v4.txt > $i.sort.CpGSNP.bed
wc -l $i.sort.CpGSNP.bed
done

90545

for i in `ls *.strand.bed`
do
wc -l $i 
done



for i in `ls *fwd.txt`
do
awk '{print $1}' $i >> Hapmap2
done
sort -u Hapmap2 > Hapmap3
mv 

bedtools intersect -wao -a /home/local/MFLDCLIN/guosa/hpc/db/Hapmap/phase2-3/oncotarget.bed -b hapmapPhaseIIISummary.txt.cpgsnp

bedtools intersect -wao -a /home/local/MFLDCLIN/guosa/hpc/db/Hapmap/phase2-3/oncotarget.bed -b allSNP150.hg19.cpg-allele-v4.txt


rs767677267                 

chr1	114446340	114451307	PTPN22
 
plink --bfile chr1 --recode 12 fastphase --make-bed --chr 1 --extract PTPN22.txt --keep CEU.txt --from-bp 114446340 --to-bp 114451307 --out PTPN22


Potenital methylation Loading(PML)

scp shg047@23.99.137.107:/home/shg047/data2/A09_S* ./


scp shg047@23.99.137.107:/home/shg047/data1/*pl ./

scp shg047@23.99.137.107:/home/shg047/pkg/*zip ./
scp shg047@23.99.137.107:/home/shg047/pkg/*bz2 ./

should_transfer_files = YES
transfer_input_files = us.dat, wi.dat
when_to_transfer_output = ON_EXIT
log = job.log
output = job.out
error = job.err
request_cpus = 1
request_memory = 20MB
request_disk = 20MB




/home/shg047/pkg


set C/G as reference allele

plink --vcf PTPN22.CPG.vcf.gz --reference-allele mylist.txt --recode vcf --real-ref-alleles



https://imputation.sanger.ac.uk/
eagle --vcf PTPN22.vcf --geneticMapFile=/home/local/MFLDCLIN/guosa/hpc/tools/Eagle_v2.4/tables/genetic_map_hg19_withX.txt.gz --outPrefix PTPN22.CPG

https://data.broadinstitute.org/alkesgroup/Eagle/downloads/
	

fastPHASE -F 5000 PTPN22.chr-1.recode.phase.inp
bedtools intersect -wao -a target.bed -b commonSNP_hg19_v1_sort_cluster.bed.summary
bedtools intersect -wao -a target2.bed -b commonSNP_hg19_v1_sort.txt 
 
rs10858022
rs1217397
rs971173
rs28381068
rs3761936
rs114661042
rs1217390
ssh nu_guos@submit-3.chtc.wisc.edu
 
ssh nu_guos@submit-3.chtc.wisc.educd 
/home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype

chr12:4,505,480-4,592,607
plink --file exomechip_SNV_PASS_BEAGLE_chr8_phased_sel2 --make-bed --out exomechip_SNV_PASS_BEAGLE_chr8_phased_sel2

bedtools sort -i commonSNP_hg19_v1.txt > commonSNP_hg19_v1_sort.txt
bedtools cluster -d 1600 -i commonSNP_hg19_v1_sort.txt > commonSNP_hg19_v1_sort_cluster.bed


grep rs12201499 allsnp150.hg19 > test.txt
grep rs6913437 allsnp150.hg19 >> test.txt
grep rs11966502 allsnp150.hg19 >> test.txt
grep rs73716691 allsnp150.hg19 >> test.txt
grep rs757671 allsnp150.hg19 >> test.txt

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/recombRate.txt.gz

awk '{print $2,$3,$4,$5,$7,$8,$10}' hg19.commonsnp150 > hg19_commonsnp150_trim.txt


12269 
gcc --version
libcurl4-openssl-dev
apt-get install libcurl4-openssl-dev
sudo apt-get install texlive
sudo apt-get install texlive-fonts-extra

source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library("GEOquery")

gse <- getGEO("GSE50579", GSEMatrix = TRUE)
show(gse)
filePaths = getGEOSuppFiles("GSE21653")
filePaths
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])
df1 <- getGSEDataTables("GSE3494")
lapply(df1, head)

cd /home/local/MFLDCLIN/guosa/hpc/pmrp/merge
# make the first time merge to find out allele need to be filped. 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --bmerge S_Hebbring_Unr  --out PMRP-Phase1-phase2-Full
# then filp anyone dataset
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP-Phase1-phase2-Full.missnp --make-bed --out 
# then filp anyone dataset and then merge again. 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip --bmerge S_Hebbring_Unr --out PMRP-Phase1-phase2-Full
# This the remaining non-merged alleles should be indels. run indel2indel.pl to change phase 2 indel mode to phase I. 
perl indel2indel.pl > S_Hebbring_Unr.bim.bim
mv S_Hebbring_Unr.bim.bim S_Hebbring_Unr.bim
# merge again for the last time 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --bmerge S_Hebbring_Unr  --out PMRP-Phase1-phase2-Full
# Now you get the merge dataset


plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip --extract PMRP-Phase1-phase2-Full.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip_INDEL
plink --bfile S_Hebbring_Unr --extract PMRP-Phase1-phase2-Full.missnp --make-bed --out S_Hebbring_Unr_INDEL

--out PMRP-Phase1-phase2-Full

cp /home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bed ./
cp /home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim ./
cp /home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.fam ./

cp /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr.fam ./
cp /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr.bim ./
cp /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr.bed ./

cp /home/guosa/hpc/pmrp/phase2/S_Hebbring_Unr.Guo.Forward.bed ./
cp /home/guosa/hpc/pmrp/phase2/S_Hebbring_Unr.Guo.Forward.bim ./
cp /home/guosa/hpc/pmrp/phase2/S_Hebbring_Unr.Guo.Forward.fam ./

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield  --list-duplicate-vars --out FinalRelease_QC_20140311_Team1_Marshfield
plink --bfile S_Hebbring_Unr.Guo.Forward  --list-duplicate-vars --out S_Hebbring_Unr.Guo.Forward

perl exm2rs.pl FinalRelease_QC_20140311_Team1_Marshfield.bim 
perl exm2rs.pl  

perl bim2indel.pl FinalRelease_QC_20140311_Team1_Marshfield.bim.bim > FinalRelease_QC_20140311_Team1_Marshfield.bim
perl bim2indel.pl S_Hebbring_Unr.Guo.Forward.bim.bim > S_Hebbring_Unr.Guo.Forward.bim

plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip

plink-merge.missnp

grep Multiple plink.log | awk -F\' '{print $2}' > M.txt


plink --make-bed --vcf FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --out FinalRelease_QC_20140311_Team1_Marshfield
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --remove exclude.txt --out FinalRelease_QC_20140311_Team1_Marshfield_Clean
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --list-duplicate-vars
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --flip PMRP-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Clean.flip

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --chr 1 --from-bp 9804693 --to-bp 9819993 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out PMRP_2018_Phase_1_and_Phase_2.Guo

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip --exclude PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip_demiss
plink --bfile  S_Hebbring_Unr.Guo.Forward --exclude PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out  S_Hebbring_Unr.Guo.Forward_demiss

plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out PMRP_2018_Phase_1_and_Phase_2.Guo

1161 genomic postion annotation mistakes in Phase 2 dataset and I replace them with the phase I annotation position. 358556 exm id were changed to rs id. 



plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out Guo

perl bim2indel.pl FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim > FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim.bim
mv FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim.bim FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out Guo

cp S_Hebbring_Unr.Guo.Forward S_Hebbring_Unr.Guo.Forward.bim
perl bim2indel.pl S_Hebbring_Unr.Guo.Forward.bim > S_Hebbring_Unr.Guo.Forward.bim.bim
mv S_Hebbring_Unr.Guo.Forward.bim.bim S_Hebbring_Unr.Guo.Forward.bim

cp FinalRelease_QC_20140311_Team1_Marshfield_Flip FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim
perl exm2rs.pl FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim 
mv FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim.bim FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim


1, first solve the indel problem in phase 2
2, change exm id to rs id in phase 1
3, unify phase 1 and phase 2 indel genotype
4, 

# phase I data
perl exm2rs.pl FinalRelease_QC_20140311_Team1_Marshfield.bim 
perl bim2indel.pl FinalRelease_QC_20140311_Team1_Marshfield.bim.bim > FinalRelease_QC_20140311_Team1_Marshfield.bim
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink
grep Multiple plink.log | awk -F\' '{print $2}' > M.txt
perl ReplaceBimPosition.pl S_Hebbring_Unr.Guo.Forward.bim M.txt FinalRelease_QC_20140311_Team1_Marshfield.bim > S_Hebbring_Unr.Guo.Forward.bim.bim
mv S_Hebbring_Unr.Guo.Forward.bim.bim S_Hebbring_Unr.Guo.Forward.bim
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip plink-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip




plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out Guo




data<-read.table("FinalRelease_QC_20140311_Team1_Marshfield.bim.bim")

for(i in 1:nrow(data)){
if(data[i,2])
}
grep "1KG"   | awk '{ print $2 }' |  xargs grep {} S_Hebbring_Unr.Guo.Forward.bim

awk 'NR==FNR{a[$1];next;}$1 in a' FinalRelease_QC_20140311_Team1_Marshfield.bim S_Hebbring_Unr.Guo.Forward.bim

CoreExome_24v1p2_A1_Anno.Guo.2018.March.short.csv


plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip.fam --extract rs7550295.txt --recode --tab --out test
plink --bfile S_Hebbring_Unr.Guo.Forward --extract rs7550295.txt --recode --tab --out test
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --merge-equal-pos --make-bed --out PMRP_2018_Phase_1_and_Phase_2.Guo

FinalRelease_QC_20140311_Team1_Marshfield.bim

/home/guosa/hpc/pmrp/phase1
check raw proble id
/home/guosa/hpc/pmrp/phase1
/home/guosa/hpc/pmrp/phase1
plink --make-bed --vcf FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --out FinalRelease_QC_20140311_Team1_Marshfield

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --list-duplicate-vars

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --flip PMRP-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Clean.flip

phas1<-read.table("/home/guosa/hpc/pmrp/phase1/myplink.bim")
phas2<-read.table("/home/guosa/hpc/pmrp/phase2/S_Hebbring_Rel.Guo.bim")
head(phas1)
head(phas2)



plink --bfile data1 --bmerge data2.ped data2.map --make-bed --out merge
plink --bfile data1 --bmerge data2.ped data2.map --make-bed --out merge



fam<-read.table("FinalRelease_QC_20140311_Team1_Marshfield_Clean.fam",sep="")
saminfo<-read.table("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",head=T,sep="\t")
head(fam)
head(saminfo)
fam[,5]=saminfo[match(fam[,1],saminfo[,1]),9]
write.table(fam,file="FinalRelease_QC_20140311_Team1_Marshfield_Clean.fam2",sep=" ",row.names=F,col.names=F,quote=F)


# phase I data

cd /home/local/MFLDCLIN/guosa/hpc/pmrp/phase1
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield
plink2 --bfile FinalRelease_QC_20140311_Team1_Marshfield --pca approx  --maf 0.05 --memory 40000 --mds-plot --out phase1.pca

--mds-plot

setwd("/home/local/MFLDCLIN/guosa/hpc/pmrp/phase1")
eigenvec<-read.table("phase1.pca.eigenvec",head=F)
saminfo<-read.table("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",head=T,sep="\t")
sam<-saminfo[match(as.character(eigenvec[,1]),saminfo$Sample.Name),]

pdf("phase1.pca1-population.pdf")
Legends<-unique(data.frame(Population=sam$PrimaryAncestry_Clean,Col=as.numeric(sam$PrimaryAncestry_Clean)))
plot(eigenvec[,4]~eigenvec[,3],cex=0.55,col=as.numeric(sam$PrimaryAncestry_Clean),pch=as.numeric(sam$PrimaryAncestry_Clean),xlab="PC1",ylab="PC2")
legend("topright",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()

pdf("phase1.pca2-population.pdf")
Legends<-unique(data.frame(Population=sam$PrimaryAncestry_Clean,Col=as.numeric(sam$PrimaryAncestry_Clean)))
plot(eigenvec[,6]~eigenvec[,5],cex=0.55,col=as.numeric(sam$PrimaryAncestry_Clean),pch=as.numeric(sam$PrimaryAncestry_Clean),xlab="PC3",ylab="PC4")
legend("bottomright",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()

pdf("phase1.pca9-sex.pdf")
Legends<-unique(data.frame(Population=sam$PrimaryAncestry_Clean,Col=as.numeric(sam$PrimaryAncestry_Clean)))
plot(eigenvec[,12]~eigenvec[,11],cex=0.55,col=as.numeric(sam$PrimaryAncestry_Clean),pch=as.numeric(sam$PrimaryAncestry_Clean),xlab="PC3",ylab="PC4")
legend("bottomright",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()


data<-read.table("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",head=T,sep="\t")
rlt<-apply(data[,23:ncol(data)],2,function(x) summary(glm(data$Sex~x))$coefficients[2,4])
rlt<-apply(data[,23:ncol(data)],2,function(x) summary(glm(data$Age_Clean~x))$coefficients[2,4])
rlt<-apply(data[,23:ncol(data)],2,function(x) summary(glm(as.numeric(data$PrimaryAncestry_Clean)~x))$coefficients[2,4])


plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --remove exclude.txt --out FinalRelease_QC_20140311_Team1_Marshfield_Clean




# check PC4< -0.1

Phase I	Eigen-value
PC1	16.6006
PC2	11.9216
PC3	8.30266
PC4	6.48869
PC5	6.23224
PC6	6.07716
PC7	5.91141
PC8	5.66893
PC9	5.40722
PC10 5.22692

HAMP
https://www.ncbi.nlm.nih.gov/nuccore/NC_000019.10?report=genbank&from=35282346&to=35285143
MALSSQIWAACLLLLLLLASLTSGSVFPQQTGQLAELQPQDRAGARASWMPMFQRRRRRDTHFPICIFCCGCCHRSKCGMCCKT
http://useast.ensembl.org/Homo_sapiens/Gene/Compara_Tree?db=core;g=ENSG00000105697;r=19:35280716-35285143;time=1530041560143.143


hg19 genomic region for HLA-B and MICA
chr6:31,232,075-31,391,038

awk '$7==CEU'


cd /home/local/MFLDCLIN/guosa/hpc/db/1000Genome

awk '$7=="CEU"' integrated_call_samples_v2.20130502.ALL.ped > CEU.txt
plink --bfile chr6 --make-bed --keep CEU.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CEU.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()


awk '$7=="CHB"' integrated_call_samples_v2.20130502.ALL.ped > CHINA.txt
awk '$7=="CHS"' integrated_call_samples_v2.20130502.ALL.ped >> CHINA.txt
plink --bfile chr6 --make-bed --keep CHINA.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CHINA.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()

awk '$7=="JPT"' integrated_call_samples_v2.20130502.ALL.ped > JPT.txt
plink --bfile chr6 --make-bed --keep JPT.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-JPT.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()


chr6:31,980,794-33,074,227



awk '$7=="CHB"' integrated_call_samples_v2.20130502.ALL.ped > CHINA.txt
awk '$7=="CHS"' integrated_call_samples_v2.20130502.ALL.ped >> CHINA.txt
plink --bfile chr6 --make-bed --keep CHINA.txt --maf 0.05 --snps-only --chr 6 --from-bp 31980794 --to-bp 33074227
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out NOTCH-HLADPQ.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("green", "red"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CHINA-NOTCH.pdf")
levelplot(input, main="LD between NOTCH4 and HLA-D-P-Q in CHINESE", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()







plink --file fA --merge-list allfiles.txt --make-bed --out mynewdata

10124 Pmrp.Steven.Study.Individual.List.txt


grep rs4475691 /home/local/MFLDCLIN/guosa/hpc/db/hg19/allsnp150.hg19
grep rs4475691 /home/local/MFLDCLIN/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim

/home/local/MFLDCLIN/guosa/hpc/pmrp/phase1

my $bim="/home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim";
my $allsnp="/home/guosa/hpc/db/hg19/allsnp150.hg19";

https://watson.hgen.pitt.edu/pub/dbvor/dbvor.1.11.tgz
https://watson.hgen.pitt.edu/pub/slink/apps
https://watson.hgen.pitt.edu/pub/haplo/haplo.HP.Z


6       exm535823       0       32546879        G       A
6       exm535982_ver2  0       32552096        A       G
6       exm2122129_ver4 0       32552121        T       A

for i in `ls *fa`
do
perl cgpositionFinder.pl $i
done


cp /home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype/*newbed

cd /home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype/

plink --file exomechip_SNV_PASS_BEAGLE_chr6_phased_sel2 --make-bed --snps-only --chr 6 --from-bp 32544312 --to-bp 32552716 --out DRB1

awk '{print $2}' plink.bim > mysnps.txt


/home/local/MFLDCLIN/guosa/hpc/db/hg19/fa/chr10.CpG.positions.txt

awk '{print $2,$3,$4,$5}' hg19.commonsnp150 | grep chr10 > chr10.hg19.snp.txt

grep rs35499948 allsnp150.hg19

chr10.hg19.snp.txt



# The easiest way to get lubridate is to install the whole tidyverse:
install.packages("tidyverse")
# Alternatively, install just lubridate:
install.packages("lubridate")
# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("tidyverse/lubridate")


lSLF.AD.Z

# The easiest way to get lubridate is to install the whole tidyverse:
install.packages("tidyverse")
# Alternatively, install just lubridate:
install.packages("lubridate")
# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("tidyverse/lubridate")


install.packages("caret")
install.packages("mlbench")
install.packages("devtools")

library("devtools")
install_github("tidyverse/lubridate")

library(mlbench)
library(caret)
library(ggpubr)
library(car)
library("BBmisc")
shapiro.test(data$GCC.FA.Z)
shapiro.test(powerTransform(data$GCC.FA.Z))
newdata <- preProcess(data[,3:ncol(data)], method=c("BoxCox"))
tidyverse/lubridate

input1<-data$GCC.FA.Z
input2<-powerTransform(data$GCC.FA.Z)
pdf("GCC.FA.Z.pdf")
par(mfrow=c(2,2))
ggqqplot(input1)
ggqqplot(input2)
dev.off()


wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

phase3.pca.eigenvec

plink --bfile ../S_Hebbring_Unr.Guo --keep parmloss.txt --chr 23 --allow-no-sex --recode --tab --out 1176608-1-0238062177

plink --bfile S_Hebbring_Unr.Guo --snp rs16964862 --allow-no-sex --recode --tab --out rs16964862
plink --bfile S_Hebbring_Unr.Guo.Forward --snp rs16964862 --allow-no-sex --recode --tab --out rs16964862.Forward
plink --bfile S_Hebbring_Unr.Guo --snp exm184366 --allow-no-sex --recode --tab --out exm184366




rs16964862
plink --bfile exomechip_SNV_PASS_BEAGLE_chr13_phased_sel2 --snp rs16964862 --allow-no-sex --recode --tab --out rs16964862
cd /home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype
plink --bfile S_Hebbring_Unr --update-alleles top_to_AB.txt –-make-bed –-out S_Hebbring_Unr.Forward
rs16964862

# hair color
plink --bfile  S_Hebbring_Unr.Guo.Forward --snp rs117322171 --allow-no-sex --recode --tab --out rs117322171
/home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype

grep T,G 

ln -s /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data  /home/local/MFLDCLIN/guosa/hpc/pmrp/phase2/rawdata/S_Hebbring_2128_Released_Data

awk '{pirnt $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$9}' allsnp150

/mqnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/Strand_Translation_Files

# raw plink file of autism 
/mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data

# copy to my folder
cd /home/local/MFLDCLIN/guosa/hpc/autism/data
cp /mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data/All_samples_Exome_QC.bed ./
cp /mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data/All_samples_Exome_QC.fam ./
cp /mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data/All_samples_Exome_QC.bim ./

# change phen xlsx to plink phenotype file
install.packages("openxlsx")

# 314 samples listed in quantitative iffusion tensor MRI data file
# 256 samples have whole-exom sequencing data
# 14 samples don't have MRI quantitative measurement 
# 19 samples were removed since quality control
# 242 samples were included for the assciation study (MRI ~ Allele + age)

# outlier samples
1, 2 low call rate ( u70704cl, u69388s)
2, Gender discrepancy(u38386cl,u65210cl,317814-UW)
3.1, Family(336051,395993 remove, keep 372278)
3.2, Family(370121 remove, keep 386915)
3.3  MZ twin(u28908s remove, keep u28906s-B-Redo)
3.4  duplicated samples (Saliva vs cell line) (u68413d remove, keep u68413s)
4.0  9 PCA outlier(u62997s,u1941001s,u90503s,u59502cl,u64061s,u65457s,u810031s,u810030s,u59504s)
totally, these samples were removed (space sparate): u62997s,u1941001s,u90503s,u59502cl,u64061s,u65457s,u810031s,u810030s,u59504s,u28908s,u68413s,370121,336051,395993,u38386cl,u65210cl,317814-UW,u70704cl,u69388s

Actually, I found the data have already removed 4.0, 3.4, 3.3, 


library("openxlsx")
phen=read.xlsx("DougsPhenotypes_with_avg_age_QTL.xlsx", sheet = 1)
fam<-read.table("All_samples_Exome_QC.fam")
rank1<-match(unlist(lapply(fam[,1],function(x) gsub("u","",strsplit(as.character(x),"-|s|c")[[1]][1]))),phen$ID)
newphen<-cbind(fam,phen[rank1,])
colnames(newphen)[1:2]<-c("FID","IID")
colnames(newphen)[11]<-c("AvgAge")

# prepare covariates plink file (All_samples_Exome_QC.cov, only have age)
cov<-newphen[,c(1,2,11)]
cov[is.na(cov)]<- -9
write.table(cov,file="DougsPhenotypes_with_avg_age_QTL.newphen.cov",sep="\t",quote=F,col.names=c("FID","IID","AvgAge"),row.names=F)

# prepare multiple phenotype plink file (All_samples_Exome_QC.phen)
mphen<-newphen[,c(1,2,13:ncol(newphen))]
mphen[is.na(mphen)]<- -9
write.table(mphen,file="All_samples_Exome_QC.phen",sep="\t",quote=F,col.names=F,row.names=T)

# use perl script to submit job and creat result file with corresponding phenotype names
exc<-read.table("excludeSample.txt")
match(exc[,1],newphen[,1])

# 1031667 variants removed due to MAF<0.01 and 1618874 variants and 249 samples pass filters and QC.  
plink2 --bfile All_samples_Exome_QC --ci 0.95 --genotypic --maf 0.01 --remove excludeSample.txt --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov --linear --out test
plink --bfile All_samples_Exome_QC --ci 0.95 --linear mperm=500  --maf 0.01 --remove excludeSample.txt --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov  --out test
plink --bfile binary_fileset --recode vcf-iid --out new_vcf
plink --bfile All_samples_Exome_QC --ci 0.95 --linear perm --aperm 10 1000000 0.0001 0.01 5 0.001  --maf 0.01 --remove excludeSample.txt --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov  --out test  


plink --bfile All_samples_Exome_QC --snp exm184366 --allow-no-sex --recode --tab --out exm184366
plink --bfile All_samples_Exome_QC --snp kgp8664031 --allow-no-sex --recode --tab --out kgp8664031
plink --bfile All_samples_Exome_QC --snp rs9786510 --allow-no-sex --recode --tab --out rs9786510

plink --bfile All_samples_Exome_QC --snp rs2761764 --allow-no-sex --recode --tab --out rs2761764
plink --bfile All_samples_Exome_QC --snp rs17118281 --allow-no-sex --recode --tab --out rs17118281





awk '$12<0.00005 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt | grep rs6500552 

awk '$12<0.000000031 {print FILENAME,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' *.assoc.linear | grep AD



# find permuation significant assocaiton
0.05/1618874=3.1*10^-8

awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm 


awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt

awk '$3<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt | grep rs6500552 


kgp8664031

plink2 --bfile All_samples_Exome_QC --impute-sex
 --impute-sex
 
cov<-newphen[,c(1,2,11)]
cov[is.na(cov)]<- -9
mphen<-newphen[,c(1,2,10,13:ncol(newphen))]
mphen[is.na(mphen[,3]),3]<- -9
mphen[mphen[,3]=="Case",3]=1
mphen[mphen[,3]=="Control",3]=0
mphen[is.na(mphen)]<- -9
write.table(mphen,file="All_samples_Exome_QC.phen",sep="\t",quote=F,col.names=T,row.names=F)

# PCA on chrX data
 
7065331-1-0238095238

# collect chrX haplotype and remove missing 
plink --bfile S_Hebbring_Unr.Guo --chr 23 --allow-no-sex --noweb --recode --tab --out chrX
plink --bfile S_Hebbring_Unr.Guo --chr 24 --allow-no-sex --noweb --recode --tab --out chrY

plink --bfile chrXY --pca approx  --maf 0.05 --memory 30000 --out chrXY.pca

plink2 --file chrX --pca approx --maf 0.05 --memory 30000 --out chrX.pca



# IDn plink 23: X, 24: Y, 25: X+Y, 26:M 
# plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 23 --allow-no-sex --noweb --recode --tab --out chrX
# plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 24 --allow-no-sex --noweb --recode --tab --out chrY
system("plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 23 --allow-no-sex --noweb --recode --tab --out chrX")
system("plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 24 --allow-no-sex --noweb --recode --tab --out chrY")
chrx<-read.table("chrX.ped",stringsAsFactors=F,colClasses = c("character"))
chrx<-chrx[,-which(apply(chrx,2,function(x) sum(x==0))>7000)]
chry<-read.table("chrY.ped",stringsAsFactors=F,colClasses = c("character"))
chry<-chry[,-which(apply(chry,2,function(x) sum(x==0))>7000)]
GenderFScore<-read.table("plink.sexcheck",head=T)
saminfo<-read.table("S_Hebbring_Release_Sample_Sheet.txt",head=T,sep="\t")
# chrx het
Tmp<-c()
for(i in 1:nrow(chrx)){
x<-as.character(unlist(chrx[i,seq(7,ncol(chrx),by=2)]))
y<-as.character(unlist(chrx[i,seq(8,ncol(chrx),by=2)]))
Tmp<-c(Tmp,sum(x==y))
print(i)
}
chrXHetRatio=1-Tmp/((ncol(chrx)-6)/2)
Tmp<-unlist(apply(chry[,7:ncol(chry)],1,function(x) sum(x=="0")))
chrYcallrate=1-(Tmp/(ncol(chry)-6))
pedigreeGender=saminfo[match(chrx[,1],saminfo$Sample_Name),]$Gender
Fscore<-GenderFScore[match(chrx[,1],GenderFScore$IID),]$F
rlt<-data.frame(IID=chrx[,1],chrXHetRatio,chrYcallrate,pedigreeGender,Fscore)
write.table(rlt,file="chrXY.gender.prediction.txt",sep="\t",quote=F)

# find permuation significant assocaiton
0.05/1618874=3.1*10^-8

awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm 
awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt
awk '$3<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt | grep rs6500552 

plink --file raw-GWA-data --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out raw-GWA-data
plink --bfile raw-GWA-data --extract raw-GWA-data.prune.in --genome --out raw-GWA-data
plink --file data --indep 50 5 2

plot(density(Tmp))

saminfo<-read.table("S_Hebbring_Release_Sample_Sheet.txt",head=T,sep="\t")
rlt<-saminfo[match(data[which(Tmp<100),2], saminfo$Sample_Name),]$Gender


mkdir data1
sudo mount -t cifs //discoverydata.file.core.windows.net/discoverydata ./data1 -o vers=3.0,username=discoverydata,password=jyZAwZ5yzqth7jwWSW+XOLCGB4Xsf4NvOwMNF8u7jwymXisTSO6MBT6YligXzmYIC/jsKuUiD84vRsYeoF5IQA==,dir_mode=0777,file_mode=0777,sec=ntlmssp
mkdir data2
sudo mount -t cifs //discoverydata.file.core.windows.net/discoverydata2 ./data2 -o vers=3.0,username=discoverydata,password=jyZAwZ5yzqth7jwWSW+XOLCGB4Xsf4NvOwMNF8u7jwymXisTSO6MBT6YligXzmYIC/jsKuUiD84vRsYeoF5IQA==,dir_mode=0777,file_mode=0777,sec=ntlmssp
mkdir data3
sudo mount -t cifs //discoverydata.file.core.windows.net/discoverydatatissue ./data3 -o vers=3.0,username=discoverydata,password=jyZAwZ5yzqth7jwWSW+XOLCGB4Xsf4NvOwMNF8u7jwymXisTSO6MBT6YligXzmYIC/jsKuUiD84vRsYeoF5IQA==,dir_mode=0777,file_mode=0777,sec=ntlmssp

sudo apt-get install libssl-dev
sudo apt-get install libxml2-dev

install.packages("openssl")
install.packages("git2r")
install.packages("httr")
install.packages("devtools") 
devtools::install_github("Azure/doAzureParallel") 



install pip 
install cutadapt 

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xzvf chromFa.tar.gz
for i in `ls *fa`
do
perl ../../bin/cgpositionFinder.pl $i &
done


for i in `ls *job`
do
sh $i &
done




tail -n 1 *.txt | awk '/^==>/ {a=substr($0, 5, length-8); next} {print a,$1}' | awk '$2>0 {if ($2 !~/chrY/) print $1}' | xargs -I {} rm {}
tail -n 1 *.hap | awk '/^==>/ {a=substr($0, 5, length-8); next} {print a,$1}' | awk '$2>0 {if ($2 !~/chrY/) print $1}' | xargs -I {} rm {}

for i in `ls *hapInfo.txt`
do
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/CpGI.hg19.bed4 --output $i.CpGI.hap
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/hg19.LINE.bed --output $i.LINE.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/NAS1/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed --output $i.mhb
done
perl ~/bin/hap2mhl.pl > mhl.txt

sudo rstudio-server stop
sudo rstudio-server start
sudo rstudio-server restart

.libPaths( c( .libPaths(), "/media/Home_Raid1/shg047/R/x86_64-pc-linux-gnu-library/3.4") )

# Estellar2016 claim GSM1279519 is a colon tissue, however, in my prediction it is Intestine, therefore, I compare his tissue with roadmap Small_Intestine and Colon
bigWigCorrelate GSM1120321_UCSD.Small_Intestine.Bisulfite-Seq.STL001.wig.gz.bw GSM1010989_UCSD.Sigmoid_Colon.Bisulfite-Seq.STL003.wig.gz.bw                  # 0.9788
bigWigCorrelate GSM1120321_UCSD.Small_Intestine.Bisulfite-Seq.STL001.wig.gz.bw ../../Estellar2016/bw/GSM1279519_CpGcontext.Colon.txt.bedgraph.sort.hg19.bw   # 0.6365
bigWigCorrelate GSM1010989_UCSD.Sigmoid_Colon.Bisulfite-Seq.STL003.wig.gz.bw ../../Estellar2016/bw/GSM1279519_CpGcontext.Colon.txt.bedgraph.sort.hg19.bw     # 0.637193

[[File:Enrichment-of-mhl-ccr.png|200px]]
[[File:Enrichment-mhl-LCP.png|200px]]

# rebuild Rstudio server in genome-miner Version 1.0.136 – © 2009-2016 RStudio, Inc.
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi rstudio-server-1.0.143-amd64.deb

.libPaths( c( .libPaths(), "/media/Home_Raid1/shg047/R/x86_64-pc-linux-gnu-library/3.4") )


Rstudio sever in genome-miner:
http://132.239.25.238:8787/

Rscript ~/bin/mhlpredict.R -f mhl.txt

for i in `ls *hapInfo.txt`
do
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/CpGI.hg19.bed4 --output $i.CpGI.hap
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/hg19.LINE.bed --output $i.LINE.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/NAS1/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed --output $i.mhb
done
perl ~/bin/hap2mhl.pl > mhl.txt

for i in `ls SRX*hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/CpGI.hg19.bed4 --output $i.CpGI.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/hg19.LINE.bed --output $i.LINE.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/NAS1/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed --output $i.mhb.hap
done

grep SINE rmsk.hg19.bed > hg19.SINE.bed
grep LINE rmsk.hg19.bed > hg19.LINE.bed
grep Simple_repeat rmsk.hg19.bed > hg19.LINE.bed

## How to Build your own R package
install.packages("devtools")     # for R package development
install.packages("roxygen2")     # for R package documents (R function => R man)
library("devtools")              # for R package development
library("roxygen2")              # for R package documents (R function => R man)
devtools::document()             # Creat tar.gz source package for distribution
devtools::build()                # Creat tar.gz source package for distribution
install.packages("../monod_1.3.tar.gz")
library("monod")
#' @export                       # before function to export to namespace
data.table::fread()              # Never use library() or require() in a R package!

install.packages("impute")

qsub SRR1035893.job # check what happened!

scp SRP* shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS1/Jenkinson2017NG/samplesheet/
du ./ -h --max-depth 1

cd /media/Home_Raid1/shg047/NAS1/Estellar2016/hapinfo
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/PRJNA229055.txt
for i in `ls SRX*`
do
perl ~/bin/hapinfo2wig.pl $i > $i.bedgraph
echo $i
done
perl ~/bin/

qstat -u shg047 | grep .read | awk '{print $1}' |xargs -I {} qdel {}
qstat -u shg047 | grep trimmed | awk '{print $1}' |xargs -I {} qdel {}
qstat -u shg047 | grep run1 | awk '{print $1}' |xargs -I {} qdel {}
qstat -u shg047 | grep run2 | awk '{print $1}' |xargs -I {} qdel {}

perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072071.txt
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072075.txt
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072078.txt
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072141.txt

wc -l ../samplesheet/SRP072071.txt
wc -l ../samplesheet/SRP072075.txt
wc -l ../samplesheet/SRP072078.txt
wc -l ../samplesheet/SRP072141.txt


grep SRX1651659 ../samplesheet/SRP072071.txt
grep SRX1651659 ../samplesheet/SRP072075.txt
grep SRX1651659 ../samplesheet/SRP072078.txt

scp /media/Home_Raid1/shg047/work/meth450/dyh/20150413/DataResults/1rawdata/idat/3999547166/*idat 


sguo@sph3736:~/Downloads/RA450$ scp shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/meth450/dyh/20150514/DataResults/1rawdata/idat/3999423021/*idat ./
+																																																																													
cp /media/Home_Raid1/dinh/RRBS_HaploInfo/Matrices/* /media/Home_Raid1/shg047/NAS1/mhl/

ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA,outlier.colour="white")+ coord_flip()


awk 'NR==7258037 {print}' 6-P-1.sorted.clipped.bam.bed

# define RRBS hotspot regions
cd /home/shg047/oasis/monod/bam
for i in `ls *bam`
do
echo $i 
bedtools bamtobed -i $i | awk '{print $1,$2,$3}' OFS="\t" | grep -v 'chrLambdaNEB' | sort -u | sort -k1,1 -k2,2n > $i.bed
bedtools merge -i $i.bed > $i.merge.bed
done
/gpfs/home/guosa/hpc/db/hg19/wgbs/slideWindow


awk '{if ($1>0.1) print $1,$2,$3}' OFS="\t" MM
awk 'i++ {if($1~/ARRS/) print i}' ../../bak/bak.db
awk 'END{if($1 !~/chrY/) print FILENAME,$i,$2,$3,$4}' 
awk 'END{if($1 !~/chrY/) print FILENAME,$i,$2,$3,$4}' *hapInfo.txt
awk -F\| '$3 > 0 { print substr($3,1,6)}' file1
tail -n 1 *.txt | awk '/^==>/ {a=substr($0, 5, length-8); next} {print a,$1}' | awk '$2>0 {if ($2 !~/chrY/) print $1}' | xargs -I {} rm {}


xargs -I {} qsub {}

cov2methcount.sh
#!/bin/bash
# Usage: sh bedgraph2amf.sh -b input.bedgraph -d input.bed 
# Extension: perl ~/bin/tab2matrix.pl > matrix.txt
while getopts i:o: option
do
 case "${option}"
 in
 b) cov=${OPTARG};;
 d) output=${OPTARG};;
 esac
done
awk '{print $1,$2,"+","CpG",$4/100,$5+$6}' $cov >  $output.methcount
pmd -o $output.pmd $output.methcount
hmr -o $output.hmr $output.methcount


cd /media/Home_Raid1/zhl002/NAS1/WGBS/analysis
bedgraph2amf bedgraph bed 
awk -F":|-" '{print ouch$1,$2,$3,$1":"$2"-"$3}' OFS="\t"

#!/bin/bash
# Usage: sh bedgraph2amf.sh -b input.bedgraph -d input.bed 
# Extension: perl ~/bin/tab2matrix.pl > matrix.txt
while getopts b:d: option
do
 case "${option}"
 in
 b) bedgraph=${OPTARG};;
 d) bed=${OPTARG};;
 esac
done
# help Alice to treat the methylation bedgraph data to AMF data
sort -k1,1 -k2,2n $bedgraph > $bedgraph.sort
awk '{print $1,$2,$3,$4}' OFS="\t" $bedgraph.sort > $bedgraph.bed4
awk '{print $1,$2,$3,$1":"$2"-"$3}' $bed > $bed.bed4
bedGraphToBigWig $bedgraph.bed4 ~/work/db/mm9/mm9.chrom.sizes $bedgraph.bw
bigWigAverageOverBed $bedgraph.bw $bed.bed4 $bedgraph.tab


cd /media/Home_Raid1/zhl002/NAS1/WGBS/analysis
for i in `ls *trim`
do
sort -k1,1 -k2,2n $i > $i.sort
awk '{print $1,$2,$3,$4}' OFS="\t" $i.sort > $i.bed4
bedGraphToBigWig $i.bed4 ~/work/db/mm9/mm9.chrom.sizes $i.bw
bigWigAverageOverBed $i.bw meth_expr_regions.bed.cor.bed $i.tab
perl ~/bin/matrixGeneration.pl > matrix.txt
print $i
done

find ./ -name "*.txt" -mtime +1 -type f | xargs -I {} scp shg047@genome-miner.ucsd.edu:/media/Home_Raid1/zhl002/NAS1/WGBS/analysis unsupervised \;

bedtools 
bedgraphtools
wigtools

methmatrix2gender.R

/home/shg047/oasis/db/FHM.bed


Genomic dataset for Normal Colon tissue
WGBS	SRX332737	SRR949213	Normal Colon
WGBS	SRX332737	SRR949214	Normal Colon
WGBS	SRX332737	SRR949215	Normal Colon
WGBS    SRX381553       GSM1279519      Colon_normal

Genomic dataset for Normal Lung tissue
WGBS      SRX1649893       SRX1649893         Normal Lung    
WGBS      SRX1651655       SRX1651655         Normal Lung
WGBS      SRX1651658       SRX1651658         Normal Lung
WGBS      GSM1279527       SRX381713          Normal Lung

WGBS	  STL001GA-01	   STL001GA-01	      Normal Lung
WGBS	  STL002LG-01	   STL002LG-01	      Normal Lung
WGBS	  N37-Lung	   N37-Lung	      Normal Lung
RRBS	  ENCFF000LVO	   ENCFF000LVO	      Normal Lung
RRBS	  ENCFF000LVR	   ENCFF000LVR	      Normal Lung

SRX1649893
SRX1651655
SRX1651658

awk '{if ($1>0.1) print $1,$2,$3}' OFS="\t" MM

/home/shg047/oasis/Ziller2013/sortbam/hapinfo
SRR949213 SRR949214 SRR949215
cd /media/NAS3_volume2/Dinh/WGBS_LTS33/Hg19/Ziller_Harvard/BAMfiles
scp Colon_Primary_Normal.chr*sorted.clipped.bam* shg047@tscc-login.sdsc.edu:/home/shg047/oasis/Ziller2013/sortbam/SRX332737
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed > saminfo.txt


for i in `ls *bed`; do awk -F"\t|\'|\/" '{print $1,$2,$3,$5/$6,$5,$6}' $i ; done

cd /home/shg047/oasis/Ziller2013/bw
for i in BiSeq_cpgMethylation_BioSam_1121_Colon_Adjacent_Normal.BiSeq.bed GSM1204465_BiSeq_cpgMethylation_BioSam_1120_Colon_Primary_Tumor.BiSeq.bed BiSeq_cpgMethylation_BioSam_157_REMC_19_colonic_mucosa.BiSeq.bed
do
#awk -F"\t|\'|\/" '{print $1,$2,$3,$5/$6}' OFS="\t" $i > $i.bedgraph
#bedGraphToBigWig $i.bedgraph ../../db/hg19/hg19.chrom.sizes $i.bw
bigWigAverageOverBed $i.bw /oasis/tscc/scratch/shg047/monod/hapinfo/Colon77CCP.txt $i.tab
done
perl ~/bin/tab2matrix.pl 



ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1279nnn/GSM1279527/suppl/GSM1279527_CpGcontext.Lung.txt.gz



DATA[match("chr9:103791378-103791447",rownames(DATA)),]
DATA[match("chr12:125051726-125051786",rownames(DATA)),]
DATA[match("chr5:174152050-174152076",rownames(DATA)),]
perl ~/bin/tab2matrix.pl | grep 'chr5:174152050-174152076'


wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1204nnn/GSM1204465/suppl/GSM1204465_BiSeq_cpgMethylation_BioSam_1120_Colon_Primary_Tumor.BiSeq.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1204nnn/GSM1204466/suppl/GSM1204466_BiSeq_cpgMethylation_BioSam_1121_Colon_Adjacent_Normal.BiSeq.bed.gz


cd /oasis/tscc/scratch/shg047/Ziller2013/fastq
perl ~/bin/smartbismark.pl --input saminfo.txt --genome hg19 --server TSCC --submit no --queue hotel
for i in SRR949210 SRR949211 SRR949212 SRR949213 SRR949214 SRR949215
do
qsub $i.pbs
done


/home/shg047/oasis/Chen2016CellResearch/bam
for i in SRR1654403 SRR1654398
do

echo " #!/bin/csh" > $i.job
echo " #PBS -N bam2sort" >> $i.job
echo " #PBS -q hotel" >> $i.job
echo " #PBS -l nodes=1:ppn=8" >> $i.job
echo " #PBS -l walltime=168:00:00" >> $i.job
echo " #PBS -V" >> $i.job
echo " #PBS -M shg047@ucsd.edu" >> $i.job
echo " #PBS -o $i.o" >> $i.job
echo " #PBS -e $i.e" >> $i.job
echo " #PBS -m abe" >> $i.job
echo " #PBS -A k4zhang-group" >> $i.job
echo "cd \$PBS_O_WORKDIR" >> $i.job
echo  samtools sort -@ 8 -o ../sortbam/$i\_bismark_bt2_pe.sort.bam ../bam/$i\_1_val_1_bismark_bt2_pe.nonCG_filtered.bam >> $i.job
echo  samtools index ../sortbam/$i\_bismark_bt2_pe.sort.bam >> $i.job
qsub $i.job
done


SRR1654396_bismark_bt2_pe.sort.bam


# How to add the rsa key of your own linux computer to remove linux server system to avoid input passwd everytime.
# first run the keygen command to generate the keys of your own computer. 
ssh-keygen -t rsa
# tell the remote servers the key of your own computer so that you don't need to input passwd everytime.
cat ~/.ssh/id_rsa.pub | ssh shg047@genome-miner.ucsd.edu 'cat >> ~/.ssh/authorized_keys'
cat ~/.ssh/id_rsa.pub | ssh shg047@tscc-login.sdsc.edu 'cat >> ~/.ssh/authorized_keys'

awk 'i++ {if($1~/ARRS/) print i}' ../../bak/bak.db
find ./ -name ".R" |xargs -I {} grep unsupervised \;
find ./ -name "*.R" | xargs grep 'unsupervised'
find ./ -name "*.R" | xargs grep 'combat'

for i in `ls *hapInfo.txt`
do
perl ~/bin/hapinfo2wig.pl $i > $i.wig
done

perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 


cd /media/Home_Raid1/shg047/work/monod/hapinfo
perl methHMH.pl 6-T-1.sorted.clipped.bam.hapInfo.txt 6-P-1.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-1 
perl methHMH.pl 6-T-2.sorted.clipped.bam.hapInfo.txt 6-P-2.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-2 
perl methHMH.pl 6-T-3.sorted.clipped.bam.hapInfo.txt 6-P-3.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-3 
perl methHMH.pl 6-T-4.sorted.clipped.bam.hapInfo.txt 6-P-4.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-4 
perl methHMH.pl 6-T-5.sorted.clipped.bam.hapInfo.txt 6-P-5.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-5 
cat caHMH* > HMH.txt
Rscript caHMH.R

#>>>>>>>>
caHMH.R
#>>>>>>>>
data<-read.table("HMH.txt",head=T,sep="\t",as.is=T)
input<-data.matrix(data[,8:ncol(data)])
caHMH<-data[which(apply(input,1,function(x) sum(x)==0)),1]
unique(as.character(caHMH))
write.table(unique(as.character(caHMH)),file="caHMH.rlt.txt",col.names=F,row.names=F,quote=F)

#>>>>>>>>
caHMH.sh
#>>>>>>>>
awk -F":|-" '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" caHMH.rlt.txt > caHMH.bed
bedtools intersect -wao -a caHMH.bed -b ~/work/db/hg19/hg19_refGene.bed | sort -k1,1n


cor2bed.sh
awk -F":|-" '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t"

chr19	58951755	58951921	chr19:58951755-58951921
scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt ~/work/db/hg19
awk '{print $1,$2,$2+1}' OFS="\t" ~/work/db/hg19/HsGenome19.CpG.positions.txt > HsGenome19.CpG.positions.bed

for i in `ls SRX*hapInfo.txt`
do
if [ -e "$i.hap" ]
then
perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 
fi
done

cp /home/shg047/oasis/Encode/hapinfo/ENC*hapInfo.txt ../../monod/hapinfo/
cp /home/shg047/oasis/Estellar2016/hapinfo/SRX*.hapInfo.txt ../../monod/hapinfo/

for i in `ls *hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 
done

for i in `ls N37*hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 
done

perl ~/bin/haptools.pl --input 6-T-1.sorted.clipped.bam.hapInfo.txt --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output test
head 6-T-1.sorted.clipped.bam.hapInfo.txt.hap
head test.hap


/media/12TB_ext/DD_Ext12T/RRBS_MONOD/Plasma_RRBS_2015/BAMfiles

addr:132.249.107.90    2017-05-07

#!/usr/bin/perl
use strict;
chdir "/home/shg047/work/monod/rrbs_kun";
my @file=glob("RRBS*P*");
foreach my $file(@file){
my (undef,$cancerType,$Samid,undef)=split /[-P]/,$file;
my $newName="$cancerType-P-$Samid";
system("cp $file $newName");
}

/media/12TB_ext/DD_Ext12T/Capture_MONOD/150209_SN216/SeqCap/BAMfiles


[59] "WB_centenarian.all_chrs"               
[60] "WB_middle-age.all_chrs"                
[61] "WB_new-born.all_chrs"                  
[62] "SRX381569_tumor_colon"                 
[63] "SRX381716_adenocarcinoma_lung"         
[64] "SRX381719_squamous_cell_tumor_lung"    
[65] "SRX381722_small_cell_tumor_lung" 

find . -name "*" -type d -exec rm -r "{}" \;

for i in `ls *bam.hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed ../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i
done

wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2600/E-MTAB-2600.sdrf.txt
perl -lane "{print @F[31]}" E-MTAB-2600.sdrf.txt | xargs -I {} wget {}

samtools tview ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam ~/oasis/db/hg19/hg19.fa chr10:100027918-100027944
samtools tview ../bam/6-P-1.sorted.clipped.bam ~/work/db/hg19/hg19.fa chr10:100027918-100027944

cp /media/12TB_ext/DD_Ext12T/RRBS_MONOD/MONOD_RRBS_2014_RemappedTogether_Jan2016/BAMfiles/BAMfiles/*bam ./
cp /media/12TB_ext/DD_Ext12T/RRBS_MONOD/Plasma_RRBS_2016/BAMfiles/*bam ./


perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed > saminfo.txt

perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt non bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt


for i in `ls *hapInfo.txt`
do
echo "perl methhaptools.pl --input $i --bed ../../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i"
done


cd /home/shg047/oasis/monod/hapinfo/bak20170501
for i in `ls *hapInfo.txt` 
do
echo " #!/bin/csh" > $i.job
echo " #PBS -N hapinfo2summary" >> $i.job
echo " #PBS -q hotel" >> $i.job
echo " #PBS -l nodes=1:ppn=1" >> $i.job
echo " #PBS -l walltime=72:00:00" >> $i.job
echo " #PBS -V" >> $i.job
echo " #PBS -M shihcheng.guo@gmail.com" >> $i.job
echo " #PBS -m abe" >> $i.job
echo " #PBS -A k4zhang-group" >> $i.job
echo " cd \$PBS_O_WORKDIR" >> $i.job
echo " perl ~/bin/methhaptools.pl --input $i --bed ../../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i" >> $i.job
qsub $i.job
done

for i in `ls *hapInfo.txt`
do
echo ""
done

grep -P "chr1" /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed

awk '{print $1,$2,$3,$1":"$2"-"$3}' /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed | grep -P "chr1" >  mm9.mhb.0.5.chr1.bed
awk '{print $1,$2,$3,$1":"$2"-"$3}' /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.4.bed | grep -P "chr1" >  mm9.mhb.0.4.chr1.bed
awk '{print $1,$2,$3,$1":"$2"-"$3}' /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.3.bed | grep -P "chr1" >  mm9.mhb.0.3.chr1.bed
perl ~/bin/haptools.pl --input SRR1248477.hapInfo.txt --bed mm9.mhb.0.5.chr1.bed --output SRR1248477.hapInfo.txt.R5
perl ~/bin/haptools.pl --input SRR1248477.hapInfo.txt --bed mm9.mhb.0.4.chr1.bed --output SRR1248477.hapInfo.txt.R4
perl ~/bin/haptools.pl --input SRR1248477.hapInfo.txt --bed mm9.mhb.0.3.chr1.bed --output SRR1248477.hapInfo.txt.R3
bigWigAverageOverBed /home/shg047/oasis/Alice/mouse/hapinfo/Mouse_ESC.meth.bw mm9.mhb.0.5.chr1.bed SRR1248477.hapInfo.txt.R5.tab
bigWigAverageOverBed /home/shg047/oasis/Alice/mouse/hapinfo/Mouse_ESC.meth.bw mm9.mhb.0.3.chr1.bed SRR1248477.hapInfo.txt.R3.tab
bigWigAverageOverBed /home/shg047/oasis/Alice/mouse/hapinfo/Mouse_ESC.meth.bw mm9.mhb.0.4.chr1.bed SRR1248477.hapInfo.txt.R4.tab 


/home/shg047/oasis/Alice/mouse/scMeth/hapinfo

/home/shg047/oasis/Roadmap/tfbs

/media/12TB_ext/DD_Ext12T/Capture_MONOD/141216_HiSeqRapidRun/BAMfiles

for i in `ls *bw`
do 
bigWigAverageOverBed $i ../tfbs/wgEncodeRegTfbsClusteredUnique.bed $i.tab
done


cd /home/shg047/oasis/Roadmap/wig
mkdir mhb
for i in `ls *bw`
do 
bigWigAverageOverBed $i /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed ./mhb/$i.mhb.tab
done

 

R CMD INSTALL openssl_0.9.6.tar.gz --configure-vars='INCLUDE_DIR=/usr/bin/pkg-config LIB_DIR=	'
libssl-dev
system("wget https://cran.r-project.org/src/contrib/xml2_1.1.1.tar.gz")
install.packages("xml2_1.1.1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/selectr_0.3-1.tar.gz")
install.packages("selectr_0.3-1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz")
install.packages("magrittr_1.5.tar.gz")
install.packages("rvest_0.3.2.tar.gz")
install.packages("BatchGetSymbols_1.1.tar.gz")

cd /home/shg047/oasis/Jenkinson2017NG/methyfreq


cd /home/shg047/oasis/Jenkinson2017NG/methyfreq
coverage2cytosine --merge_CpG --minDepth 6 --gzip --zero_based --genome_folder /home/shg047/oasis/db/hg19/ -o SRR3263674.test SRR3263674_1_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz


coverage2cytosine --merge_CpG --minDepth --zero_based --gzip --genome_folder $BismarkRefereDb -o SRR3263674.mergeCpG.bed $sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz



for i in {1..22} X Y M
do
cd /home/shg047/oasis/db/mm9
perl ~/bin/chrosomeCut.pl chr$i 10000 >> mm9.cut10k.bed
done 
wc -l /home/shg047/oasis/db/mm9/mm9.cut10k.bed
cd /home/shg047/oasis/Alice/mouse/scMeth/sortbam
perl /home/shg047/bin/samInfoPrep4Bam2Hapinfo.pl /home/shg047/oasis/db/mm9/mm9.cut10k.bed > saminfo.txt
perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt NON bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/mm9.CpG.positions.txt


# touch all the files 
find . -exec touch {} \;

cp /media/Home_Raid1/shg047/work/Alice/mouse/mhb/mm9.mhb.0.5.bed

for i in `ls *hapInfo.txt`
do
name=${i%%.*}
echo $name $i > $name.input
perl ../hapinfo2mhl.pl $name.input mm9.mhb.0.5.bed > ./MHB/$name.mhb.mhl
done




for i in `ls *hapInfo.txt`
do
name=${i%%.*}
done


scp *.hapInfo.txt shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS1/Alice/mouse/hapinfo/public/
scp SRX209459.hapInfo.txt
scp SRX209458.hapInfo.txt


perl ~/bin/smartMethSRR.pl SRP072071.txt 33 submit
perl ~/bin/smartMethSRR.pl SRP072075.txt 33 submit
perl ~/bin/smartMethSRR.pl SRP072078.txt 33 submit
perl ~/bin/smartMethSRR.pl SRP072141.txt 33 submit

perl ~/bin/smartMethSRR.pl SRP072071.txt 33 none
perl ~/bin/smartMethSRR.pl SRP072075.txt 33 none
perl ~/bin/smartMethSRR.pl SRP072078.txt 33 none
perl ~/bin/smartMethSRR.pl SRP072141.txt 33 none



bedtools shuffle -i /media/Home_Raid1/zhl002/NAS1/WGBS/MHB/Galonska_bivalentdomain.bed -excl

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from  hg19.chromInfo" > hg19.chrom.sizes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from  mm9.chromInfo" > mm9.chrom.sizes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from  mm10.chromInfo" > mm10.chrom.sizes

perl hapinfo2mhl.pl input.txt interest2432.bed > output24321.mhl &
perl hapinfo2mhl.pl input.txt Galonska_bivalentdomain.bed > output3004.mhl &
perl methHMHSum.pl input.txt interest2432.bed output24321 &
perl methHMHSum.pl input.txt Galonska_bivalentdomain.bed output3004 &

cd /home/shg047/oasis/Alice/mouse/kun
perl hapinfo2mhl.pl input.txt interest2432.bed > output24321.mhl 

cd /home/shg047/oasis/Alice/mouse/kun
perl hapinfo2mhl.pl input.txt /media/Home_Raid1/zhl002/NAS1/WGBS/MHB/Galonska_bivalentdomain.bed > output3004.mhl &


cron
# sudo vim /etc/crontab
# 17 18 * * * root bash /media/Home_Raid1/shg047/bak/bak/update.sh
# method 1
sudo vim /etc/crontab
37 16   * * *   root    bash /media/Home_Raid1/shg047/bak/bak/bak.sh
# method 2
sudo cp /media/Home_Raid1/shg047/bak/bak/bak.sh /etc/cron.daily/bak
# then restart 
sudo service cron restart
# check cron auto-run list


17 18 * * * root bash /media/Home_Raid1/shg047/bak/bak/update.sh

01 * * * * root echo "This command is run at one min past every hour"
17 8 * * * root echo "This command is run daily at 8:17 am"
17 20 * * * root echo "This command is run daily at 8:17 pm"
00 4 * * 0 root echo "This command is run at 4 am every Sunday"
* 4 * * Sun root echo "So is this"
42 4 1 * * root echo "This command is run 4:42 am every 1st of the month"
01 * 19 07 * root echo "This command is run hourly on the 19th of July"

load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/morehiseq_020617/ipsnt_hiseq_serum.RData")
load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/mips_nt/seurat_allsample_2.RData")
load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/hiseq_020617/seurat_analysis/ipsnt_hiseq_3.RData")
data<-pbmc@scale.data
id1=unlist(lapply(colnames(data),function(x) unlist(strsplit(x,"-"))[2]))
id2=pbmc@ident
group1<-c(which(id2==1),which(id2==2),which(id2==4),which(id2==6))
group2<-c(which(id2==3),which(id2==4),which(id2==7),which(id2==8))
data1=data[,group1]
data2=data[,group2]
Entropy1<-apply(data1,1,function(x) rnashannon(x))
Entropy2<-apply(data2,1,function(x) rnashannon(x))


R2<-read.table("/media/Home_Raid1/zhl002/NAS1/WGBS/ntips.pvalue.R2.txt",head=T,row.names=1)
bed<-data.frame(cor2bed(rownames(R2)),ipsR2[,1],ipsR2[,3])
bed2<-bed2gene(bed,refbed)[,c(1,2,3,4,5,11,12)]
colnames(bed2)<-c("CHR","START","END","iPS-R2","SCNT-R2","Gene","Group")
write.table(bed2,file="LD-R2-table.txt",col.names =NA,row.names =T,sep="\t",quote=F)
LDR2<-read.table("/media/NAS3_volume2/shg047/Alice/mouse/hapinfo/kun/plan3/LD-R2-table.txt")

/media/Home_Raid1/shg047/work/db/mm9/whyte_EN.bed
setwd("/media/NAS1/ZL_NAS1/WGBS")
data=read.table("markers_mhl_regions.txt",head=F)
head(data)
par(mfrow=c(3,1))
boxplot(subset(data,data[,20]=="pluripotent-serum")[,c(4:7)])
boxplot(subset(data,data[,20]=="diff-serum")[,c(4:7)])
boxplot(subset(data,data[,20]=="diff-mix")[,c(4:7)])


source("http://bioconductor.org/biocLite.R")
biocLite("DeconRNASeq")
library("DeconRNASeq")

ES Enhancer

for i in `ls *job`
do
sh $i &
done

perl ~/bin/smartMethSRR.pl SRP072071.txt 33 non
perl ~/bin/smartMethSRR.pl SRP072075.txt 33 non
perl ~/bin/smartMethSRR.pl SRP072078.txt 33 non
perl ~/bin/smartMethSRR.pl SRP072141.txt 33 non

for i in `ls *hapInfo.txt`
do
grep 34303551 $i > ./test/$i.head
done

perl methHMHSum.pl input.txt interest2432.bed output24321
/media/Home_Raid1/zhl002/NAS1/RNA_seq/hiseq_020617/seurat_analysis
/home/shg047/oasis/db/mm9
ipsnt_33k_4.RData
pbmc33k.merged@data
ID=1,2,3,4
indx= 4 5 9 10
SRR3269805_2.fastq.gz

/media/Home_Raid1/shg047/NAS1/Alice/mouse/hapinfo/jpg

perl ~/bin/bam2hapInfo2PBS.pl  saminfo.txt submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

File:CapseqSaminfoConfigBAM2hapinfo2017.txt 

awk '{ sum += $3-$2; n++ } END { if (n > 0) print sum / n; }' CapSeq.bed
samtools view PCT-6.sorted.clipped.bam |awk '{print $3,$4,$4+100}' OFS="\t" | bedtools merge -d 100 -i  - > CapSeqqMerge.bed &

for i in ls *bam
do
samtools view $i | awk "{print $3,$4,$4+100}" 
done

samtools view 


sudo vim /etc/init.d/lampp
#!/bin/bash
/opt/lampp/lampp start
sudo update-rc.d lampp defaults
insserv: warning: script 'K01lampp' missing LSB tags and overrides
insserv: warning: script 'lampp' missing LSB tags and overrides

/media/12TB_ext/DD_Ext12T/Capture_MONOD/141216_HiSeqRapidRun/BAMfiles

cd ..
for i in `ls *bam.hapInfo.txt`
do
head -n 5000 $i > ./test/$i.head
done
cd ./test
perl hapinfo2bed.pl Indx04.sortc.bam.hapInfo.txt.head | sort -k1,1 -k2,2n | bedtools merge -i - > Indx04.bed
bedtools intersect -wa -a ~/work/db/mm9/mm9.refGene.bed -b Indx04.bed > Indx04.input.bed
perl ../methHMHSum.pl input.txt Indx04.input.bed output 


perl methHMHSum.pl input.txt interest12.bed output12-2 &
perl methHMHSum.pl input.txt interest-151.bed output151 &
perl methHMHSum.pl input.txt interest-548.bed output548 &


for i in `ls *hapInfo.txt`
do
grep 34303551 $i > ./test/$i.head
done

34303551

/media/Home_Raid1/shg047/NAS1/Alice/mouse/hapinfo
my ($bin,$NM,$chr,$strand,$txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds,undef,$genesymbol,undef) = split /\t/;
	
chr10:100000000-100010000       CCTCTCT 1       100000149,100000169,100000204,100000238,100000301,100000306,100000361
/media/Home_Raid1/shg047/NAS1/RRBS
scp guoshicheng@gate1.picb.ac.cn:/picb/humpopg6/guoshicheng/lvzhen/*.gz ./

the website:
https://plus.google.com/112774990052239691835
Dongqing Tan
593 6TH AVE
SAN FRANCISCO, CA 94118-3816
United States
Phone number: 415-373-7338

for i in `ls *bam.hapInfo.txt`
do
perl hapinfosummary.pl $i >> result.txt
done

for i in {1..11}
do
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw ./
done

for i in {1..2}
do
http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas$i/tracks_hg19/Human_NormalPancreas$i.meth.bw
done


for i in `ls *bam.hapInfo.txt`
do
head -n 50 $i > ./test/$i.head
done


/media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/public
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDzCsAX+XaSJsDJR71BlYeOLPGgmz5JtpuvL1mYwOsvEZGNA7Y7ufxZE3gDHTZ/fKjfVbewH2FLu+2It6d9LayEy8JQYCcb0oqY1nXHZMdf8bpIfPGz9GJabOHLV6iqx0avpVg64gLtLiWsNV9AfBWd8gjjBTPK33M69SH4/QXAi9zglu6Y325M7FdhyM/7GH/5rIEjd6P5GM6LWJQQkAKpGJf6RcYg9bSyy239SlAFkEJMnBTYCo9xOJYHcyPQAuX6WKiZd/3Sxlqce2zpOIC6O7Q6Q07HWa5aO54R9YXy1pz5NBbn9EtAOZIcS1OuklD2zs8hsjRkijWHlX20Ahbn shg047@tscc-login1.sdsc.edu
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC33azQGufroerEJp+PuQmrskhbFa6n9sxfXDoKeNI0yPiJCMwvSdxlsZ6k+BqB+84tcTssbHPSVW1EJ5bVRFt6GuDbZoaewHffCkd7hVYZCvZFFtU/o69O+3i0MzWpfO7nwY5Dfhg4rO33Xh8ijKTL4PdAHeoqQ3Cf5aYRfe3f71MEdS7ObfC1FL5h1vKtIsL/IReJ6EUp6h37wyqtnkIQqmkIhJHYc/zCh+CJVgOU6RiCo7rbJBFPozHEITBfrDzhe0azO3+eGdzcerLoKok8eqSc9AKCmSYUW4kdhnuN1j/nOGUaf3i0jZZUWyBldxykxh8Tob1xrgSnL9gA9PfD shg047@genomeMiner

for i in {1..22} X Y
do
grep chr$i HsGenome19.CpG.positions.txt > HsGenome19.CpG.positions.chr$i.txt
done

chr10:116387101-116387201

for i in `ls *hapInfo.txt`
do
perl ~/bin/hapinfo2mf.pl $i > ./methyfreq/$i.bedgraph 
echo $i
done

find . -name "*" -type d -exec rm -r "{}" \;

/media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/public
You can also make a payment via IVR (Interactive Voice Response) by calling 1-800-892-4357 without any fees or charge.

cp /media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/mESC_enhancer.bed ~/work/db/mm9
cp /media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/Galonska_bivalentdomain.bed ~/work/db/mm9

for i in mESC_enhancer.bed Galonska_bivalentdomain.bed mm9.Enhancer.bed mm9.Promoter.bed mm9.Exon.bed mm9.Intron.bed
do
for j in miPS_B3 miPS_1E12P20 SCNT_NB3.BED SCNT_B12 B6_mESC HA.129.mESC ST_mESC_mix
do
echo $i $j
perl /media/Home_Raid1/shg047/work/Alice/bin/AverageValueNearbyFunctionRegion.pl /media/Home_Raid1/shg047/work/db/mm9/$i  $j.BED.txt.trim  $j.$i.dist.txt
done
done

/media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/Galonska_bivalentdomain.bed

miPS_B3.BED.txt.trim

ls miPS_B3.BED.txt.trim
ls miPS_1E12P20.BED.txt.trim
ls SCNT_NB3.BED.txt.trim
ls SCNT_B12.BED.txt.trim
ls B6_mESC.BED.txt.trim
ls HA.129.mESC.BED.txt.trim
ls ST_mESC_mix.BED.txt.trim

echo "# Twin Methylation Dataset Analysis" >> README.md
git init
git add README.md
git commit -m "project start"

/media/Home_Raid1/zhl002/NAS3/ZL_LTS33/mouse_WGBS/bedfiles

for i in `ls *.bed`
do
awk '{if ($4>5) print $1,$2,$3}' OFS="\t" $i > $i.txt
done

#!/usr/bin/perl
use strict;
my @file=glob("*.bed");
foreach my $file(@file){
open F,$file;
while(<F>){
my($chr,$start,$end,$num)=split/\s+/;
print "$chr\t$start\t$end\n" if $num>10;
}
}

/home/shg047/bak/plink/china/gzip/input.bed
cd /home/shg047/bak/plink/china/gzip/
perl -p -i -e 's/\s+/\t/g' /home/shg047/bak/plink/china/gzip/input.bed

# plan 1: 2i vs serum in SCNT
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun
cat Indx06.sortc.bam.hapInfo.txt Indx07.sortc.bam.hapInfo.txt > nt2i.HapInfo.txt
cat Indx09.sortc.bam.hapInfo.txt Indx10.sortc.bam.hapInfo.txt > ntserum.HapInfo.txt
perl ~/bin/hapinfoMerge.pl nt2i.HapInfo.txt > nt2i.HapInfo.Uni.txt
perl ~/bin/hapinfoMerge.pl ntserum.HapInfo.txt > ntserum.HapInfo.Uni.txt
http://www.rencai8.com/job_info?action=view&job_position_id=476884
mv nt2i.HapInfo.Uni.txt ./plan1
mv ntserum.HapInfo.Uni.txt ./plan1
cd ./plan1 
for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl nt2i.HapInfo.Uni.txt > nt2i.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl nt2i.HapInfo.Uni.txt.$i > R2.nt2i.HapInfo.Uni.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl ntserum.HapInfo.Uni.txt > ntserum.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ntserum.HapInfo.Uni.txt.$i > ntserum.HapInfo.Uni.txt.$i
done

# plan 2: ips (1,2,4,5) vs scnt (6,7,9,10)
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun
touch  ips.HapInfo.txt
for i in 01 02 04 05 
do
cat Indx$i.sortc.bam.hapInfo.txt >> ips.HapInfo.txt
done

touch  scnt.HapInfo.txt
for i in 06 07 09 10 
do
cat Indx$i.sortc.bam.hapInfo.txt >> scnt.HapInfo.txt
done

perl ~/bin/hapinfoMerge.pl ips.HapInfo.txt > ips.HapInfo.Uni.txt
perl ~/bin/hapinfoMerge.pl scnt.HapInfo.txt > scnt.HapInfo.Uni.txt

mv ips.HapInfo.Uni.txt ./plan2
mv scnt.HapInfo.Uni.txt ./plan2
cd ./plan2

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl ips.HapInfo.Uni.txt > ips.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ips.HapInfo.Uni.txt.$i > R2.ips.HapInfo.Uni.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl scnt.HapInfo.Uni.txt > scnt.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl scnt.HapInfo.Uni.txt.$i > R2.scnt.HapInfo.Uni.txt.$i
done

# plan 3: ips(4,5) vs scnt(9,10)
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun
touch ips.HapInfo.txt
for i in 04 05 
do
cat Indx$i.sortc.bam.hapInfo.txt >> ips.HapInfo.txt
done

touch scnt.HapInfo.txt
for i in 09 10 
do
cat Indx$i.sortc.bam.hapInfo.txt >> scnt.HapInfo.txt
done

perl ~/bin/hapinfoMerge.pl ips.HapInfo.txt > ips.HapInfo.Uni.txt
perl ~/bin/hapinfoMerge.pl scnt.HapInfo.txt > scnt.HapInfo.Uni.txt

mv ips.HapInfo.Uni.txt ./plan3
mv ips.HapInfo.Uni.txt ./plan3
cd ./plan3

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl ips.HapInfo.Uni.txt > ips.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ips.HapInfo.Uni.txt.$i > R2.ips.HapInfo.Uni.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl scnt.HapInfo.Uni.txt > scnt.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl scnt.HapInfo.Uni.txt.$i > R2.scnt.HapInfo.Uni.txt.$i
done

cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun/plan1
perl ../../../../bin/R2matrix.pl > plan1.R2.txt
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun/plan2
perl ../../../../bin/R2matrix.pl > plan2.R2.txt
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun/plan3
perl ../../../../bin/R2matrix.pl > plan3.R2.txt


 for i in `ls *.bam`; 
 do 
 touch $i.bam2mf.job
 echo '#!/bin/csh' >$i.bam2mf.job
 echo "#PBS -N $i" >> $i.bam2mf.job
 echo "#PBS -q hotel" >> $i.bam2mf.job
 echo "#PBS -l nodes=1:ppn=1" >> $i.bam2mf.job
 echo "#PBS -l walltime=168:00:00" >> $i.bam2mf.job
 echo "#PBS -M shihcheng.guo@gmail.com" >> $i.bam2mf.job
 echo "#PBS -m abe" >> $i.bam2mf.job
 echo "#PBS -A k4zhang-group" >> $i.bam2mf.job
 echo "cd $(pwd)" >> $i.bam2mf.job
 echo bismark_methylation_extractor --single-end --bedGraph --cutoff 1 --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../sortBam/$i >> $i.bam2mf.job; 
 qsub $i.bam2mf.job
 done
 
perl ~/bin/hapinfo2bedgraph.pl CTR97.hapInfo.txt > CTR97.bedgraph

head /home/shg047/oasis/DennisLo2015/methyfreq/CTR97.merged_CpG_evidence.cov

chr10:100004428-100004429       80.000000
chr10	100004428	100004429	chr10:100004428-100004429

grep 100004429 /home/shg047/oasis/DennisLo2015/methyfreq/CTR97.merged_CpG_evidence.cov
grep 100004428 CTR97.bedgraph

samtools tview ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam ~/oasis/db/hg19/hg19.fa
samtools view -bh ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam chr10:100004428-100004429 -o CTR97.chr10-100004428-100004429.bam
samtools view ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam chr10:100004428-100004429

perl ~/bin/mergedBam2hapInfoV2.pl test.bed CTR97.chr10-100004428-100004429.bam bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt > CTR97.hapInfo.test.txt


/home/shg047/oasis/DennisLo2015/sortBam

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed /oasis/tscc/scratch/shg047/Alice/mouse/sortBam/SRX1019866.sortc.bam bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/mm9.CpG.positions.txt > ../hapinfo/SRX1019866.sortc.bam.hapInfo.txt

cd /media/Home_Raid1/shg047/oasis/Alice/mouse/hapinfo/
bedtools intersect -wa -a mm9.mhb.0.5.bed -b /media/Home_Raid1/shg047/work/db/mm9/CpGI.mm9.txt  | wc -l


/media/Home_Raid1/shg047/oasis/Alice/mouse/hapinfo/mm9.mhb.0.5.bed

/home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed


/media/Home_Raid1/shg047/work/db/mm9

表观遗传学科研资讯
crystalqing@yahoo.com
perl hapinfo2mhb.pl -f input.txt -t 0.5 -o mm9.mhb.multinput.bed &
perl ~/bin/hapinfo2mhb.pl mm9.HapInfo.txt.SumUniq 0.5 > mm9.mhb.0.5.bed &
perl ~/bin/hapinfo2mhb.pl mm9.HapInfo.txt.SumUniq 0.3 > mm9.mhb.0.3.bed &


ln -s  /opt/biotools/htseq/bin/htseq-count  htseq-count
ln -s  /opt/biotools/htseq/bin/htseq-qa htseq-qa

/home/shg047/oasis/Alice/mouse/hapinfo


sftp -o Port:8777 'user@domain.com'@domain.com
sftp -o Port:8777 -o User=user@domain.com domain.com
sshpass -p 'Gsc$$8343383' ssh shg047@tscc-login.sdsc.edu

wget --ftp-user='zhang' --ftp-password='HDIIWpP' ftp://169.228.63.66/
wget -r --user='zhang' --password='HDIIWpP' ftp://169.228.63.66/

/home/shg047/software/gmp-4.2.4/build/.libs
/home/shg047/software/gmp-4.2.4
/home/shg047/software/mpfr-2.4.1/build/.libs
/home/shg047/software/mpfr-2.4.1/build

/home/shg047/software/mpfr-3.1.4

export C_INCLUDE_PATH=/home/shg047/software/mpfr-2.4.1/build
export LD_LIBRARY_PATH=/home/shg047/software/mpfr-2.4.1/build/.libs
export LIBRARY_PATH=$LD_LIBRARY_PATH

cd /home/shg047/software/mpc-1.0.3/build
export LD_LIBRARY_PATH=/home/shg047/software/mpfr-2.4.1/build/.libs
../configure --with-gmp-lib=/home/shg047/software/gmp-4.2.4/build/.libs


configure: error: libmpfr not found or uses a different ABI (including static vs shared).


perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/db/mm9/mm9.cut10K.bed

cd /home/shg047/oasis/Ziller2013/methyfreq
coverage2cytosine --merge_CpG  --genome_folder  -o 


perl ~/bin/smartbismark.pl --input saminfo.txt --server TSCC --queue hotel --genome mm9
for i in `ls *pbs`; do qsub $i; done
-rw-r--r--  1 shg047 k4zhang-group  1.1G Feb  6 00:30 Indx06.merged.bam_trimmed.fq
-rw-r--r--  1 shg047 k4zhang-group  39G Jan 27 13:13 Indx05.merged.bam.fastq


/home/shg047/oasis/DennisLo2015/methyfreq

coverage2cytosine --merge_CpG --genome_folder ~/oasis/db/hg19/align/bismark/ -o $sample.tmp $file


wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xzvf chromFa.tar.gz
for i in `ls *fa`
do
perl ../../bin/cgpositionFinder.pl $i &
done

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFaMasked.tar.gz

trim_galore --phred33 --fastqc --illumina --rrbs SRR1286404.fastq.gz --output_dir ../fastq_trim
bismark --bowtie2 --phred33-quals --multicore 2 --fastq -L 32 -N 0 -D 5 -R 1 /home/shg047/oasis/db/hg19/align/bismark ../fastq_trim/SRR1286404_trimmed.fq.gz -o ../bam
bismark_methylation_extractor --multicore 3 --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../bedgraph  ../bam/SRR1286404_trimmed_bismark_bt2.bam
samtools sort ../bam/SRR1286404_trimmed_bismark_bt2.bam -o ../sortbam/SRR1286404_trimmed_bismark_bt2.sort.bam
samtools index ../sortbam/SRR1286404_trimmed_bismark_bt2.sort.bam

bismark --bowtie2 --phred33-quals --fastq /media/Home_Raid1/shg047/work/db/mm9/bismark  SRX1091397.fastq

cd /media/Home_Raid1/shg047/oasis/Alice/mouse/test

bismark --bowtie2 --phred33-quals --fastq /media/Home_Raid1/shg047/work/db/mm9/bismark  SRX080191.fastq.gz

chr1:98109663-98109763
# check chr1:100006294-100006356 in SRX080191 what happened, why 
# SRX080191.hapInfo.txt:chr1:100000000-100010000  CC      3       100006294,100006294
samtools view SRX080191.sort.bam | less -S 

cd /home/shg047/oasis/Alice/mouse/hapinfo
grep 100006294 SRX080191.hapInfo.txt
samtools view -b -o SRX080191.short.bam SRX080191.sort.bam chr1:100006294-100006356 
samtools tview SRX080191.sort.bam /home/shg047/db/mm9/mm9.fa
samtools tview SRX080191.sort.bam /home/shg047/oasis/db/mm9/mm9.fa
samtools tview Indx01.merged.bam /home/shg047/oasis/db/mm9/mm9.fa


samtools tview SRX080191.sort.bam /home/shg047/db/hg19/hg19.fa
samtools tview SRX080191.sort.bam /home/shg047/oasis/db/hg19/hg19.fa

perl R2matrix.pl > MatrixR2.txt


chr1:100000000-100010000        CC      3       100006294,100006294
chr1:100000000-100010000        CCCC    1       100006294,100006294,100006350,100006350
chr1:100000000-100010000        CCCCCC  1       100006294,100006294,100006350,100006350,100006356,100006356
chr1:100000000-100010000        TT      1       100006294,100006294

samtools tview SRX080191.short.bam ~/work/db/


for i in {1..19} X Y
do
grep -w chr$i Mouse.mm9.HapInfo.txt > Mouse.mm9.chr$i.Hapinfo.txt 
perl ~/bin/hapinfoMerge.pl Mouse.mm9.chr$i.Hapinfo.txt 2> chr$i.err 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt 0.3 > mm9.mhb.chr$i.0.3.txt 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt 0.5 > mm9.mhb.chr$i.0.5.txt  
done

for i in {1..19} X Y
do
grep -w chr$i mm9.HapInfo.txt > Mouse.mm9.chr$i.Hapinfo.txt 
perl ~/bin/hapinfoMerge.pl Mouse.mm9.chr$i.Hapinfo.txt 2> chr$i.err 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt.uni.txt 0.3 > mm9.mhb.chr$i.0.3.txt 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt.uni.txt 0.5 > mm9.mhb.chr$i.0.5.txt  
done


# scRNA-seq could explain methylation haplotype variation

iPS<-read.table("/oasis/tscc/scratch/shg047/Alice/mouse/alice/scRNA/ipsc_R2high.bed")


# 01/27/2017
for i in `ls *pbs`; do qsub $i; done
perl ~/bin/smartbismark.pl --input saminfo.txt --genome mm9 --server TSCC --queue glean
for i in {1..22} X Y M
do
perl ~/bin/genomecut.pl chr$i 10000 /home/shg047/db/mm9/mm9.chrom.sizes > mm9.$i.10k.bed
done
cat mm9.chr*.10k.bed > mm9.chr*.10K.bed

/home/shg047/db/mm9/mm9.cut10K.bed
perl R2matrix.pl > MatrixR2.txt
cat Indx01.hapInfo.txt Indx02.hapInfo.txt Indx04.hapInfo.txt Indx05.hapInfo.txt > iPS.HapInfo.txt
cat Indx06.hapInfo.txt Indx07.hapInfo.txt Indx09.hapInfo.txt Indx10.hapInfo.txt > SCNT.HapInfo.txt
perl ~/bin/hapinfoMerge.pl iPS.HapInfo.txt
perl ~/bin/hapinfoMerge.pl SCNT.HapInfo.txt

mv SCNT.HapInfo.txt.SumUniq SCNT.hapinfo.txt
mv iPS.HapInfo.txt.SumUniq iPS.hapinfo.txt

perl ~/bin/hapinfo2BlocAvgR2.pl iPS.HapInfo.txt.SumUniq > iPS.MHB.R2.txt
perl ~/bin/hapinfo2BlocAvgR2.pl SCNT.HapInfo.txt.SumUniq > SCNT.MHB.R2.txt

# 01/27/2017
for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl SCNT.hapinfo.txt > SCNT.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl SCNT.hapinfo.txt.$i > R2.SCNT.hapinfo.txt.$i
perl ~/bin/randomSampleFromHaploInfo.pl iPS.hapinfo.txt > iPS.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl iPS.hapinfo.txt.$i > R2.iPS.hapinfo.txt.$i
done


# compare two file list in fold A and fold B by diff and ls
diff <(ls ../methyfreq/*cov.gz | awk -F"/" '{print $3}'| awk -F"_" '{print $1}') <(ls *sort.bam | awk -F"_" '{print $1}')

# compare two file list in fold A and fold B by diff and ls
ls ../methyfreq/*cov.gz | awk -F"/" '{print $3}'| awk -F"_" '{print $1}' > a
ls *sort.bam | awk -F"_" '{print $1}' > b
grep -v -f a b 

# compare two file list in fold A and fold B by diff and ls and then qsub job which are not finished. 
grep -v -f a b | awk '{print "qsub "$1"\*.job"}'

diff <(ls ../methyfreq/*cov.gz | awk -F"/" '{print $3}'| awk -F"_" '{print $1}') <(ls *sort.bam | awk -F"_" '{print $1}')
CTR154_trimmed.fq.gz_bismark_bt2.sort.bam

for i in `ls Indx*bam`
do
samtools fastq $i > $i.fastq &
done

gzip $i &

/media/Home_Raid1/dinh/NAS3_volume1_mnt/shg047/GEO/GSE63123_RAW.tar
tar xvf GSE63123_RAW.tar

1, cluster all scRNA by subset genes plu-diff
2, select good cells and separate ips and scnt
3, find differential express genes between ips and scnt
4, show previous knowledge base on this data
5, check Methylation to these DEG
6, find methylation biomarker/matric to show the difference between ips and scnt

/home/shg047/work/DennisLo2015/bam/
 perl ~/software/SNPsplit/SNPsplit --bisulfite --SNP_file CTR151_trimmed.fq.gz_bismark_bt2.sort.bam
samtools view HOT162_trimmed.fq.gz_bismark_bt2.sort.bam | head -n 8000 > ./test/HOT162_trimmed.fq.gz_bismark_bt2.sort.sam
samtools view Pregnancy.8.run1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam | head -n 8000 > ./test/Pregnancy.8.run1.read1_val_1.fq.gz_bismark_bt2_pe.sort.sam

sh CTR84_trimmed.fq.gz_bismark_bt2.sort.bam.bam2mf.job &
sh LTP1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam.bam2mf.job &


CTR132_trimmed.fq.gz_bismark_bt2.sort.sortn.bam

perl ~/bin/bismarkbam2methyfreq.pl --input saminfo.txt --server TSCC --genome hg19 --queue glean 

bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 8 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/HOT215_trimmed.fq.gz_bismark_bt2.sort.bam
perl ~/bin/hapinfo2mhb.pl merge.txt.SumUniq 0.3 > mm9.MHB.0.3.txt   
perl ~/bin/hapinfo2mhb.pl merge.txt.SumUniq 0.5 > mm9.MHB.0.5.txt 
mkdir ../mhb
mv mm9.MHB.0.3.txt ../mhb
mv mm9.MHB.0.5.txt ../mhb/
Added a short File description table in Methylation Haplotype Analysis User Guide.docx
sudo ufw allow 8787
lynx http://<ip>:8787
change chrome to firfox
http:132.239.25.238:8787
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949201/SRR949201_1.fastq.gz 
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949201/SRR949201_2.fastq.gz 
2017-01-17 13:28:06 (0.00 B/s) - "SRR949201_1.fastq.gz" saved [16274079739]
14679170302 Jan 17 14:46 SRR949201_2.fastq.gz
15350187717 Jan 17 15:07 SRR949201_1.fastq.gz
15550948069 Jan 16 18:50 SRR949201_2.fastq.gz
16274079739 Jan 16 18:50 SRR949201_1.fastq.gz
How to prepare bed format for human and mouse refGene

wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
gunzip -f *.gz
 
fetchChromSizes hg19 > hg19.chom.sizes
fetchChromSizes hg38 > hg38.chrom.sizes
fetchChromSizes mm9 > mm9.chrom.sizes

 
 
#PBS -N Indx04.job
#PBS -q glean
#PBS -l nodes=1:ppn=6
#PBS -l walltime=168:00:00
#PBS -o Indx04.log
#PBS -e Indx04.err
#PBS -V
#PBS -M shicheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group

 
 
 fetchChromSizes mm10 > mm10.chrom.sizes
 
 perl ../mm9/refGene2bed.pl mm9.refGene.txt mm9.chrom.sizes > mm9.refGene.bed
 perl ../mm9/refGene2bed.pl mm10.refGene.txt mm10.chrom.sizes > mm10.refGene.bed
 perl ../mm9/refGene2bed.pl hg19.refGene.txt hg19.chrom.sizes > hg19.refGene.bed
 perl ../mm9/refGene2bed.pl hg38.refGene.txt hg38.chrom.sizes > hg38.refGene.bed
 mv mm9.refGene.bed ../mm9
 cd ../mm9
 
 
 wc -l *.refGene.bed
 
  826320 hg19_refGene.bed
  1372108 hg19.refGene.bed
  1450953 hg38.refGene.bed
   798143 mm10.refGene.bed
   797710 mm9.refGene.bed

  1372104 hg19.refGene.bed
  1449295 hg38.refGene.bed
   798075 mm10.refGene.bed
   797710 mm9.refGene.bed
  4417184 total
 
  1372108 hg19.refGene.bed
  1450953 hg38.refGene.bed
   798143 mm10.refGene.bed
   797710 mm9.refGene.bed
  4418914 total
  
    1372071 hg19.refGene.bed
  1447078 hg38.refGene.bed
   797990 mm10.refGene.bed
   797657 mm9.refGene.bed
  4414796 total

  14679170302 Jan 17 14:46 SRR949201_2.fastq.gz
-rw-r--r-- 1 shg047 k4zhang-group 15350187717 Jan 17 15:07 SRR949201_1.fastq.gz



 wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
 wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
 wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
 wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

 for i in `ls *.gz`; 
 do 
 SRR=${i%%_*}; 00
 echo $SRR >> list.txt; 
 done
 
 for j in `sort -u list.txt`
 do
 vdb-validate $j
 done
 
 for i in `ls *.gz`; 
 do 
 gunzip -t $i 2 > $i.err &
 done
 find . -name "*gz.err" -type f -size +0c -exec ls -larth {} \;
 

cutadapt: error: In read named 'SRR949201.44244409 D1JR8ACXX130107:7:1108:14746:100581 length=99': length of quality sequence (90) and length of read (99) do not match

less SRR949201_1.fastq.gz |grep -n -A4 'SRR949201.44244409 D1JR8ACXX130107:7:1108:14746:100581 length=99' 
	
perl ~/bin/smartbismark.pl --input PRJNA201480.txt --submit no --genome hg19 --server TSCC
qsub SRR949201.pbs
qsub SRR949202.pbs
qsub SRR949210.pbs
qsub SRR949211.pbs

sh SRR949201.job &
sh SRR949202.job &
sh SRR949210.job &
sh SRR949211.job &

for i in `ls *.gz`
do
md5sum $i >> md5sum.txt
md5sum $i >> md5sum.txt
done

SRR949202_1 missed in fastq_trim folder
_2 missed in fastq_trim folder
SRR949210_1 missed in fastq_trim folder
_2 missed in fastq_trim folder
SRR949211_1 missed in fastq_trim folder
SRR949211_2 missed in fastq_trim folder


coverage2cytosine --merge_CpG SRR949207_1_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz --genome_folder ~/oasis/db/hg19/align/bismark/ -o SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph

cd /home/shg047/software/R-3.3.2
./configure --prefix=/home/shg047/software/R-3.3.2 '--with-cairo' \
 '--with-jpeglib' '--with-readline' '--with-tcltk' '--with-x=no'\
 '--with-blas' '--with-lapack' '--enable-R-profiling' '--with-tiff=yes'\
 '--enable-R-shlib'\
 '--enable-memory-profiling'
 make clean
 make
 
 
Name is names what the name 
 
Clear idea, clean story, base on truth and get some another truth. (drugs and interaction)
Our value is fix some gap in the academic fileds.
We would better do some thing beyond context dependent. how to interpret the result or study with context dependent. 
 => human mutation. 
1, small group and train and fight yourself
2, train yourself to find the bar.
3, you should know how to fight alone.  
4, build your own contribution. also you need your project huge. that means you need more energy to have some light on
5, people or paper or boht. 

cd /media/Home_Raid1/shg047/software/InfiniumPurify
for i in `ls jhu*`
do
python InfiniumPurify.py -f $i -c LUAD 
done


Bowtie with parameters "-q --phred33-quals -n 1 -e 99999999 -l 25 -I 1 -X 2000 -a -m 15 -S -p 6",
E(i,j)=log2(TPMi,j/10+1), where TPMi,j refers to transcript-per-million (TPM) for gene i in sample j, as calculated by RSEM


head -n 20 SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
head -n 20 SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov


wc -l  SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov


1731880
15609274

grep chrUn_gl000214 SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
grep SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov

head SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
head SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov

system("cat $cov >> $SRX.bedgraph");
/home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

 #!/bin/csh
 #PBS -N MergeCOVbySRX
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group 
cd /home/shg047/oasis/Estellar2016/methyfreq
coverage2cytosine --merge_CpG SRR1035731_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.head --genome_folder ~/oasis/db/hg19/align/bismark/ -o SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph


Chase Sapphire Preferred (CSP)
Amex Centurion
Citi Prestige
Amex Premier Rewards Gold (PRG)
Amex Platinum
Chase Ritz Carlton
Chase Palladium
Amex SPG
Discover it


 perl ~/bin/samInfoPrep4Bam2Hapinfo.pl /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | grep -v PPP3CA
 
qsub SRR949205.pbs
qsub SRR949212.pbs
qsub SRR949215.pbs

curl "https://gdc-api.nci.nih.gov/legacy/files/4a8ffe0d-d7e6-4712-ad04-472955c84c77?fields=cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type&format=tsv"
curl "https://gdc-api.nci.nih.gov/legacy/files/087ec4fb-a621-4fcf-8276-1c74782bcc2c?fields=cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type&format=tsv"
 
du ./ -h --max-depth 1

for i in BISULFITE EXTENSION HYBRIDIZATION NEGATIVE NON-POLYMORPHIC NORM_A NORM_G NORM_C, NORM_T RESTORATION STAINING SPECIFICITY TARGET
do
grep $i GPL21145-48548_EPIC.txt | wc -l
done

grep NEGATIVE GPL21145-48548_EPIC.txt 

# 2016-12-30
cd /media/Home_Raid1/shg047/work/Roadmap/bw
for i in GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig.gz.bw GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.gz.bw GSM1120336_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw GSM1120337_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/CpGSF.hg19.sort.bed4 $i.tab
done

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/CpGSF.hg19.sort.bed4 $i.tab
mv $i.tab /media/Home_Raid1/shg047/work/Roadmap/bw
done

cd /media/Home_Raid1/shg047/work/Roadmap/bw
perl ~/bin/bigWigAverageOverBed2Matrix.pl > CpGI.txt

data<-read.table("CpGI.txt")
group<-unlist(lapply(rownames(data),function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(Pos=rownames(data),data,group)
head(input)
library("reshape2")
colnames(input)<-c("POS","Kidney","Right_Ventricle 1","Right_Ventricle 2","Kidney_Tumor","Kidney_Normal","Kidney_Tumor","Kidney_Normal","Esophagus","group")
input.long<-melt(input, id.vars=c("group","POS"))
library(ggplot2)
png("rmsk.CpGI.methylation.png")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA)+ coord_flip()
dev.off()
tapply(input.long$value,input.long$group,function(x) median(x,na.rm=T))
tapply(input.long$value,input.long$variable,function(x) mean(x,na.rm=T))
head(subset(input.long,group=="CpGI"))



# 2016-12-30
cd /media/Home_Raid1/shg047/work/Roadmap/bw
for i in GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig.gz.bw GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.gz.bw GSM1120336_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw GSM1120337_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
done

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
mv $i.tab /media/Home_Raid1/shg047/work/Roadmap/bw
done

cd /media/Home_Raid1/shg047/work/Roadmap/bw
perl ~/bin/bigWigAverageOverBed2Matrix.pl > Repeat.txt

data<-read.table("repeat.txt")
group<-unlist(lapply(rownames(data),function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(Pos=rownames(data),data,group)
head(input)
library("reshape2")
colnames(input)<-c("POS","Kidney","Right_Ventricle 1","Right_Ventricle 2","Kidney_Tumor","Kidney_Normal","Kidney_Tumor","Kidney_Normal","Esophagus","group")
input.long<-melt(input, id.vars=c("group","POS"))
library(ggplot2)
png("rmsk.repeat.methylation.png")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA)+ coord_flip()
dev.off()
tapply(input.long$value,input.long$group,function(x) median(x,na.rm=T))
tapply(input.long$value,input.long$variable,function(x) mean(x,na.rm=T))
head(subset(input.long,group=="CpGI"))


# 2016-12-30
cd /media/Home_Raid1/shg047/work/db/hg19
perl rmsk.pl -i rmsk.hg19 -o rmsk.hg19.bed
sort -u rmsk.hg19.bed > rmsk.hg19.bed.temp
sort -k1,1 -k2,2n rmsk.hg19.bed.temp > rmsk.hg19.bed
rm rmsk.hg19.bed.temp
awk '{print $12}' rmsk.hg19 | sort -u

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
done

perl ~/bin/ Repeat.txt

data<-read.table("Repeat.txt")
info<-read.table("/media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed")
group<-unlist(lapply(info[,4],function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(data,group)

library("reshape2")
colnames(input)<-c("T","N","T","N","group")
input.long<-melt(input, id.vars=c("group"))

library(ggplot2)
pdf("rmsk.methylation.pdf")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA,outlier.colour="white")+ coord_flip()
dev.off()
##

git pull
git status
git rebase --continue

#

http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://132.239.25.238/shg047/myhub.txt
 
ln -s /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/bw/GSM1546663_5mC-P1-T-corrected.txt.bw

/media/Home_Raid1/shg047/work/h
UCSC Over
mv /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/bw
myhub.txt

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/
for i in `ls *bw`
do
ln -s /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/$i /media/Home_Raid1/shg047/work/hub/hg19/$i
done

sudo ln -s /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/GSM1546663_5mC-P1-T-corrected.txt.bw /media/Home_Raid1/shg047/work/hub/hg19/GSM1546663_5mC-P1-T-corrected.txt.bw



# 2016-12-29
 liftOver GSE17972.hg18.bedgraph ~/work/db/hg18/hg18ToHg19.over.chain GSE17972.hg19.bedgraph tmp
 sort -u -k1,1 -k2,2n GSE17972.hg19.bedgraph > GSE17972.hg19.bedgraph.sort
 bedGraphToBigWig GSE17972.hg19.bedgraph.sort ~/work/db/hg19/hg19.chrom.sizes GSE17972.hg19.bw
 track type=bigWig color=0,0,255 visibility=2 maxHeightPixels=128:30:11 smoothingWindow=16 windowingFunction=mean name="PBMC" description="Yanhuang-methylome" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Yanhuang2010/GSE17972.hg19.bw

 liftOver GSE17972.hg18.bedgraph ~/work/db/hg18/hg18ToHg38.over.chain GSE17972.hg38.bedgraph tmp
 sort -u -k1,1 -k2,2n GSE17972.hg38.bedgraph > GSE17972.hg38.bedgraph.sort
 bedGraphToBigWig GSE17972.hg38.bedgraph.sort /media/Home_Raid1/shg047/work/db/hg38/hg38.chrom.sizes GSE17972.hg38.bw
 track type=bigWig color=0,0,255 visibility=2 maxHeightPixels=128:30:11 smoothingWindow=16 windowingFunction=mean name="Yanhuang-methylome" description="PBMC" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Yanhuang2010/GSE17972.hg38.bw
 
 GSE17972.hg19.bw
 GSE17972.hg19.bedgraph.sort
 
for i in `ls *job`
do

if [! -e ]

for i in SRR1035834 SRR1035835 SRR1035844 SRR1035831 SRR1035784 SRR1035893 SRR1035884 SRR1035843 SRR1035895 SRR1035832 SRR1035882 SRR1035881 SRR1035845 SRR1035833
do
qsub $i.fastq.download.job
done

for i in `ls *bam`
do
touch $i
done


 install CPAN
 reload cpan

 # qmap to bedgraph, liftOver hg18 to hg19 and hg38, finally sort the bedgraph 

 for i in `ls *.txt`
 do
 echo $i
 perl qmap2bedgraph.pl $i > $i.bedgraph
 done
 
 # merge liftOver 
 for i in `ls *.bedgraph`
 do
 echo $i
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.hg19 tmp
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.hg38 tmp
 done
 
 # Sort liftOver bedgraph
 for i in `ls *.bedgraph.hg19`
 do
 echo $i
 sort -k1,1 -k2,2n $i > $i.sort
 done

 for i in `ls *.bedgraph.hg38`
 do
 echo $i
 sort -k1,1 -k2,2n $i > $i.sort
 done

 # merge bedgraph by chrosome
 cat *hg18.sort > GSE17972.YanHuang.hg18.bedgraph 
 cat *hg19.sort > GSE17972.YanHuang.hg19.bedgraph 
 cat *hg38.sort > GSE17972.YanHuang.hg38.bedgraph 

 # bedgraph to bigwig 
 for i in `ls *.bedgraph`
 do
 echo $i
 bedGraphToBigWig $i ~/work/db/hg18/hg18.chrom.sizes $i.bw
 bedGraphToBigWig $i ~/work/db/hg19/hg19.chrom.sizes $i.bw
 bedGraphToBigWig $i ~/work/db/hg38/hg38.chrom.sizes $i.bw
 done

 
 
 
 
 for i in `ls *chr10*txt`
 do
 sort -k1,1 -k2,2n $i.bedgraph.hg19 > $i.bedgraph.hg19.sort
 bedGraphToBigWig $i.bedgraph.hg19 ~/work/db/hg19/hg19.chrom.sizes $i.hg19.bw
 done



 
liftOver GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph.hg19 tmp
liftOver GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph.hg38 tmp


http://132.239.25.238/shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546666_5mC-P2-N-corrected.txt.bw

shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546663_5mC-P1-T-corrected.txt.bw
shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546664_5mC-P1-N-corrected.txt.bw
shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546665_5mC-P2-T-corrected.txt.bw
shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546666_5mC-P2-N-corrected.txt.bw

http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl= 
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://132.239.25.238/shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/myhub.txt

track type=bigWig name="XX" description="YY" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546663_5mC-P1-T-corrected.txt.bw


for i in `ls *5mC-P*txt`
do
sort -k1,1 -k2,2n $i >$i.sort
perl tobedgraph.pl $i.sort 
bedGraphToBigWig $i.sort.bedgraph ~/db/hg19/hg19.chrom.sizes $i.bw 
rm $i.sort 
rm $i.sort.bedgraph
done


for i in `ls *.job`
do
sh $i &
done



GetOptions ( "input=s"   => \$input,           # string
             "submit=s"  => \$submit,          # flag
             "genome=s" => \$genome,          # string
             "server=s" => \$server)          # flag
o
perl ~/bin/smartbismark.pl --input PRJNA201480.txt --submit no --genome hg19 --server TSCC

bismark --bowtie2 --multicore 4 --phred33-quals --fastq -L 25 -N 1 /home/shg047/db/hg19/bismark/ -1 ../fastq_trim/SRR949204_1_val_1.fq.gz -2 ../fastq_trim/SRR949204_2_val_2.fq.gz -o ../bam &


Memo for today's discussion: 
Human brain includes two broad classes of cells: neurons and glial cells. 
New understanding to AS based on deconvolution:
Region projection for AD based on region reference => which regions => Glial cells  or Neurons cells? 8224138
Glial cells  => which subtype ? => inflammatory related pathway? => eQTL => AD
Neurons cells =>  which subtype ? => which pathways ? => eQTL => AD
New Method based on deconvolution:
1,  Brain region (reference) and cell type (reference) projection for AD
2, Differential gene expression after deconvolution identify 
3, Causal network after deconvolution identify cell type (reference) specific interaction? 
Thanks. 

Shicheng


cp /media/Home_Raid1/zhl002/NAS3/WGBS/permutation/*bed  /opt/lampp/htdocs/


cd /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/hg19
for i in `ls *bed`
do
genome="hg19"
wget -O $i.results.tsv "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=$genome&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2F132.239.25.238%2Fshg047%2FNAS3%2FAlice%2FWGBS%2Fpermutation%2Fhg19%2F$i"
done

cd /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/mm9
for i in `ls *bed`
do
genome="mm9"
wget -O $i.results.tsv "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=$genome&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2F132.239.25.238%2Fshg047%2FNAS3%2FAlice%2FWGBS%2Fpermutation%2Fmm9%2F$i"
done


/opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/hg19

/opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/

报案：受到台湾公民（蘇勝慧）的骚扰，诽谤，损害名誉权，人身健康威胁及恐吓案件
chr14:29318659
fastq-dump --maxSpotId 100 --minSpotId 200 SRR1648428

samtools tview N22_bismark_bt2_pe.sort.bam ~/db/hg19/

qstat -u shg047 | grep condo | awk '{print $1}' | xargs -I {} qdel {}﻿
for i in 0 T U  
do 
rm *$i.bam
done

trim_galore --paired --phred33 --fastqc --illumina N21_R1.fastq.gz N21_R2.fastq.gz --output_dir ../fastq_trim
bismark --bowtie2 --multicore 4 --phred33-quals --fastq -L 21 -N 1 ~/db/hg19/align/bismark/ -1 ../fastq_trim/N21_R1_val_1.fq.gz -2 ../fastq_trim/N21_R2_val_2.fq.gz -o ../bam
filter_non_conversion --paired ../bam/N21_R1_val_1_bismark_bt2_pe.bam
samtools sort ../bam/N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam ../sortbam/N21_bismark_bt2_pe.sort
samtools index ../sortbam/N21_bismark_bt2_pe.sort.bam
bismark_methylation_extractor --cutoff 10 --paired-end --bedGraph --ignore 1 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam
bedtools intersect -wa -a N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov -b ../bed/target.bed > ../bedgraph/N21.bedgraph

N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.

chr8:67344678-67344679

cd /home/shg047/work/Minghua2016/methyfreq
for i in `ls *bismark.zero.cov`
do
bedtools intersect -wa -a $i -b ../target.bed > $i.bedgraph
done

bedtools intersect -wa -a T9_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov -b ../target.bed

trim_galore --paired --phred33 --fastqc --illumina SRR949193_1.fastq.gz SRR949193_2.fastq.gz --output_dir ../../fastq_tri
CU1_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam
CU1_val_1_bismark_bt2_pe.nonCG_filtered.bam
trim_galore --paired --phred33 --fastqc --illumina SRR949197_1.fastq.gz SRR949197_2.fastq.gz --output_dir ../fastq_trim
SCNT.wnt5a_mouse.sort.bam
SRR949205_1_val_1_bismark_bt2_pe.bam

/media/Home_Raid1/shg047/db/aligndb/hg19/bismark
scp SRR949206_1_trimmed.fq.gz shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/Ziller2013/fastq_trim 
scp SRR949206_2_trimmed.fq.gz shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/Ziller2013/fastq_trim
scp SRR949206.bismark.job shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/Ziller2013/fastq


zcat ../fastq_trim/SRR949205_1_val_1.fq.gz | head -n 4000 > ../SRR949205_1_val_1.fq
zcat ../fastq_trim/SRR949205_2_val_2.fq.gz | head -n 4000 > ../SRR949205_2_val_2.fq

gzip ../SRR949205_1_val_1.fq
gzip ../SRR949205_2_val_2.fq

bismark --bowtie2 --phred33-quals --fastq -L 25 -N 1 --multicore 6 /home/shg047/db/hg19/bismark/ -1 SRR949205_1_val_1.fq.gz -2 SRR949205_2_val_2.fq.gz -o ./

find . -exec touch {} \;


Dear Sir/Madam,
This user use my photo and my name as its profile and post kinds of things full of defamation and insult and post my private things through google plus to public. It break the California laws and the profiles are with full of defamation and insult. Please remove all the photos and my private things. Please warn the user it is illegal. 

https://plus.google.com/100087098265453308274
 
https://plus.google.com/u/0/108664216580872773217
 
https://plus.google.com/collection/UpxPlB?hl=en-US
 
https://plus.google.com/collection/04KQlB?hl=en-US
 
https://plus.google.com/117041654657062516807?hl=en-US
 
https://plus.google.com/100087098265453308274?hl=en-US

If you have any questions, please call me.  my phone number is: 281-685-5882

By the way, This is not the first time she do that,  I have noticed google plus week. google plus delete her profiles. However, she create a fake profile to impersonate me again. Please remove the profile as soon as possible. 

Thanks. 




find . -maxdepth 5 -type f -exec touch {} +  &


cp /home/shg047/work/Minghua2016/fastq/N53_R1.fastq.gz /home/shg047/work/Minghua2016/fastq/N53_R2.fastq.gz /home/shg047/software/Bismark/test

mv /home/shg047/work/DennisLo2015/sortBam/HPLC.bam.bai ./HPLC.bismark.bam.bai
mv /home/shg047/work/DennisLo2015/sortBam/HPLC.bam ./HPLC.bismark.bam


### Produce Run Time
my $end_run = time();
my $run_time = $end_run - $start_run;
my $days  = int($run_time/(24*60*60));
my $hours = ($run_time/(60*60))%24;
my $mins  = ($run_time/60)%60;
my $secs  = $run_time%60;

warn "Bismark completed in ${days}d ${hours}h ${mins}m ${secs}s\n";
print REPORT "Bismark completed in ${days}d ${hours}h ${mins}m ${secs}s\n";

chr14:95233733-95237069

samtools view -b CTR107_trimmed.fq.gz_bismark_bt2.sort.bam chr14:95233733-95237069

samtools view -H file.bam > header.sam  # extract header only
samtools reheader header.sam file.unique.bam


/home/shg047/work/Alice/WNT5A/mouse/bam

samtools view -H SCNT.wnt5a_mouse.sort.bam > header.txt
samtools view SCNT.wnt5a_mouse.sort.bam | tail -n 300 > test.sam
samtools view -b -T ~/db/mm9/mm9.fa test.sam -o test.bam
samtools sort test.bam -o test.sort.bam
samtools index test.sort.bam
samtools view test.sort.bam

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/test.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt 

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/test.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > ./SCNT.hapInfo.txt
/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/SCNT.wnt5a_mouse.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes MM9.CpG.positions.plus1.txt > SCNT.hapInfo.plus1.txt
/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/SCNT.wnt5a_mouse.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes MM9.CpG.positions.minus1.txt > SCNT.hapInfo.minus.txt

chr14:29321905

samtools tview SCNT.wnt5a_mouse.sort.bam ~/db/mm9/mm9.fa

Edit: An example script would be something of the form: 

samtools view -h foo.bam | awk '{if(found==0) {if($2=="SN:MT") {found=1; getline;}} print $0}' | samtools view -bSo foo.reformated.bam -


/home/shg047/oasis/db/mm9/MM9.CpG.positions.txt

samtools tview test.sort.bam ~/db/mm9/mm9.fa
samtools tview SCNT.wnt5a_mouse.sort.bam ~/db/mm9/mm9.fa


awk '{print $1,$2+1}' OFS="\t" /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > MM9.CpG.positions.plus1.txt
awk '{print $1,$2-1}' OFS="\t" /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > MM9.CpG.positions.minus1.txt



125194864
 
perl ~/bin/hapinfo2LDR2.pl ipsc.wnt5a chr14:29318659-29338701 < iPSC.hapInfo.txt
perl ~/bin/hapinfo2LDR2.pl scnt.wnt5a chr14:29318659-29338701 < SCNT.hapInfo.txt

bedtools coverage -d -abam SCNT.wnt5a_mouse.sort.bam -b wnt5a.bed > scnt.wnt5a.cov
bedtools coverage -d -abam iPSC.wnt5a_mouse.sort.bam -b wnt5a.bed > ipsc.wnt5a.cov

bismark_methylation_extractor --comprehensive -s --bedGraph SCNT.wnt5a_mouse.sort.bam -o /home/shg047/work/Alice/WNT5A/mouse/bam/rlt

PileOMeth extract -p 1 -q 1 --minDepth 1 --mergeContext ~/db/mm9/mm9.fa SCNT.wnt5a_mouse.sort.bam
 
 
 bam2fastq SCNT.wnt5a_mouse.sort.bam  
  bam2fastq iPSC.wnt5a_mouse.sort.bam 
 mv s_6_M_sequence.txt iPSC.fastq
 mv s_5_M_sequence.txt SCNT.fastq

 for i in iPSC.fastq SCNT.fastq
 do
 bismark  /home/shg047/oasis/db/mm9/ $i
 done
 
 
wc -l SCNT.wnt5a_mouse.sort_CpG.bedGraph

# sort the new bam files
for i in `ls *.WNT5A_human.bam`
do
# samtools sort $i -o $i.sort.bam
samtools index $i.sort.bam
done


# get the haplotype
perl ~/bin/bam2hapinfo

# calculate LDR


# get 



# get methylation level



cd /oasis/tscc/scratch/shg047/Alice/sortBam
/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/sortBam/human.bed /home/shg047/work/Alice/sortBam/SRR1286725_trimmed_bismark_bt2.sortc.bam.WNT5A_human.bam bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt > ../hapinfo/SRR1286725.hapInfo.txt

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/SCNT.wnt5a_mouse.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > ../hapinfo/SCNT.hapInfo.txt

sh SRR1286295.job &
sh SRR1286296.job &
sh SRR1286297.job &
sh SRR1286298.job &
sh SRR1286299.job &

Last time you mentioned: We should focus on two aspects: (i) the set of marker regions that you identified for mapping tissue-of-origin and tumors; (ii) the method/algorithm that you take to use these markers for making prediction.


# assign to the nearest promoter
bedtools sort -i /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.bed.txt > /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.sort.bed.txt
bedtools intersect -wo -a /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.sort.bed.txt -b /media/Home_Raid1/shg047/db/hg19/hg19_refGene.bed > /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.sort.bed.reference.txt


wget -r --user=zhang --password='HDIIWpP' ftp://igm-storage1.ucsd.edu/161121_D00611_0401_BH3J7CBCXY_Zhang10X/

bismark --bowtie2 --phred33-quals --fastq -L 25 -N 1 --multicore 6 /home/shg047/db/hg19/bismark/ -1 ../fastq_trim/T9_R1_trimmed.fq.gz -2 ../fastq_trim/T9_R2_trimmed.fq.gz -o ../bam --basename T9_R1

gem3-mapper -p --bisulfite-mode -I GCA_000001405.15_GRCh38_no_alt_analysis_set_BS.gem -s 1 -p -M 4

for i in `ls *.fastq`
do
gzip $i
done
 
perl ~/bin/pairAutoBismark.pl ../PRJNA201480.txt 33 non


RASSF1	/home/shg047/oasis/monod/Fig1B/bam/RASSF1.61-WGBS.bam	/home/shg047/oasis/monod/Fig1B/bam/mhb.bed
SHOX2	/home/shg047/oasis/monod/Fig1B/bam/SHOX2.61-WGBS.bam	/home/shg047/oasis/monod/Fig1B/bam/mhb.bed

perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

samtools sort RASSF1.61-WGBS.bam -o RASSF1.61-WGBS.sort.bam
samtools sort SHOX2.61-WGBS.bam -o SHOX2.61-WGBS.sort.bam

samtools index RASSF1.61-WGBS.sort.bam
samtools index SHOX2.61-WGBS.sort.bam

mv RASSF1.61-WGBS.sort.bam RASSF1.61-WGBS.bam
mv SHOX2.61-WGBS.sort.bam SHOX2.61-WGBS.bam
mv RASSF1.61-WGBS.sort.bam.bai RASSF1.61-WGBS.bam.bai
mv SHOX2.61-WGBS.sort.bam.bai SHOX2.61-WGBS.bam.bai

chr3	50373531	50384063	chr3:50373531-50384063	
chr3	157811261	157826490	chr3:157811261-157826490

RASSF1: chr3:50373531-50384063 
SHOX2: chr3:157811261-157826490

#!/bin/csh
#PBS -q pdafm
#PBS -l nodes=1:ppn=1
#PBS -l walltime=8:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd $PBS_O_WORKDIR
perl ~/bin/hapinfo2LDR2.pl SHOX2 chr3:157811261-157826490 < SHOX2.hapInfo.txt &
perl ~/bin/hapinfo2LDR2.pl RASSF1 chr3:50373531-50384063 < RASSF1.hapInfo.txt &

 
perl ~/bin/bam2hapInfo2PBS.pl mhb.bed submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

/media/Home_Raid1/shg047/monod/Fig1b/bam/mhb.bed

scp shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/monod/Fig1b/bam/* ./


shg047@genomeMiner:~/work/db/mm9$ /media/Home_Raid1/shg047/work/db/mm9/mm9.refGene.bed^C



samtools view -bS -t ~/db/hg19/hg19.fa.fai RASSF1.61-WGBS.sam -o ~/RASSF1.61-WGBS.bam
samtools view -bS -t ~/db/hg19/hg19.fa.fai SHOX2.61-WGBS.sam -o ~/SHOX2.61-WGBS.bam

mv /media/Home_Raid1/shg047/RASSF1.61-WGBS.bam ./
mv /media/Home_Raid1/shg047/SHOX2.61-WGBS.bam ./



zcat SRR1648425_1.fastq.gz | head -n 40000 > ./test/SRR1648425_1.fastq.gz
zcat SRR1648425_2.fastq.gz | head -n 40000 > ./test/SRR1648425_2.fastq.gz

for i in {1..5}
do
perl HMHDection.pl 6-P-$i.hapInfo.txt 6-T-$i.hapInfo.txt excl2.txt 6-$i.txt
perl HMHDection.pl 7-P-$i.hapInfo.txt 7-T-$i.hapInfo.txt excl2.txt 7-$i.txt
done

cd /home/shg047/work/monod/hapinfo/december/mix
cat ../6-*.txt.txt > CRC-MHM.txt
cat ../7-*.txt.txt > LC-MHM.txt

perl -lane 'print if ! /LOC/' CRC-MHM.txt > CRC-MHM2.txt
perl -lane 'print if ! /LOC/' LC-MHM.txt > LC-MHM2.txt

awk 'NR>1 {print $1}' CRC-MHM2.txt > CRC-MHM3.txt
awk 'NR>1 {print $1}' LC-MHM2.txt > LC-MHM3.txt

perl -p -i -e "s/[:-]/\t/g" CRC-MHM3.txt
perl -p -i -e "s/[:-]/\t/g" LC-MHM3.txt
 
sort -u CRC-MHM3.txt > CRC-MHM4.txt
sort -u LC-MHM3.txt > LC-MHM4.txt

#!/bin/csh
#PBS -q pdafm
#PBS -l nodes=1:ppn=1
#PBS -l walltime=8:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd $PBS_O_WORKDIR
# merge hapinfo for samples from same group
perl ~/bin/hapinfoMerge.pl NT.hapInfo.txt > NT.HapInfo.txt 
perl ~/bin/hapinfoMerge.pl CT.hapInfo.txt > CT.HapInfo.txt 
perl ~/bin/hapinfoMerge.pl H1.hapInfo.txt > H1.HapInfo.txt 
# then to do mhb calling
perl ~/bin/hapinfo2mhb.pl NT.hapInfo.txt.SumUniq 0.3 > NT.mhb &
perl ~/bin/hapinfo2mhb.pl H1.hapInfo.txt.SumUniq 0.3 > H1.mhb &
perl ~/bin/hapinfo2mhb.pl CT.hapInfo.txt.SumUniq 0.3 > CT.mhb &

perl ~/bin/hapinfo2mhb.pl W.hapInfo.txt.SumUniq 0.3 > W.mhb &

/media/Home_Raid1/shg047/work/monod/hapinfo/mix

cat H1.mhb NT.mhb > N.MHB
cat CT.mhb > C.MHB

bedtools sort -i N.MHB | wc -l
bedtools sort -i C.MHB | wc -l
bedtools intersect -wa -v -a C.MHB -b N.MHB | sort -u | wc -l 

bedtools intersect -wa -a C.MHB -b N.MHB | sort -u | wc -l 
bedtools intersect -wa -v -a C.MHB -b N.MHB | sort -u | wc -l 
bedtools intersect -wa -v -a N.MHB -b C.MHB | sort -u | wc -l 

cd /home/shg047/work/monod/hapinfo/June
load("monod.mhl.22July.RData")

NT.hapInfo.txt
-rwxrwxrwx 1 dinh sambashare  629 Nov 23 13:49 hapinfo2mhb.job
-rwxrwxrwx 1 dinh sambashare 255M Nov 23 13:49 CT.hapInfo.txt
-rwxrwxrwx 1 dinh sambashare    0 Nov 23 13:49 NT.mhb
-rwxrwxrwx 1 dinh sambashare 508M Nov 23 13:49 H1.hapInfo.txt

cd /home/shg047/work/Alice/sortBam
perl ~/bin/SaminfoPre4hapinfo.pl > sampleconfig.txt
perl ~/bin/bam2hapInfo2PBS.pl sampleconfig.txt non bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt


ls SRR1286395*
/home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed


#!/usr/bin/sh
START=$(date +%s)
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "It takes DIFF=$DIFF seconds to complete this task..."

ls *gz | awk -F. '{print $1}'
 
 diff <(grep "run complete" *.err | awk -F: '{print $1}'|sort) <(grep kill *.err | awk -F: '{print $1}' | sort)
 
grep "run complete" *.err | awk -F: '{print $1}'|sort > a
ls ../bam/*bam | grep -v temp | awk -F[/_] '{print $3".err"}' | sort | wc -l> b 
diff <(grep "run complete" *.err | awk -F: '{print $1}'|sort) <(ls ../bam/*bam | grep -v temp | awk -F[/_] '{print $3".err"}' )
paste a b > c


SRR1286589.pbs



xx_trimmed_bismark_bt2.bam

# bash date
START=$(date +%s)
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "It takes DIFF=$DIFF seconds to complete this task..."﻿

# perl time
my $start_run = time();
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";

Note: setting of bismark alignment in PBS system (DNA methylation, BS-seq, RRBS)

multicore=1, ppn=8, memory=3.2G*6
multicore=2, ppn=16, memory=3.2G*12

note: '-p 1' will already use 
* 4 threads/cores for Bowtie2 (3.2G/per) plus 1 additional core for Bismark itself
* 2 threads/cores for Perl (3.2G/per)
* 2 threads/cores for SAMTOOLS

note: '-p 2' will already use 
* 8 threads/cores for Bowtie2 (3.2G/per) plus 1 additional core for Bismark itself
* 4 threads/cores for Perl (3.2G/per)
* 4 threads/cores for SAMTOOLS

62067 shg047    20   0 3353m 3.1g  812 R 100.0  5.0  32:53.22 perl
61902 shg047    20   0 3353m 3.1g 2168 R 100.0  5.0  33:48.64 perl
62048 shg047    20   0 3353m 3.1g 2168 R 100.0  5.0  33:29.31 perl
62126 shg047    20   0 3353m 3.1g  812 R 99.7  5.0  32:47.26 perl
62155 shg047    20   0 3357m 3.2g 1968 S 59.3  5.0  18:00.17 bowtie2-align-s
62154 shg047    20   0 3357m 3.2g 1968 S 58.9  5.0  17:55.70 bowtie2-align-s
62169 shg047    20   0 3357m 3.2g 1968 S 58.6  5.0  17:46.12 bowtie2-align-s
62146 shg047    20   0 3357m 3.2g 1968 S 58.3  5.0  17:55.60 bowtie2-align-s
62167 shg047    20   0 3357m 3.2g 1968 S 57.6  5.0  17:45.36 bowtie2-align-s
62171 shg047    20   0 3357m 3.2g 1968 S 57.6  5.0  17:50.46 bowtie2-align-s
62149 shg047    20   0 3357m 3.2g 1968 S 55.6  5.0  18:00.27 bowtie2-align-s
62168 shg047    20   0 3357m 3.2g 1968 S 55.3  5.0  17:49.37 bowtie2-align-s
62178 shg047    20   0 22884 5388  780 S 13.9  0.0   4:07.75 samtools
62183 shg047    20   0 22884 5388  780 S 13.6  0.0   4:08.99 samtools
62180 shg047    20   0 22884 5388  780 S 13.2  0.0   4:08.05 samtools
62184 shg047    20   0 22884 5388  780 S 13.2  0.0   4:08.89 samtools



/home/shg047/oasis/db/hg19/align/bismark

/media/Home_Raid1/shg047/git/Bismark/bismark ~/db/hg19/align/bismark/ --multicore=4 --non_directional -1 N36_R1.fastq -2 N36_R2.fastq

perl /media/Home_Raid1/shg047/git/Bismark/


# How to install git with adminstration 
sudo apt-get update
sudo apt-get install git
git config --global user.name "Shicheng-Guo"
git config --global user.email "shicheng.guo@hotmail.com"
git config --list
git commit --amend --reset-author

# How to install git without admin (TSCC in UCSD)
cd software
wget https://www.kernel.org/pub/software/scm/git/git-2.10.2.tar.xz
tar xf xzvf git-2.10.2.tar.xz
cd git-2.10.2
./configure --prefix="$HOME/git"
make install
git config --list

# Go to home and bulid the databse for git
cd
mkidr git



	
cd 
scp /media/Home_Raid1/shg047/NAS3/Alice/human/fastq
scp ../human/SRR1286624.fastq.gz.bismark.pbs shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS3/Alice/human/fastq

ls *fastq | echo

ls *fastq | paste - -




wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/pysam/pysam-0.6.tar.gz

cd /media/Home_Raid1/shg047/NAS3/Minghua2016/fastq
bedtools unionbedg -i *bedgraph -filler NA -header -names CT1.bam.bedgraph CT2.bam.bedgraph CT3.bam.bedgraph CT.bam.bedgraph CU1.bam.bedgraph CU2.bam.bedgraph CU3.bam.bedgraph CU.bam.bedgraph ET1.bam.bedgraph ET2.bam.bedgraph ET3.bam.bedgraph ET.bam.bedgraph EU1.bam.bedgraph EU2.bam.bedgraph EU3.bam.bedgraph EU.bam.bedgraph HE.bam.bedgraph N100.bam.bedgraph N101.bam.bedgraph N102.bam.bedgraph N103.bam.bedgraph N10.bam.bedgraph N11.bam.bedgraph N12.bam.bedgraph N13.bam.bedgraph N14.bam.bedgraph N15.bam.bedgraph N16.bam.bedgraph N17.bam.bedgraph N18.bam.bedgraph N19.bam.bedgraph N1.bam.bedgraph N20.bam.bedgraph N21.bam.bedgraph N22.bam.bedgraph N23.bam.bedgraph N24.bam.bedgraph N25.bam.bedgraph N26.bam.bedgraph N27.bam.bedgraph N28.bam.bedgraph N29.bam.bedgraph N2.bam.bedgraph N30.bam.bedgraph N31.bam.bedgraph N32.bam.bedgraph N33.bam.bedgraph N34.bam.bedgraph N35.bam.bedgraph N36.bam.bedgraph N37.bam.bedgraph N38.bam.bedgraph N39.bam.bedgraph N3.bam.bedgraph N40.bam.bedgraph N41.bam.bedgraph N42.bam.bedgraph N43.bam.bedgraph N44.bam.bedgraph N45.bam.bedgraph N46.bam.bedgraph N47.bam.bedgraph N48.bam.bedgraph N49.bam.bedgraph N4.bam.bedgraph N50.bam.bedgraph N51.bam.bedgraph N52.bam.bedgraph N53.bam.bedgraph N54.bam.bedgraph N55.bam.bedgraph N56.bam.bedgraph N57.bam.bedgraph N58.bam.bedgraph N59.bam.bedgraph N5.bam.bedgraph N60.bam.bedgraph N61.bam.bedgraph N62.bam.bedgraph N63.bam.bedgraph N64.bam.bedgraph N65.bam.bedgraph N66.bam.bedgraph N67.bam.bedgraph N68.bam.bedgraph N69.bam.bedgraph N6.bam.bedgraph N70.bam.bedgraph N71.bam.bedgraph N72.bam.bedgraph N73.bam.bedgraph N74.bam.bedgraph N75.bam.bedgraph N76.bam.bedgraph N77.bam.bedgraph N78.bam.bedgraph N79.bam.bedgraph N7.bam.bedgraph N80.bam.bedgraph N81.bam.bedgraph N82.bam.bedgraph N83.bam.bedgraph N84.bam.bedgraph N85.bam.bedgraph N86.bam.bedgraph N87.bam.bedgraph N88.bam.bedgraph N89.bam.bedgraph N8.bam.bedgraph N90.bam.bedgraph N91.bam.bedgraph N92.bam.bedgraph N93.bam.bedgraph N94.bam.bedgraph N95.bam.bedgraph N96.bam.bedgraph N97.bam.bedgraph N98.bam.bedgraph N99.bam.bedgraph N9.bam.bedgraph T100.bam.bedgraph T101.bam.bedgraph T102.bam.bedgraph T103.bam.bedgraph T10.bam.bedgraph T11.bam.bedgraph T12.bam.bedgraph T13.bam.bedgraph T14.bam.bedgraph T15.bam.bedgraph T16.bam.bedgraph T17.bam.bedgraph T18.bam.bedgraph T19.bam.bedgraph T1.bam.bedgraph T20.bam.bedgraph T21.bam.bedgraph T22.bam.bedgraph T23.bam.bedgraph T24.bam.bedgraph T25.bam.bedgraph T26.bam.bedgraph T27.bam.bedgraph T28.bam.bedgraph T29.bam.bedgraph T2.bam.bedgraph T30.bam.bedgraph T31.bam.bedgraph T32.bam.bedgraph T33.bam.bedgraph T34.bam.bedgraph T35.bam.bedgraph T36.bam.bedgraph T37.bam.bedgraph T38.bam.bedgraph T39.bam.bedgraph T3.bam.bedgraph T40.bam.bedgraph T41.bam.bedgraph T42.bam.bedgraph T43.bam.bedgraph T44.bam.bedgraph T45.bam.bedgraph T46.bam.bedgraph T47.bam.bedgraph T48.bam.bedgraph T49.bam.bedgraph T4.bam.bedgraph T50.bam.bedgraph T51.bam.bedgraph T52.bam.bedgraph T53.bam.bedgraph T54.bam.bedgraph T55.bam.bedgraph T56.bam.bedgraph T57.bam.bedgraph T58.bam.bedgraph T59.bam.bedgraph T5.bam.bedgraph T60.bam.bedgraph T61.bam.bedgraph T62.bam.bedgraph T63.bam.bedgraph T64.bam.bedgraph T65.bam.bedgraph T66.bam.bedgraph T67.bam.bedgraph T68.bam.bedgraph T69.bam.bedgraph T70.bam.bedgraph T71.bam.bedgraph T72.bam.bedgraph T73.bam.bedgraph T74.bam.bedgraph T75.bam.bedgraph T76.bam.bedgraph T77.bam.bedgraph T78.bam.bedgraph T79.bam.bedgraph T7.bam.bedgraph T80.bam.bedgraph T81.bam.bedgraph T82.bam.bedgraph T83.bam.bedgraph T84.bam.bedgraph T85.bam.bedgraph T86.bam.bedgraph T87.bam.bedgraph T88.bam.bedgraph T89.bam.bedgraph T8.bam.bedgraph T90.bam.bedgraph T91.bam.bedgraph T92.bam.bedgraph T93.bam.bedgraph T94.bam.bedgraph T95.bam.bedgraph T96.bam.bedgraph T97.bam.bedgraph T98.bam.bedgraph T99.bam.bedgraph T9.bam.bedgraph > esca.bg
python ~/software/BSseeker2/bs_seeker2-align.py -1 N36_R1.fastq -2 N36_R2.fastq --aligner=bowtie2 -o ./96.bam -f bam -g /media/Home_Raid1/shg047/NAS3/Minghua2016/fa/target.fa
python ~/software/BSseeker2/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i 96.bam -o 96 --txt --db /media/Home_Raid1/shg047/software/BSseeker2/bs_utils/reference_genomes/target.fa_bowtie2


bedtools interact /home/shg047/db/hg19/lncRNA.hg19.bed
132.239.25.238

awk '{print $1,$3,$4,$10}' -OFS="\t" gencode.v19.long_noncoding_RNAs.gtf | head

bismark -q --phred33-quals --multicore=1 -n 1 -l 20 --non_directional /media/Home_Raid1/shg047/db/aligndb/hg19/bismark N36.R1.fastq -o ../bam/

chr1	35037044	35037796	chr1:35037044-35037796

Sylvia Schmalz 

for i in `ls *bam`
do

python /media/Home_Raid1/shg047/software/BSseeker2/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i $i -o $i --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
done




ls -larth SRR1286295.fastq.gz
ls -larth SRR1286296.fastq.gz 
ls -larth SRR1286297.fastq.gz
ls -larth SRR1286298.fastq.gz
ls -larth SRR1286299.fastq.gz


ENSG00000270170.1

GSM978967
cat cerebellum_MethylC-seq_chr*.BED > GSM978968.Brain.BED

for i in {1..11}
do
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw
done
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas1/tracks_hg19/Human_NormalPancreas1.meth.bw
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas2/tracks_hg19/Human_NormalPancreas2.meth.bw
wget http://smithlab.usc.edu/methbase/data/Gao-Human-2015/Human_BloodHealthy/tracks_hg19/Human_BloodHealthy.meth.bw

Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex






/media/Home_Raid1/shg047/NAS3/tcga/chol/lncrna/bam
cd /home/shg047/oasis/Minghua2016/fastq
bismark -q --phred33-quals -n 1 -l 20 --non_directional /home/shg047/oasis/db/hg19 -1 CT1_R1.fastq -2 CT1_R2.fastq -o ../

bismark -q --phred33-quals --multicore=8 -n 1 -l 20 --non_directional /media/Home_Raid1/shg047/db/aligndb/hg19/bismark -1 T96.R1.fastq -2 T96.R2.fastq -o ./
bismark -q --phred33-quals --multicore=8 -n 1 -l 20 --non_directional /media/Home_Raid1/shg047/db/aligndb/hg19/bismark -1 N36.R1.fastq -2 N36.R2.fastq -o ./



wget http://pellegrini.mcdb.ucla.edu/BS_Seeker2/sup_data/DS2_PE_simu_perfect_end1.fa.gz
wget http://pellegrini.mcdb.ucla.edu/BS_Seeker2/sup_data/DS2_PE_simu_perfect_end2.fa.gz

pip install 'pysam=6,<7'

cd /home/shg047/oasis/Minghua2016/fastq
python /home/shg047/software/BSseeker2-master/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i N46.bam -o N46.mr --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
bs_seeker2


cd /home/shg047/db/hg19/
REF="/home/shg047/db/hg19/hg19.fa";
python /home/shg047/software/bwa-meth/bwameth.py index $REF
REF="/home/shg047/db/hg19/hg19.fa";
python /home/shg047/software/bwa-meth/bwameth.py --reference $REF N28_R1.fastq N28_R2.fastq | samtools view -b - > N28.bwameth.bam
samtools view -cf2 N28.bwameth.bam
samtools view -cF4 N28.bwameth.bam
samtools flagstat N28.bwameth.bam

/home/shg047/software/Python-2.7.4/python /home/shg047/software/BSseeker2_v2.0.9/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i N46.bam -o N46.mr --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
bs_seeker2


	 samtools index T19.bam
	 samtools sort T19.bam -o T19.sort
	 

wget ftp://ftp.gnu.org/pub/gnu/emacs/emacs-25.1.tar.xzvf
tar xf emacs-25.1.tar.xzvf
./configure
make
export PATH=$PAHT:

cd /home/shg047/oasis/Minghua2016/fastq

#!/usr/bin/sh
mkdir $HOME/db
cd $HOME/db
for i in hg18 hg19 hg38 mm9 mm10 
do
[! -d "$i"] && mkdir $i 
cd $i
mkdir "fa"
cd "fa"
wget -r -np http://hgdownload.soe.ucsc.edu/goldenPath/$i/bigZips/$i.chromFa.tar.gz  
wget -r -np http://hgdownload.soe.ucsc.edu/goldenPath/$i/bigZips/chromFa.tar.gz     
wget -r -np http://hgdownload.soe.ucsc.edu/goldenPath/$i/bigZips/$i.chrom.sizes
tar xzvf *.gz
rm *.gz
cat *fa ../$i.fa
cd ../../
done

wget -r http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/$i.chromFa.tar.gz  
wget -r http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz  
tar xzvf *.gz
rm *.gz
cat *fa ../$i.fa


python ~/software/BSseeker2/bs_seeker2-build.py -f /home/shg047/db/hg19/hg19.fa --aligner=bowtie2
/home/shg047/software/Python-2.7.12/python /home/shg047/software/BSseeker2/bs_seeker2-align.py -1 CT1_R1.fastq -2 CT1_R2.fastq -g /home/shg047/oasis/db/hg19/hg19.fa -t Y
cd 

python ~/software/BSseeker2/bs_seeker2-build.py -f /home/puweilin/Software/bowtie2-2.2.9/ChrM/Homo_sapiens.GRCh38.dna.chromosome.MT.fa --aligner=bowtie2

python ~/software/BSseeker2/bs_seeker2-align.py -1 T86_R1.fastq -2 T86_R2.fastq --aligner=bowtie2 -o ./T86.bam -f bam -g /home/shg047/db/hg19/hg19.fa

# genome-miner
/home/shg047/software/Python-2.7.4/python ~/software/BSseeker2/bs_seeker2-build.py -f /home/shg047/db/hg19/hg19.fa --aligner=bowtie2


Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex





bsrate -c ~/oasis/db/hg19/hg19.fa -o T84.mr.bsrate T84.mr

# Line-1 Methylation Primer Set-1
5′-GGACGTATTTGGAAAATCGGG-3′ 
5′-AATCTCGCGATACGCCGTT-3′ 
5′-TCGAATATTGCGTTTTCGGATCGGTTT-3′ 

# Line-1 Methylation Primer Set-2
5′-TTGAGTTGTGGTGGGTTTTATTTAG-3′ 
5′-TCATCTCACTAAAAAATACCAAACA-3′.

# Line-1 Methylation Primer Set-3 
5’-TTGGTTAGGTGTGGGATATAGTT-3’
5’-CAAAAAATCAAAAAATTCCCTTTCC-3’

mkdir pip
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py --prefix="./"
pip install --user --upgrade cutadapt

~/.local/bin/cutadapt --help
cp ~/.local/bin/cutadapt ~/bin/
cp ~/.local/bin/cutadapt /media/Home_Raid1/shg047/software/trim_galore_zip

git clone https://github.com/BSSeeker/BSseeker2.git

git clone https://github.com/marcelm/cutadapt.git
bismark -q --phred33-quals -n 1 -l 40 --non_directional ~/NAS3/
bismark -q --phred33-quals --bowtie2 -N 1 -L 30 ~/NAS3/db/hg19/ -1 T84_R1.fastq -2 T84_R2.fastq
bismark -q --phred33-quals -n 1 -l 30 ~/NAS3/db/hg19/ -1 T84_R1_val_1.fq -2 T84_R2_val_2.fq
bismark --bowtie2 --phred33-quals --RRBS -- --fastq -L 30 -N 1 --multicore 4 ~/NAS3/db/hg19/ -1 T84_R1_val_1.fq -2 T84_R2_val_2.fq  -o ./

trim_galore --phred33 --paired --fastqc --non_directional --rrbs --illumina -o ./ T84_R1.fastq T84_R2.fastq
trim_galore --paired -a GATCGGAAGAGCA -a2 GCTCTTCCGATCT --retain_unpaired  --trim1  T21.5.read1.fq.head T21.5.read2.fq.head  

chr3:142839987
chr17:59534637
chr2:74742664
chr19:58446159-58447359
19	58446371
chr12:95942887

walt -i ~/NAS3/db/hg19/chrosome/methpipe.hg19.dbindex -1 T84_R1.fastq -2 T84_R2.fastq -o T84.mr



for i in `ls *gz`
do
gunzip -f $i
done

bismark --bowtie2 --phred33-quals --fastq -L 30 -N 1 --multicore 4 ~/NAS3/db/hg19/ -1 T83_R1.fastq.gz -2 T83_R2.fastq.gz  -o ../bam
chr10:118810290
chr13:103498159
chr19:37341777
chr19:12163473
chr10:123922851
chr12:128751882

makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex

# 2016-11-03
cat cerebellum_MethylC-seq_chr*.BED > GSM978968.Brain.BED

for i in {1..11}
do
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw
done
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas1/tracks_hg19/Human_NormalPancreas1.meth.bw
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas2/tracks_hg19/Human_NormalPancreas2.meth.bw
wget http://smithlab.usc.edu/methbase/data/Gao-Human-2015/Human_BloodHealthy/tracks_hg19/Human_BloodHealthy.meth.bw

# 2016-11-03
Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex



for i in `ls *bw`
do
bigWigAverageOverBed $i ES.bed $i.es.out
done

# download histone modification from roadmap/encode
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29611/
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56712
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56712/suppl/GSE56712_RAW.tar
tar xvf GSE56712_RAW.tar
rm GSE56712_RAW.tar
gunzip *.gz
rm *bam
rm *bigWig

# In terms of Super-Enhancer(H3k27ac)
cat *H3k27ac*broadPeak > Enhancer.txt
awk '{if ($7>5 && $5>100) print $1,$2,$3}' OFS="\t" Enhancer.txt > Enhancer.bed
sort -k1,1 -k2,2n Enhancer.bed > EnhancerSort.bed
bedtools merge -i EnhancerSort.bed > Human-Enhancer-hg19-H3k27ac-PMID24119843.bed
mv Human-Enhancer-hg19-PMID24119843.bed Human-Enhancer-hg19-H3k27ac-PMID24119843.bed
wc -l Human-Enhancer-hg19-H3k27ac-PMID24119843.bed

perl bed3Tobed4.pl Human-hg19-H3k27ac-PMID24119843.bed > Human-hg19-H3k27ac-PMID24119843.bed4

bedtools intersect -u -wa -a ~/NAS3/db/GPL13534.map -b superEnhancerCluster.bed | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/GPL13534.map | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/hg19/CpGI.hg19.bed | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/hg19/CpG.Shore.hg19.bed | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/hg19/CpG.Shelf.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shore.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpGI.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpGI.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shore.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shelf.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/Enhancers.Fantom.hg19.bed -b superEnhancer.bed | wc -l

# In terms of Promoter(H3k4me3)
zcat *H3k4me3*broadPeak.gz > promter.txt
awk '{if ($7>5) print $1,$2,$3}' OFS="\t" promter.txt > promter.bed
sort -k1,1 -k2,2n promter.bed > promterSort.bed
bedtools merge -i promterSort.bed > promterSortMerge.bed
wc -l promterSortMerge.bed
bedtools intersect -u -wa -a ~/NAS3/db/GPL13534.map -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpGI.hg19.bed -b promterSortMerge.bed| wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shore.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpGI.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpGI.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shore.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shelf.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/Enhancers.Fantom.hg19.bed -b superEnhancer.bed | wc -l

TCGA-ESCA-PBMC.bed
bedtools intersect -u -wa -a TCGA-ESCA-PBMC.bed -b ~/NAS3/db/hg19/ 




ls -larth *H3k27ac*broadPeak | wc -l
ls -larth *H3k27me3*broadPeak | wc -l
ls -larth *H3k79me2*broadPeak | wc -l
ls -larth *H3k9me3*broadPeak | wc -l
ls -larth *H3k4me2*broadPeak | wc -l

cd 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56712/suppl/GSE56712_RAW.tar
tar xvf GSE56712_RAW.tar
rm GSE56712_RAW.tar
rm *bam
rm *bigWig
gunzip *.gz
for i in H3k4me1 H3k4me2 H3k4me3 H3k9ac H3k9me1 H3k9me3 H3k27ac H3k27me3 H3k36me3 H3k79me2 H4k20me1 
do
cat *$i*broadPeak > promter.txt
awk '{if ($7>5) print $1,$2,$3}' OFS="\t" promter.txt > promter.bed
sort -k1,1 -k2,2n promter.bed > promterSort.bed
bedtools merge -i promterSort.bed > Human-hg19-$i-PMID24119843.bed
done
rm promter.txt
rm promter.bed
rm promterSort.bed


TCGA-ESCA-PBMC.bed

cd /media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca



#!/bin/csh
#PBS -q hotel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=168:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shicheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /oasis/tscc/scratch/shg047/Roadmap/wig
bigWigCorrelate bwList.txt > bwListCorrelationResult.txt


for i in chr10:79857968-79858049 chr10:79858088-79858205 chr10:79858211-79858298 chr10:79858319-79858347 chr6:122658641-122658674 chr6:122659252-122659297
do
perl ~/bin/hapinfo2LDR2.pl $i.tissue.R2.txt $i < tissue.hapinfo.txt 
perl ~/bin/hapinfo2LDR2.pl $i.ips.R2.txt $i < ipsNT.hapinfo.txt
perl ~/bin/hapinfo2LDR2.pl $i.esc.R2.txt $i < ESC.hapinfo.txt
done

file<-list.files(pattern="*.rsq$")
for(i in file){
print(i)
M<-read.table(i,head=T,row.names=1,as.is=T)
library("grDevices")
col=colorRampPalette(c("white", "red"))(10) 
M[lower.tri(M)] <- NA
pdf(paste(i,"pdf",sep="."))
image(data.matrix(M),col = col,frame=F,xaxt="n",yaxt="n")
dev.off()
}

#!/bin/csh
#PBS -q hotel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=9:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shicheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /oasis/tscc/scratch/shg047/Roadmap/wig
wigToBigWig GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz ~/oasis/db/hg19/hg19.chrom.sizes GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.bw


/home/shg047/oasis/Holger2016/bedGraph

cd ~/bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/addCols                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ameme                   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/autoDtd                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/autoSql                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/autoXml                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ave                     
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/aveCols                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtChain                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtSort                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtSwap                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtToMaf                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtToPsl                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedCommonRegions        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedCoverage             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedExtendRanges         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGeneParts            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphPack            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedIntersect            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedItemOverlapCount     
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedPileUps              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedRemoveOverlap        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedRestrictToPositions  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToExons              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToPsl                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedWeedOverlapping      
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedInfo              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedNamedItems        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedSummary           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigPslToPsl             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed    
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCluster           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCorrelate         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigInfo              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blastToPsl              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blastXmlToPsl           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/                   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/calc                    
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/catDir                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/catUncomment            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainAntiRepeat         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainFilter             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainNet                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainPreNet             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSort               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSplit              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainStitchId           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSwap               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainToAxt              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainToPsl              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainToPslBasic         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/checkAgpAndFa           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/checkCoverageGaps       
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/checkHgFindSpec         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/checkTableCoords        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chopFaLines             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chromGraphFromBin       
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chromGraphToBin         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/colTransform            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/countChars              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/crTreeIndexBed          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/crTreeSearchBed         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/dbSnoop                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/dbTrash                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/estOrient               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faAlign                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faCmp                   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faCount                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faFilter                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faFilterN               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faFrag                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faNoise                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faOneRecord             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faPolyASizes            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faRandomize             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faRc                    
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToFastq               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTab                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faTrans                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fastqStatsAndSubsample  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fastqToFa               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/featureBits             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/findMotif               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gapToLift               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredCheck           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredFilter          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredHisto           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredSingleCover     
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBigGenePred   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToFakePsl       
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToMafFrames     
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToProt          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/getRna                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/getRnaPred              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToPsl               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gmtime                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/headRest                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgFindSpec              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgGcPercent             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadBed               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadChain             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadMaf               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadMafSummary        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadNet               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadOut               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadOutJoined         
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgLoadWiggle            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgSpeciesRna            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgTrackDb               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgWiggle                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgsql                   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hgsqldump               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/htmlCheck               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hubCheck                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hubPublicCheck          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ixIxx                   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/lavToAxt                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/lavToPsl                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ldHgGene                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOverMerge           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftUp                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/linesToRa               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/localtime               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafAddIRows             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafAddQRows             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafCoverage             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafFetch                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafFilter               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafFrag                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafFrags                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafGene                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafMeFirst              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafOrder                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafRanges               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafSpeciesList          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafSpeciesSubset        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafSplit                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafSplitPos             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafToAxt                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafToPsl                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafToSnpBed             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafsInRegion            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/makeTableList           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/maskOutFa               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mktime                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mrnaToGene              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netChainSubset          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netClass                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netFilter               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netSplit                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netSyntenic             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netToAxt                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netToBed                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/newProg                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/newPythonProg           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/nibFrag                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/nibSize                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/oligoMatch              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/overlapSelect           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/paraFetch               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/paraSync                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/positionalTblCheck      
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCat                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCheck                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslDropOverlap          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslFilter               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslHisto                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslLiftSubrangeBlat     
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslMap                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslMrnaCover            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslPairs                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslPartition            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslPosTarget            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslPretty               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslRecalcMatch          
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslReps                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslScore                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslSelect               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslSort                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslStats                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslSwap                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslToBed                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslToBigPsl             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslToChain              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslToPslx               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslxToFa                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/qaToQac                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/qacAgpLift              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/qacToQa                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/qacToWig                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/raSqlQuery              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/raToLines               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/raToTab                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/randomLines             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/rmFaDups                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/rowsToCols              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/sizeof                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/spacedToTab             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/splitFile               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/splitFileByColumn       
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/sqlToXml                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/stringify               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/subChar                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/subColumn               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/tailLines               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/tdbQuery                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/textHistogram           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/tickToDate              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/toLower                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/toUpper                 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/transMapPslToGenePred   
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/trfBig                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitDup               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitMask              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa              
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/udr                     
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/validateFiles           
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/validateManifest        
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate            
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigEncode               
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig             
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wordLine                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/xmlCat                  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/xmlToSql                
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/


for i in `ls GSM*`
do
perl GSE52270.pl $i > $i.bedgraph &
done

for i in `ls *bedgraph`
do
sort -k1,1 -k2,2n $i > $i.sort &
done

for i in `ls *bedgraph`
do
bedGraphToBigWig $i.sort ~/oasis/db/hg19/hg19.chrom.sizes $i.sort.hg19.bw &
done

for i in `ls *.wig.gz`
do
wigToBigWig $i ~/oasis/db/hg19/hg19.chrom.sizes $i.bw &
done





for i in `ls *hapinfo.txt`
do
perl ~/bin/hapinfo2BlocAvgR2.pl $i > R2.$i
done
perl ../

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl tissue.hapinfo.txt > tissue.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl tissue.hapinfo.txt.$i > R2.tissue.hapinfo.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl tissue.hapinfo.txt > tissue.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl tissue.hapinfo.txt.$i > R2.tissue.hapinfo.txt.$i
perl ~/bin/randomSampleFromHaploInfo.pl ipsNT.hapinfo.txt > ipsNT.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ipsNT.hapinfo.txt.$i > R2.ipsNT.hapinfo.txt.$i
perl ~/bin/randomSampleFromHaploInfo.pl ESC.hapinfo.txt > ESC.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ESC.hapinfo.txt.$i > R2.ESC.hapinfo.txt.$i
done

cp /oasis/tscc/scratch/zhl002/WGBS_mouse/HapInfo/ESC.hapinfo.txt.SumUniq  /home/shg047/oasis/mouse/hapinfo/group/ESC.hapinfo.txt 
cp /oasis/tscc/scratch/zhl002/WGBS_mouse/HapInfo/ipsNT.hapinfo.txt.SumUniq  /home/shg047/oasis/mouse/hapinfo/group/ipsNT.hapinfo.txt 
cp /oasis/tscc/scratch/zhl002/WGBS_mouse/HapInfo/tissue.hapinfo.txt.SumUniq  /home/shg047/oasis/mouse/hapinfo/group/tissue.hapinfo.txt

perl 


cd /home/shg047/oasis/mouse/sortBam
perl ../saminfoPre4bam2hapinfo.pl > Samconfig2.txt
perl /home/shg047/bin/bam2hapInfo2PBS.pl Samconfig2.txt submit bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt

for i in `ls *fastq`
do
zip $i.gz $i
done


perl ../saminfoPre4bam2hapinfo.pl > SamConfig.txt

SRX209451
SRX271137
SRX271142


 #!/bin/csh
 #PBS -N mf2bedGraph
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 perl ~/bin/hapinfo2LDR2.pl  rlt chr10:110930425-110930618 Mouse.Merge.Hapinfo.txt.SumUniq
 

 #!/bin/csh
 #PBS -N mf2bedGraph
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 perl ~/bin/hapinfo2mhb.pl Mouse.Merge.Hapinfo.txt.SumUniq 0.5 > MHB.RD90up80.r0.5-2.bed

 
 cd /home/shg047/oasis/mouse/hapinfo
perl /home/shg047/oasis/mouse/mergeHapinfo/hapinfoMerge.pl Mouse.Merge.Hapinfo.txt Mouse.Merge.Hapinfo.Merge.txt

 #!/bin/csh
 #PBS -N mf2bedGraph
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 perl /home/shg047/oasis/mouse/mergeHapinfo/hapinfoMerge.pl Mouse.Merge.Hapinfo.txt Mouse.Merge.Hapinfo.Merge.txt
 
 #!/bin/csh
 #PBS -N hapinfo2mhl
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
perl ~/bin/MethylFreq2wig.pl 
 
 
cd /oasis/tscc/scratch/zhl002/WGBS_mouse/bam
chr1:3000017-3000744
   
   view -b Indx01.merged.bam.sorted.bam chr1:3000017-3000744 > xx.bam
   
   
qsub SRX080192.job
qsub SRX1091397.job
qsub SRX209448.job
qsub SRX209453.job
qsub SRX209458.job
qsub SRX271136.job



/home/shg047/oasis/mouse/mergeHapinfo/hapinfo/MHB.Mouse.R2-0.5.bed



for i in `ls *job`
do
qsub $i
done

rm SRX1019864.job
rm SRX1019865.job
rm SRX1019866.job
rm SRX1019867.job

for i in {4..7}
do
qsub SRX101986$i.job
done




samtoos view SRR630257_1_trimmed_bismark_bt2.bam
SRR630258_1_trimmed_bismark_bt2.bam

SRR630257.12464425_HWI-ST1113_0085:5:1308:18692:107279_length=100       16      chr4    2999991 8       95M     *       0       0       TCACCAACAAAAATTCACTCGAACAACTAAATTCTTCTTTTTTTTTTTCCATT
SRR630257.5426488_HWI-ST1113_0085:5:1203:8069:184455_length=100 0       chr4    3001064 42      95M     *       0       0       GATTTTTAGAGTGGTTGTATAAGTTTAGAATTTTATTAATAATGGAAGAGTGTTTTTTTTT

SRR630258.15688879_HWI-ST1113_0085:6:2106:2276:6097_length=100  16      chr4    3000069 39      95M     *       0       0       CATTTCCAATACTATACCAAAAATCCCCCATACCTACCCACCCGCACTCCCCTACCCACCT
SRR630258.8741011_HWI-ST1113_0085:6:1208:15012:167343_length=100        0       chr4    3000448 32      95M     *       0       0       TTCGATAAAATTTTGTTAGTGTATGTAATGGTGTTAGCGTTTGGATGTTGATT
(
SRR630257.12464425_HWI-ST1113_0085:5:1308:18692:107279_length=100       16      chr4    2999991 8       95M     *       0       0       TCACCAACAAAAATTCACTCGAACAACTAAATTCTTCTTTTTTTTTTTCCATT
SRR630258.15688879_HWI-ST1113_0085:6:2106:2276:6097_length=100  16      chr4    3000069 39      95M     *       0       0       CATTTCCAATACTATACCAAAAATCCCCCATACCTACCCACCCGCACTCCCCTACCCACCT
SRR630258.8741011_HWI-ST1113_0085:6:1208:15012:167343_length=100        0       chr4    3000448 32      95M     *       0       0       TTCGATAAAATTTTGTTAGTGTATGTAATGGTGTTAGCGTTTGGATGTTGATT
SRR630257.5426488_HWI-ST1113_0085:5:1203:8069:184455_length=100 0       chr4    3001064 42      95M     *       0       0       GATTTTTAGAGTGGTTGTATAAGTTTAGAATTTTATTAATAATGGAAGAGTGTTTTTTTTT
(


cat *hapinfo.txt > Merge.hapInfo.txt
perl ../hapinfoMerge.pl Merge.hapInfo.txt

perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.3 > MHB.Mouse.R2-0.3.bed &
perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.5 > MHB.Mouse.R2-0.5.bed &
perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.7 > MHB.Mouse.R2-0.7.bed &
perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.9 > MHB.Mouse.R2-0.9.bed &


qsub SRR2011294_1.fastq.gz.job
for i in 4 5 6 7
do
qsub SRR201129$i\_1.fastq.gz.job
done


perl hapinfo2BlocAvgR2.pl mESC.sort.hapinfo.txt.SumUniq > mESC.sort.hapinfo.R2.txt
perl hapinfo2BlocAvgR2.pl Adult.sort.hapinfo.txt.SumUniq > Adult.sort.hapinfo.R2.txt

grep chr10:100004267-100004288 *sort.hapinfo.R2.txt
grep chr10:100004267-100004288 *hapinfo.txt.SumUniq
grep chr10:100004267-100004288 output.mf

grep chr10:68979553-68979781 *sort.hapinfo.R2.txt
grep chr10:68979553-68979781 *hapinfo.txt.SumUniq
grep chr10:68979553-68979781 output.mf

# R2 for each block
/home/shg047/oasis/mouse/mergeHapinfo/mESC.sort.hapinfo.R2.txt
/home/shg047/oasis/mouse/mergeHapinfo/Adult.sort.hapinfo.R2.txt
# 5me level for each block
/home/shg047/oasis/mouse/mergeHapinfo/mf/output.mf


修改hapinfo2mhb.pl的程序，输出R2,D,D-primer.
写一个程序，进行haploinfo.merge

sort Adult.hapinfo.txt > Adult.sort.hapinfo.txt
sort Adult.hapinfo.txt > Adult.sort.hapinfo.txt

data<-read.table("output.mf")
subset(data,Adult>0.6 & mESC<0.3)
 
qsub hapin2R2-mESC.job
qsub hapin2R2-adult.job


#!/bin/csh
#PBS -n hapinfo2r2
#PBS -q pdafm
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o Adult.R2.log
#PBS -e Adult.R2.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /home/shg047/oasis/mouse/mergeHapinfo
perl ~/bin/hapinfo2LDR2ByBed.pl /home/shg047/oasis/mouse/RD/hapinfo/Mouse.MHB.Alice.RD90_80up_R0.5.bed Adult.hapinfo.txt > Adult.R2.txt
cd /home/shg047/oasis/mouse/mergeHapinfo
perl ~/bin/hapinfo2LDR2ByBed.pl /home/shg047/oasis/mouse/RD/hapinfo/Mouse.MHB.Alice.RD90_80up_R0.5.bed mESC.hapinfo.txt > mESC.R2.txt
cd /home/shg047/oasis/mouse/mergeHapinfo/mf
perl /home/shg047/bin/hapinfo2mf.pl /home/shg047/oasis/mouse/mergeHapinfo/mf > output.mf

sort Adult.hapinfo.txt > Adult.sort.hapinfo.txt
sort mESC.hapinfo.txt > mESC.sort.hapinfo.txt & 
perl ~/bin/hapinfo2LDR2.pl rlt chr1:10006430-10006441 < Adult.sort.hapinfo.txt


reformat.sh in=test.fq out=test.phred33.fq qin=64 qout=33

perl ~/bin/hapinfo2LDR2.pl rlt chr10:107810947-107810953  < test.sort

chr1:103067130-103067142

perl ~/bin/hapinfo2LDR2.pl rlt chr1:50329971-50335398 < Mouse.MHB.Alice.chr1.bam.hapInfo.txt
? get line-req
chr1:50329971-50335398

perl /home/shg047/bin/bam2hapInfo2PBS.pl bam2hapinfo.config  submit bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt



/home/shg047/oasis/mouse/hapinfo/Mouse.MHB.Alice.RD90_80up_R0.5.bed

/home/shg047/oasis/mouse/RD/hapinfo/



mergedBam2hapInfoV2.pl
ls SRR2011294

qsub SRR2011294

qsub SRR2011294_1.fastq.gz.job
qsub SRR2011295_1.fastq.gz.job
qsub SRR2011296_1.fastq.gz.job
qsub SRR2011297_1.fastq.gz.job


for i in `ls *cov`
do
awk '{sum +=$4;n++}END{print sum/n}' $i
done

ls -l  | awk -F : '{sum+=$5} END {print "AVG=",sum/NR}'


for i in {1..19} X Y M
do
perl ~/bin/hapinfo2mhb.pl Mouse.MHB.Alice.chr$i.bam.hapInfo.txt 0.3 >> Mouse.MHB.Alice.RD90_80up_R0.3.bed
done

hapinfo2mld.pl


perl hapinfo2mld.pl Mouse.MHB.Alice.chr$i.bam.hapInfo.txt >> Mouse.MHB.r0.25.bed

perl ~/bin/hapinfo2LDR2.pl rlt  chr14:100392580-100392941 < Mouse.MHB.Alice.chr14.bam.hapInfo.txt
samtools tview -p chr14:32546776 Mouse.MHB.Alice.chr14.bam /home/shg047/oasis/db/mm9/mm9.fa

chr14:32546776-32546898


perl /home/shg047/bin/bam2hapInfo2PBS.pl bam2hapinfo.config  submit bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt


for i in {1..19} X Y M
do
samtools index Mouse.MHB.Alice.chr$i.bam &
done

cd /home/shg047/oasis/db/mm9/mm9

MM9.CpG.positions.txt

/home/shg047/oasis/db/mm9/MM9.CpG.positions.txt


cat *cg.pos > MM9.CpG.positions.txt




/home/shg047/oasis/mouse/RD/MouseRD90UP80.bed

source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Mmusculus.UCSC.mm9")

chrs <- names(BSgenome.Mmusculus.UCSC.mm9)[!grepl("random", names(BSgenome.Mmusculus.UCSC.mm9))] #filter out the "upstream" chromosomes
CGs <- lapply(chrs, function(x) start(matchPattern("CG", BSgenome.Mmusculus.UCSC.mm9[[x]])))
names(CGs) <- chrs
seqlengths(genome)
genome$chr1 # same as genome[["chr1"]]


for i in {1..19} X Y M
do
awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD30_80up.bed
# awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD50_80up.bed
# awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD70_80up.bed
# awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD90_80up.bed
done

for i in {1..19} X Y M
do
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' Mouse.MHB.Alice.chr$i.bam.gencov.bed 
done


awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD30_80up.bed
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD50_80up.bed
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD70_80up.bed
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD90_80up.bed




for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$4>29 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD30.bed &
# awk '$4>49 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD50.bed &
# awk '$4>69 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD70.bed &
# awk '$4>89 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD90.bed &
done

for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD30.bed > Mouse.MHB.Alice.chr$i.bam.RD30_80up.bed &
#awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD50.bed > Mouse.MHB.Alice.chr$i.bam.RD50_80up.bed &
#awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD70.bed > Mouse.MHB.Alice.chr$i.bam.RD70_80up.bed &
#awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD90.bed > Mouse.MHB.Alice.chr$i.bam.RD90_80up.bed &
done



for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD50.bed > Mouse.MHB.Alice.chr$i.bam.RD50_80up.bed &
done


#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=1
#PBS -o output.file
#PBS -me
#PBS -S /bin/bash

set -e

module add apps/vcftools-0.1.12b
export PERL5LIB=/cm/shared/apps/VCFTOOLS-0.1.12b/vcftools_0.1.12b/perl/
module add apps/tabix-0.2.6

wd="/projects/UK10K/GOYA/imputed"
snplist="/panfs/panasas01/sscm/eprcr/GOYA/speliotes_snps.txt"
tempdir="${HOME}/tempdir"
out="${HOME}/output.vcf.gz"

#  Run the code

cd ${wd}
mkdir -p ${tempdir}


find . -maxdepth 1 -name '*.vcf.gz' | while read f; do

	# Output filenames - strip path and extensions
	fn=$(basename "${f}")
	fn="${fn%.*.*}"

	# Create copy of file with tab delimiter instead of space delimiter
	# Note - you may want to just do this to all the files because vcftools won't recognise them as they currently are
	zcat ${f} | tr ' ' '\t' | sed 's/exm-//g' | gzip -c > ${tempdir}/temp.vcf.gz

	# Extract SNPs
	vcftools --gzvcf ${tempdir}/temp.vcf.gz --snps ${snplist} --recode --keep-INFO-all --out ${tempdir}/${fn}

	# Prepare for merging
	bgzip ${tempdir}/${fn}.recode.vcf
	tabix ${tempdir}/${fn}.recode.vcf.gz
done
rm ${tempdir}/temp.vcf.gz

# Merge all output files
# note: maybe better to just select files that have an rs in them:
# zgrep rs tempdir/*.vcf.gz
vcf-merge ${tempdir}/*vcf.gz | bgzip -c > ${out}


samtools view -b -q 20 Mouse.MHB.Alice.MergeBam.sort.bam chr1 > Mouse.MHB.Alice.chr1.bam &
samtools view -b -q 20 Mouse.MHB.Alice.MergeBam.sort.bam chr2 > Mouse.MHB.Alice.chr1.bam &
awk '$4>29 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr1.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr1.bam.RD10.bed &
awk '$4>29 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr2.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr2.bam.RD10.bed &


./configure --prefix=/home/shg047/software/R-3.3.1 '--with-cairo' \
 '--with-jpeglib' '--with-readline' '--with-tcltk' '--with-x=no'\
 '--with-blas' '--with-lapack' '--enable-R-profiling' '--with-tiff=yes'\
 '--enable-R-shlib'\
 '--enable-memory-profiling'
 
SRR1248497_1.fas
SRR630206
SRR833525_1


samtools index Mouse.MHB.Alice.MergeBam.sort.bam
for i in {1..22} X Y M
do
samtools view -b -q 20 Mouse.MHB.Alice.MergeBam.sort.bam chr$i > Mouse.MHB.Alice.chr$i.bam &
done

for i in {01,02,04,05,06,07,09,10,11,12}
do
echo $i
done

for i in {01,02,04,05,06,07,09,10,11,12}
do
cat s_*_1_Indx$i.txt.gz > Indx$i.read1.fq.gz &
cat s_*_2_Indx$i.txt.gz > Indx$i.read2.fq.gz &
done

 
zcat *1_Indx02*.gz > Indx02.s.1.fq.gz 
zcat *2_Indx02*.gz > Indx02.s.2.fq.gz 
Indx10

perl ../../fastq/bismarkPBS.pl SamConfig.txt submit

mv	s_1_Indx01.txt.gz	Indx01_s_1.txt.gz
mv	s_1_Indx02.txt.gz	Indx02_s_1.txt.gz
mv	s_1_Indx04.txt.gz	Indx04_s_1.txt.gz
mv	s_1_Indx05.txt.gz	Indx05_s_1.txt.gz
mv	s_1_Indx06.txt.gz	Indx06_s_1.txt.gz
mv	s_1_Indx07.txt.gz	Indx07_s_1.txt.gz
mv	s_1_Indx09.txt.gz	Indx09_s_1.txt.gz
mv	s_1_Indx10.txt.gz	Indx10_s_1.txt.gz
mv	s_1_Indx11.txt.gz	Indx11_s_1.txt.gz
mv	s_1_Indx12.txt.gz	Indx12_s_1.txt.gz
mv	s_2_Indx01.txt.gz	Indx01_s_2.txt.gz
mv	s_2_Indx02.txt.gz	Indx02_s_2.txt.gz
mv	s_2_Indx04.txt.gz	Indx04_s_2.txt.gz
mv	s_2_Indx05.txt.gz	Indx05_s_2.txt.gz
mv	s_2_Indx06.txt.gz	Indx06_s_2.txt.gz
mv	s_2_Indx07.txt.gz	Indx07_s_2.txt.gz
mv	s_2_Indx09.txt.gz	Indx09_s_2.txt.gz
mv	s_2_Indx10.txt.gz	Indx10_s_2.txt.gz
mv	s_2_Indx11.txt.gz	Indx11_s_2.txt.gz
mv	s_2_Indx12.txt.gz	Indx12_s_2.txt.gz


for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$4>29 { print $1"\t"$2"\t"$3}' merge.chr$i.bam.gencov.bed | bedtools merge -d 10 -i - > merge.chr$i.bamRD10.bed
awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' merge.chr$i.bamRD10.bed > merge.chr$i.bam.RD10_80up.bed
done


/home/shg047/oasis/mouse/RD
Mouse.MHB.Alice.MergeBam.bam

/oasis/tscc/scratch/zhl002/WGBS_mouse

Previous Script:
1, Download or Prepare SRA Download Configure  (http://www.ebi.ac.uk/ena/data/view/SRP028600)
2, perl fastqDownload.pl SamConfig.txt 8 submit   (fastq-dump --split-files --gzip)
3, trim_galore --phred33 --fastqc --stringency 3 -q 20 --trim1 --length 20 --gzip --clip_R1 2 --three_prime_clip_R1 2 --illumina SRR299055_1.fastq.gz --output_dir ../fastq_trim
4, bismark_genome_preparation ./     # merge all the fa to one file (mm9.fa or hg19.fa)
5, 
--------------------------------------------------------------------------------------------------------------






gawk '{ sum += $1 }; END { print sum }' file
gawk -F: '{ print $1 }' /etc/passwd

ls *_report.txt | awk -FS="_" 'BEGIN { FS = "_" } {print $1}' | sort -u | wc -l
ls *_report.txt | awk -F_ '{print $1}' | sort -u | wc -l

perl -p -i -e 's/walltime=167/walltime=48/g' *job
perl -p -i -e 's/ppn=16/ppn=4/g' *job
perl -p -i -e 's/multicore 6/multicore 2/g' *job

perl bismarkPBS.pl FastqMatchConfig.txt submit
bismark --bowtie2 --phred33-quals --fastq -L 30 -N 1 --multicore 6 /home/shg047/db/mm9 -1 ../fastq_trim/SRR833525_1.fastq.gz_val_1.fq.gz -2 ../fastq_trim/SRR833525_2.fastq.gz_val_2.fq.gz  -
qsub SRR833525_1.fastq.gz.job
qsub SRR833526_1.fastq.gz.job

/home/shg047/oasis/db/mm9
/home/shg047/oasis/mouse/aliceSRR2876131_1_val_1.fq.gz

for i in 01 02 04 05 06 07 09 10 11 12 
do 
cat s_*_1_Indx$i.txt.gz > s_1_Indx$i.txt.gz &
cat s_*_2_Indx$i.txt.gz > s_2_Indx$i.txt.gz &
done 

cat s_*_4_Indx01.txt.gz > s_4_Indx01.txt.gz &
cat s_*_5_Indx01.txt.gz > s_5_Indx01.txt.gz &
cat s_*_6_Indx01.txt.gz > s_6_Indx01.txt.gz &
cat s_*_7_Indx01.txt.gz > s_7_Indx01.txt.gz &
cat s_*_9_Indx01.txt.gz > s_9_Indx01.txt.gz &
cat s_*_10_Indx01.txt.gz > s_10_Indx01.txt.gz &
cat s_*_11_Indx01.txt.gz > s_11_Indx01.txt.gz &
cat s_*_12_Indx01.txt.gz > s_12_Indx01.txt.gz &

gzcat file1.fastq.gz file2.fastq.gz | gzip > merged.fastq.gz

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz

/media/LTS_60T/SeqStore2016/130104_SN1001
/media/LTS_60T/SeqStore2016/130104_SN1001


1, SRR2136776 require trim and alignment

/home/shg047/oasis/db/mm9

wget http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.16.3.tar.gz
cd /home/shg047/oasis/db/mm9

# 2016-08-22
perl -p -i -e 's/walltime=7/walltime=8/g' *job
cd /home/shg047/oasis/mouse
perl fastqDownload.pl SamConfig.txt 8 submit



# 2016-08-17

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("tutorial.Rmd")'
R -e 'library("markdown");rpubsUpload("tutorial","tutorial.html")'

library("knitr")
library("rmarkdown")
pandoc("tutorial.md",format="MediaWiki")


rmarkdown::render('tutorial.Rmd')  # output: md_document


# 2016-08-16
R -e 'library("rmarkdown");library("knitr");rmarkdown::render("NormalDevconJuly.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","NormalDevconJuly.html")'
library("knitr")
library("rmarkdown")
rmarkdown::render('DeconvolutionMixture.Rmd')
pandoc("DeconvolutionMixture.md",format="MediaWiki")

1. Download geneSCF: 
2, tar zxvf 
./geneSCF -m=update -i=test/H0.list -o=test/output/ -t=sym -db=GO_MF -bg=20000 --plot=yes -org=goa_human

ls *job |cut -b 1-11 | sort -u | wc -l
# 2016-08-17
perl bam2amfBybedByPileOMethPBS.pl mhb.bed /home/shg047/oasis/N37/sortBam
perl bam2amfBybedByPileOMethPBS.pl mhb.bed /home/shg047/oasis/SALK/bam
#2016-08-12
cd /home/shg047/oasis/N37/sortBam
cd /home/shg047/oasis/SALK/bam

PileOMeth extract -q 10 -p 5 --minDepth 5 -l mhb.bed ~/oasis/db/hg19_lambda.fa /home/shg047/oasis/N37/sortBam/Indx01.sort.bam


cd /home/shg047/oasis/monod/Figure3
PileOMeth extract -q 10 -p 5 --minDepth 5 -l mhb.bed ~/oasis/db/hg19_lambda.fa
PileOMeth extract  -q 10 -p 5 --minDepth 5 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o 7-P-3.bedGraph 
#2016-08-9
data<-read.table("xx",head=T)
y<-paste(data[,5],data[,6],sep="-")
table(y)
bedtools sort -i xx.bed > xx.sort.bed
bedtools closest -fu -D "ref" -a xx.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9}' | sort -u > xx.anno.bed
# 
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Figure-gs.pdf Figure*.pdf
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Supp.Figure-gs.pdf Supp*Fig*.pdf

pdftk Figure*.pdf cat output Figure.pdf
pdftk Supp*Fig*.pdf cat output Supp.Figure.pdf

# 2016-07-29

bedtools sort -i HumanTissueGSI.bed > HumanTissueGSI.sort.bed
bedtools closest -fd -D ref -a HumanTissueGSI.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed


# 2016-07-21
bedtools closest -a 
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=Colon.Cancer.Hg19  
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=ESCA-MH450
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=mm9
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=MONOD

# 2016-07-20
sudo /etc/init.d/apache2 stop
sudo /etc/init.d/mysql stop
sudo /etc/init.d/proftpd stop
sudo /opt/lampp/lampp start
sudo service vsftpd stop
sudo /opt/lampp/lampp stop


cd ~/Downloads
sudo chmod 777 -R bitnami-wordpress-3.9.1-1-module-linux-x64-installer.run
./bitnami-wordpress-3.9.1-1-module-linux-x64-installer.run



install.packages('gsheet')
library("gsheet")
gsheet2tbl('https://drive.google.com/open?id=0B2TJ0NCGGrdpVlBXcWJOY2hRdFk')


screen -ls | grep pts | cut -d. -f1 | awk '{print $1}' | xargs kill

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("Deconv.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","Deconv.html")'
R -e 'library("markdown");rpubsUpload("normalDev","Deconv.html")'
R -e 'library("rmarkdown");library("knitr");rmarkdown::render("RandomForst.Rmd")'


 #!/bin/csh
 #PBS -N bedgraph2matrix
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/hapinfo/June
 R -e 'library("rmarkdown");library("knitr");rmarkdown::render("RandomForst.Rmd")'



R -e 'library("rmarkdown");library("knitr");rmarkdown::render("RandomForest.Rmd")'



---
title: "Habits"
output:
  md_document:
    variant: markdown_github
---
---
title: "Planets"
author: "Manoj Kumar"
date: "March 3, 2016"
output: v
  html_document:
    toc: true # table of content true
    depth: 3  # upto three depths of headings (specified by #, ## and ###)
    number_sections: true  ## if you want number sections at each table header
    theme: united  # many options for theme, this one is my favorite.
    highlight: tango  # specifies the syntax highlighting style
---

2016-07-11
exampl.Rmd

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("exampl.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","ColonDevconJuly.Rmd.html")'

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("ColonDevconJuly.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","ColonDevconJuly.Rmd.html")'

cmake .. -DRSTUDIO_TARGET=Server -DCMAKE_BUILD_TYPE=Release  CMAKE_INSTALL_PREFIX=./
W: GPG error: http://ppa.launchpad.net trusty InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY 4F191A5A8844C542
sudo mv /var/lib/dpkg/lock
sudo mv /var/cache/apt/archives/lock 
sudo mv /var/cache/apt/archives/lock_bak
sudo rm /var/lib/apt/lists/* -vf
sudo apt-get update && sudo apt-get upgrade
sudo add-apt-repository --remove ppa:kernel-ppa/ppa

/etc/apt/sources.list
deb http://archive.ubuntu.com/ubuntu trusty main
lsb_release -sc
sudo apt-get update
sudo apt-get install gedbi-core

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("NormalDevconJuly.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","NormalDevconJuly.html")'

R CMD BATCH RunchampDMR.R &
track type=bigWig name=proteinA smoothingWindow=4 color=123,100,50 autoScale=on viewLimits=1:200 visibility=full windowingFunction=maximum bigDataUrl=https://projects/files/file.bw
48502
chr1:23730471-23730518
#!/usr/bin/sh
for i in `ls *.bedGraph`
do 
(awk 'NR!=1' $i| sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4}') > $i.dehsort 
bedGraphToBigWig $i.dehsort ~/oasis/db/hg19.chrom.sizes $i.dehsort.bw 
bigWigAverageOverBed $i.dehsort.bw  /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed  $i.dehsort.bw.mf
rm *.dehsort
done
head -n 20 /home/shg047/oasis/monod/BAM/MF_PileOMeth/CRC-P-002_CpG.bedGraph
head -n 20 /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/6-P-2.sorted.clipped_CpG.bedGraph
 samtools tview /home/shg047/oasis/monod/BAM/rename/CRC-P-024.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:23730471-23730518
cp /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/mf2matrix.pl  ./
cp /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/bedgraph2matrix.pl ./
cp /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/bw2mf.pl ./
 
 
 
 samtools tview RRBS-7P23.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr10:100027957
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p chr1:17092527-17092562
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p 
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:10496-10498

samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:935977-936090
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066


 
 #!/bin/csh
 #PBS -N bedgraph2matrix
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

perl bedgraph2matrix.pl > bedgraphMatrix.txt
 
 

 
 
 

 

 #!/bin/csh
 #PBS -N deconv
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/hapinfo/June
 #Rscript --vanilla colon.updateGSI.R
 #Rscript --vanilla lung.updateGSI.R
 Rscript --vanilla normal.updateGSI.R

 
 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

(tail -n +2 6-T-4.sorted.clipped_CpG.bedGraph| sort -k 1,1 -k2,2n) > 6-T-4.sorted.clipped_CpG.bedGraph


(awk 'NR!=1' 6-T-4.sorted.clipped_CpG.bedGraph| sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4}') > 6-T-4.sorted.clipped_CpG.sort.bedGraph
bedGraphToBigWig 6-T-4.sorted.clipped_CpG.sort.bedGraph ~/oasis/db/hg19.chrom.sizes 6-T-4.sorted.clipped_CpG.sort.bedGraph.bw

 
 

# copy simulated CCTmix and LCTmix to hapinfo fold, run hapinfo2mhl perl script
cp
cp 

cd /home/shg047/oasis/monod/hapinfo/hapinfo
ls *hapInfo.txt > Hapinfo_File_list

qsub hapinfo2mhl.pbs




perl CCTmixturePBS.pl 6-P-Merge.hapinfo.txt WB.hapInfo.txt
perl LCTmixturePBS.pl 7-P-Merge.hapinfo.txt WB.hapInfo.txt
  
echo "perl /home/shg047/bin/hapinfo2mhl.pl Hapinfo_File_list > MHL.output.txt" |qsub -N hapinfo2MHL -q pdafm -l walltime=8:00:00
echo "perl /home/shg047/bin/hapinfo2mhl.pl Hapinfo_File_list > MHL.output.txt" |qsub -N h2m

hapInfo.txt

 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group 
 perl /home/shg047/bin/hapinfo2mhl.pl Hapinfo_File_list > MHL.mixture.simulation.txt


#!/usr/bin/env perl
use strict;
my @f=qw /0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/;
foreach my $f(@f){
my $cmd1="perl ../randomSampleFromHaploInfo.pl int(24850020*f) WB.hapinfo.merge > WB";
my $cmd2="perl ../randomSampleFromHaploInfo.pl int(1386115*(1-f)) 6-T-Merge.hapinfo.txt > 6T";
my $cmd3="cat WB 6T > Wcmix.f.hapinfo";
print "$cmd1\n$cmd2\n$cmd3\n";
}

/home/shg047/oasis/monod/deconvolution/simulation
# 2016-06-21

 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/colon.updateGSI.R
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lung.updateGSI.R
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lungcancer.updateGSI.withoutLCT.R
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/coloncancer.updateGSI.withoutCCT.R

echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/colon.updateGSI.R" | qsub -q pdafm
echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lung.updateGSI.R" | qsub -q pdafm
echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lungcancer.updateGSI.withoutLCT.R" | qsub -q pdafm
echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/coloncancer.updateGSI.withoutCCT.R" | qsub -q pdafm
perl -p -i -e 's/\>0.01/\>0.05/' colon.updateGSI.R
perl -p -i -e 's/\>0.01/\>0.05/' lung.updateGSI.R
perl -p -i -e 's/\>0.01/\>0.05/' lungcancer.updateGSI.withoutLCT.R
perl -p -i -e 's/\>0.01/\>0.05/' coloncancer.updateGSI.withoutCCT.R
perl format.pl lsMHL.cor.txt | bedtools sort > tsMHL.cor.sort.bed
perl format.pl lsMHL.cor.txt | bedtools sort >lsMHL.cor.sort.bed
bedtools closest -fu -D "ref" -a lsMHL.cor.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$11}' | sort -u > lsMHL.cor.sort.anno.bed
perl format.pl php.bed |bedtools sort > php.sort.bed
bedtools closest -fu -D "ref" -a php.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$9}' | sort -u > php.bed.anno.bed
bedtools intersect -a php.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$9}' | sort -u > php.bed.anno.bed
perl format.pl | bedtools sort > tsMHL.cor.sort.bed


bedtools closest -fu -D "ref" -a WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9}' | sort -u > WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.gene.anno.bed
bedtools closest -fu -D "ref" -a tsMHL.cor.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9}' | sort -u > tsMHL.cor.sort.anno.bed

# 2016-06-13

system("wget https://cran.r-project.org/src/contrib/pbkrtest_0.4-6.tar.gz")
system("wget https://cran.r-project.org/src/contrib/quantreg_5.26.tar.gz")
system("wget https://cran.r-project.org/src/contrib/lme4_1.1-12.tar.gz")
system("wget https://cran.r-project.org/src/contrib/minqa_1.2.4.tar.gz")
system("wget https://cran.r-project.org/src/contrib/nloptr_1.0.4.tar.gz")
system("wget https://cran.r-project.org/src/contrib/RcppEigen_0.3.2.8.1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/SparseM_1.7.tar.gz")
system("wget https://cran.r-project.org/src/contrib/MatrixModels_0.4-1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/ggplot2movies_0.0.1.tar.gz")
install.packages("minqa_1.2.4.tar.gz")
install.packages("nloptr_1.0.4.tar.gz")
install.packages("RcppEigen_0.3.2.8.1.tar.gz")
install.packages("lme4_1.1-12.tar.gz")
install.packages("pbkrtest_0.4-6.tar.gz")
install.packages("SparseM_1.7.tar.gz")
install.packages("MatrixModels_0.4-1.tar.gz")
install.packages("quantreg_5.26.tar.gz")
install.packages("car_2.1-2.tar.gz")
install.packages("ggplot2movies_0.0.1.tar.gz")

mix="simulation.mixture.inp.txt"
pure="simulation.signature.inp.txt"
class="simulation.lab.inp.txt"
Rscript fix_matrices.R $pure $mix $class
#R CMD ~/R/x86_64-pc-linux-gnu-library/3.2/Rserve/libs/Rserve --no-save
R CMD /home/shicheng/R/x86_64-pc-linux-gnu-library/3.2/Rserve/libs --no-save
java -Xmx3g -Xms3g -jar CIBERSORT.jar -M $mix.mod -P $pure.mod -c $class > outputs


# 2016-06-13
cd /home/shg047/deconvolution/cibersort
java -Xmx3g -Xms3g -jar /home/shg047/software/CIBERSORT/CIBERSORT.jar -M 

R CMD /home/shg047/software/R-3.3.0/library/Rserve/R/Rserve --no-save
R CMD /home/shg047/software/R-3.3.0/library/Rserve/libs/Rserve --no-save
R CMD ~/R/x86_64-pc-linux-gnu-library/3.2/Rserve/libs/Rserve --no-save

java -Xmx3g -Xms3g -jar /home/shg047/software/CIBERSORT/CIBERSORT.jar -M 

cp signatures.lung.data.txt.mod ~
cp test.lung.data.txt.mod ~


java -Xmx3g -Xms3g -jar /home/ddiep/softwares/CIBERSORT/CIBERSORT.jar -M colon.mixture.inp.txt -B colon.signature.inp.txt.mod

# 2016-06-13
the difference between nucleo and whole-cell 

perturbation screen 
genotype to gene expression profile
how to select postive cells (50% percent)
permuation to remove false postive gene set and check the remain set in the comparison





## Use IGV in Windows connecting with Linux server. 

In windows, you don't have terminal, you can use putty to connect Linux (installing IGV) server.
Be careful, You can not use Winscp to use Putty, since you can not set X11 forwarding for putty in Winscp. 

# before connect to linux sever with putty, please install and open Xming (Windows)

in the windows: install winscp and Xming
X11 forwarding in putty: https://wiki.utdallas.edu/wiki/display/FAQ/X11+Forwarding+using+Xming+and+PuTTY
in the linxu server install IGV for linux
Download (Binary Distribution): http://www.broadinstitute.org/igv/download
open Xming in windowns and Go to server IGV fold and run: java -Xmx750m -jar igv.jar


## Use IGV in Linux or Mac connecting with Linux server. 

ssh -X shg047@genome-miner.ucsd.edu


# 2016-06-01

/home/shg047/oasis/WGBSDB/PBMC/WGBS.select.1501.txt



http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer8/tracks_hg19/

# PBMC and Blood cells
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_BCell/tracks_hg19/Human_BCell.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_CD133HSC/tracks_hg19/Human_CD133HSC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_HSPC/tracks_hg19/Human_HSPC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_Neut/tracks_hg19/Human_Neut.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-100yr/tracks_hg19/Human_CD4T-100yr.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-Newborn/tracks_hg19/Human_CD4T-Newborn.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Li-PBMC-2010/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-DNMT3BMut-2012/Human_BCell-Healthy/tracks_hg19/Human_BCell-Healthy.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_Macrophage/tracks_hg19/Human_Macrophage.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_Tcell/tracks_hg19/Human_Tcell.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_NK/tracks_hg19/Human_NK.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_HSC/tracks_hg19/Human_HSC.meth.bw & 
# Lung
wget http://smithlab.usc.edu/methbase/data/Xie-Human-2013/Human_IMR90/tracks_hg19/Human_IMR90.meth.bw
wget https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Lung/UCSD.Lung.Bisulfite-Seq.STL002.wig.gz
wget https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/IMR90_Cell_Line/UCSD.IMR90.Bisulfite-Seq.combined.wig.gz

# Colon
wget http://smithlab.usc.edu/methbase/data/Berman-Human-2012/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Ziller-Human-2013/Human_Colon_Tumor_Primary/tracks_hg19/Human_Colon_Tumor_Primary.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hansen-Human-2011/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw &

# Total 45 samples 
https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Adipose_Tissue/UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.wig.gz

# bedgraph to bigwig
cd /home/shg047/oasis/Estellar2016/MF
for i in `ls *bedGraph`
do
(head -n 1 $i && tail -n +2 $i | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $i.sort
bedGraphToBigWig $i.sort ~/oasis/db/hg19/hg19.chrom.sizes $i.bw 
done


# bigWigToBedGraph
for i in `ls *bw`
do
bigWigToBedGraph $i $i.bg &
done

# bigWigAverageOverBed
for i in `ls *bw`
do
bigWigAverageOverBed $i /home/shg047/oasis/WGBSDB/PBMC/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.renew.bed $i.aob &
done

chr18:74692045–74692213
chr10:79936919–79937042




# 2016-05-29

[shg047@tscc-login1 MF_PileOMeth]$ wc -l CTR106.bedGraphV4
20653494 CTR106.bedGraphV4
[shg047@tscc-login1 MF_PileOMeth]$ wc -l CTR107.bedGraphV4
22745883 CTR107.bedGraphV4



cd /home/shg047/oasis/DennisLo2015/MF_PileOMeth

cp /home/shg047/oasis/monod/wig/bedGraph2V4.pl /home/shg047/oasis/DennisLo2015/MF_PileOMeth

CTR101_trimmed.fq.gz_bismark_bt2.sort_CpG.bedGraph

track type="bedGraph" description="/home/shg047/oasis/DennisLo2015/MF_PileOMeth/CTR108_trimmed.fq.gz_bismark_bt2.sort CpG merged methylation levels"
chr10   61056   61058   0       0       1
chr10   61332   61334   100     1       0
chr10   61652   61654   100     2       0
chr10   63584   63586   75      3       1
chr10   64773   64775   100     9       0

vim bedGraph2V4.pl





# check what happened to the MF and MHL result

chr10:100027918-100027944       TTTTTTT 2       100027918,100027922,100027925,100027927,100027930,100027938,100027944
chr10:100027957-100027992       C       1       100027957
chr10:100027957-100027992       T       1       100027957
chr10:100121047-100121066       TTT     1       100121047,100121049,100121066
chr10:100121047-100121066       CCC     3       100121047,100121049,100121066
chr10:100121047-100121066       TT      1       100121047,100121049
chr10:100151033-100151068       T       2       100151068
chr10:100206398-100206510       TTTTTT  1       100206398,100206407,100206455,100206475,100206479,100206487
chr10:100227616-100227677       TTTT    3       100227616,100227620,100227625,100227651
chr10:100227616-100227677       TTTTTT  1       100227616,100227620,100227625,100227651,100227666,100227677
chr10:100227698-100227747       TTTTTTTT        1       100227698,100227707,100227710,100227714,100227719,100227724,100227727,100227740
chr10:100227781-100227823       TT      1       100227814,100227823
chr10:100227781-100227823       TTT     1       100227781,100227787,100227814
chr10:100992364-100992404       CT      1       100992364,100992404
chr10:100992364-100992404       C       3       100992364
chr10:100995252-100995487       T       1       100995487
chr10:100995942-100996028       T       3       100996028
chr10:101089155-101089163       TTTT    1       101089155,101089157,101089160,101089163
chr10:101089204-101089305       TTTT    2       101089204,101089209,101089243,101089250
chr10:101089204-101089305       TTTTTTTCTTT     1       101089204,101089209,101089243,101089250,101089254,101089257,101089260,101089279,101089284,101089303,101089305
chr10:101089204-101089305       TTTTCTTT        2       101089250,101089254,101089257,101089260,101089279,101089284,101089303,101089305
chr10:101089382-101089424       T       3       101089382
chr10:101089443-101089461       TTTT    2       101089443,101089454,101089456,101089461
chr10:101089481-101089508       T       2       101089481
chr10:101089481-101089508       TTTTTT  1       101089481,101089485,101089491,101089494,101089506,101089508
chr10:101089526-101089565       TTTT    1       101089532,101089537,101089554,101089565
chr10:101089526-101089565       T       1       101089526
chr10:101089526-101089565       TTTTT   2       101089526,101089532,101089537,101089554,101089565
chr10:101089569-101089596       TTT     3       101089569,101089591,101089596
chr10:101089668-101089771       TTTTTT  1       101089746,101089759,101089761,101089765,101089769,101089771
chr10:101089816-101089850       TTTTT   2       101089816,101089820,101089824,101089829,101089850
chr10:101089870-101089906       TTTT    1       101089870,101089886,101089895,101089906
chr10:101089870-101089906       T       1       101089906
chr10:101089870-101089906       TTT     1       101089870,101089886,101089895
chr10:101089928-101089935       TTT     2       101089928,101089932,101089935
chr10:101089949-101089965       TTTTT   2       101089949,101089951,101089954,101089956,101089965
chr10:101190105-101190169       TTTT    1       101190105,101190120,101190141,101190169
chr10:101190105-101190169       TTTTT   2       101190105,101190120,101190141,101190157,101190169
chr10:101190297-101190505       TTT     1       101190297,101190310,101190324
chr10:101190297-101190505       TTTTTTTTTTTTTTTT        1       101190297,101190310,101190324,101190420,101190423,101190429,101190435,101190437,101190451,101190456,101190458,101190475,101190477,101190489,101190502,101190505
chr10:101190297-101190505       TTTTTTTTTTTTT   2       101190420,101190423,101190429,101190435,101190437,101190451,101190456,101190458,101190475,101190477,101190489,101190502,101190505
chr10:101190297-101190505       TTT     1       101190384,101190389,101190420
chr10:101190297-101190505       TTTT    1       101190360,101190384,101190389,101190420
chr10:101190297-101190505       TTTTTTTTTTT     1       101190420,101190423,101190429,101190435,101190437,101190451,101190456,101190458,101190475,101190477,101190489
chr10:101190527-101190864       TTTTTTTTTTTTTTT 1       101190595,101190604,101190610,101190619,101190635,101190646,101190674,101190691,101190695,101190699,101190702,101190727,101190734,101190751,101190763
chr10:101190527-101190864       TTTTTTTTTT      1       101190527,101190558,101190565,101190595,101190604,101190610,101190619,101190635,101190646,101190674
chr10:101190527-101190864       TTTTT   1       101190776,101190798,101190801,101190815,101190827
chr10:101190527-101190864       TT      1       101190841,101190864
chr10:101190527-101190864       TTT     4       101190527,101190558,101190565


RRBS-7P23.hapInfo.txt
ls RRBS-7P23*
less RRBS-7P23.sorted.clipped_CpG.bedGraph

grep 540439 RRBS-7P23.hapInfo.txt | grep chr1
grep 101190298 RRBS-7P23.sorted.clipped_CpG.bedGraph
grep 101190120 RRBS-7P23.sorted.clipped_CpG.bedGraph
grep 101190120 RRBS-7P23.sorted.clipped_CpG.bedGraph



less -S RRBS-7P23.hapInfo.txt

less 

samtools tview RRBS-7P23.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr10:100027957
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p chr1:17092527-17092562
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p 
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:935977-936090
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066





samtools tview RRBS-7P23.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:850885


chr1    850885

chr10:100027957-100027992

/home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/6-P-9.sorted.clipped_CpG.bedGraph
wc -l 6-P-9_CpG.bedGraph

cd oasis/monod/bam/RRBS1/bam/

chr1:10524-10526
samtools tview 6-P-9.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:10524
samtools tview 6-P-9.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:724245

chr1:724245


wc -l ../../../wig/6-P-9_CpG.bedGraph
1178863

head ../../../bedGraph/6-P-9_CpG.bedGraph

head 6-P-9_CpG.bedGraph
1394266 6-P-9.sorted.clipped_CpG.bedGraph





http://smithlab.usc.edu/methbase/data/Hon-Human-2012/Human_HCC1954/tracks_hg19/Human_HCC1954.meth.bw
# 2016-05-27
Install BioPerl Without Root Privileges in Ubuntu/Linxu
https://www.biostars.org/p/193668/

< How To Install BioPerl Without Root Privileges in Ubuntu/Linxu>

# Install certain nessary library
perl -MCPAN -Mlocal::lib -e 'CPAN::install(local::lib)'

# Download latest bioperl
git clone https://github.com/bioperl/bioperl-live.git
cd bioperl-live
perl Build.PL
./Build test 

# biuld test wrong?? then test them one by on, such as 
./Build test --test_files t/LocalDB/Taxonomy/sqlite.t --verbose
./Build test --test_files t/Root/RootIO.t --verbose

# install the failed library according bulid test result
perl -MCPAN -Mlocal::lib -e 'CPAN::install(DBI)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(DBD::SQLite)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(LWP::UserAgent)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(LWP::UserAgent)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::Phylo)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::DB::Sam)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Graph::Directed)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Twig)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::Ext::Align)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Parser::PerlSAX)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Simple)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(HTML::TableExtract)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(IO::Scalar)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::SAX)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::SeqIO::staden::read)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Parser::PerlSAX)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::SeqIO::staden::read)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::LibXML)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Convert::Binary::C)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::SAX::Writer)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Writer)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::SeqIO::staden::read)'

perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::DB::Sam)'

Bio::DB::Sam
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::DB::Sam)'

# test again
./Build test 
## Test bioperl installation (you should get a version number)
perl -MBio::Root::Version -le 'print $Bio::Root::Version::VERSION'



wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
tar xjf samtools-0.1.18.tar.bz2 && cd samtools-0.1.18
make CFLAGS=-fPIC
export SAMTOOLS=`pwd`
cpanm Bio::DB::Sam




# 2016-05-26

cd /home/shg047/oasis/SALK/bam
# only call SNP
samtools mpileup --skip-indels -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f <reference genome.fa> -o <output.bcf> <list of input bam files>
bcftools index  <output.bcf> <indexed.bcf>
bcftools call --skip-variants indels --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>
# only call Indel
samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f <reference genome.fa> -o <output.bcf> <list of input bam files>
bcftools index  <output.bcf> <indexed.bcf>
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>
# call indel and SNP
samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f <reference genome.fa> -o <output.bcf> <list of input bam files>
bcftools index  <output.bcf> <indexed.bcf>
bcftools call --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>

samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,AD,ADF,SP -f /home/shg047/oasis/db/hg19/hg19.fa -o STL001RV-01.chr6.sorted.clipped.bcf  STL001RV-01.chr6.sorted.clipped.bam
bcftools index  STL001RV-01.chr6.sorted.clipped.bcf  
bcftools call --skip-variants indels --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>

bcftools call --multiallelic-caller --variants-only  -O v STL001RV-01.chr6.sorted.clipped.bcf  -o STL001RV-01.chr6.sorted.clipped.vcf
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v STL001RV-01.chr6.sorted.clipped.bcf  -o STL001RV-01.chr6.sorted.clipped.vcf

cd  /oasis/tscc/scratch/ddiep/Working/160428_RapidRun
You must learn the pipeline for the methylation calling

rm * chrLambdaNEB.methylFreq

cd /home/shg047/oasis/monod/hapinfo
R CMD BATCH trim.R

/home/shg047/oasis/monod/HMHPlasma/Excl/excl.txt
742746b3f183
xse251<-unique(xse)
match(xse251,xse)

Writing reviews of academic papers: https://github.com/jtleek/reviews
A guide to reading scientific papers: GitHub - jtleek/readingpapers: 
Good Habit for Bioinformatics Analyst or Scientist: https://www.biostars.org/p/190366/
Lab notebook for bioinformatics/data analysis work: https://www.biostars.org/p/191984/
The Most Common Stupid Mistakes In Bioinformatics: https://www.biostars.org/p/7126/#191667

cp ../Bladder.hapInfo.txt ./ 
cp ../Brain.hapInfo.txt ./ 
cp ../Colon.hapInfo.txt ./ 
cp ../Esophagus.hapInfo.txt ./ 
cp ../Intestine.hapInfo.txt ./ 
cp ../Kidney.hapInfo.txt ./ 
cp ../Liver.hapInfo.txt ./ 
cp ../Lung.hapInfo.txt ./ 
cp ../Pancreas.hapInfo.txt ./ 
cp ../Stomach.hapInfo.txt ./ 
cp ../NC-P.hapInfo.txt ./ 
cp ../WB.hapInfo.txt ./ 


#!/usr/bin/perl
use CWD;
my $dir=getcwd;
my @file=glob("NC-P*hapInfo.txt");
foreach my $file(@file){
my ($sam,undef)=split/\./,$file;
open OUT,">$file.job";
print OUT "#!/bin/csh\n";
print OUT "#PBS -n $file\n";
print OUT "#PBS -q hotel\n";
print OUT "#PBS -l nodes=1:ppn=1\n";
print OUT "#PBS -l walltime=1:00:00\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shicheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $dir\n";
print OUT "perl HMHInPlasmaTest.pl $file excl.txt $sam.HMH\n";
}

system("wget https://cran.r-project.org/src/contrib/gdata_2.17.0.tar.gz")
system("wget https://cran.r-project.org/src/contrib/caTools_1.17.1.tar.gz")

Step 1. Download Roadmap Data
cd /home/shg047/oasis/Roadmap
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/EG.mnemonics.name.xls
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/FractionalMethylation.tar.gz
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/FractionalMethylation.tar.gz.md5sum
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/header
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/ReadCoverage.tar.gz
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/ReadCoverage.tar.gz.md5sum

Step 1. Repalce -1 to space, so that perl could calculate the mean and SD based on that table
for i in {1..22} "X" "Y" "M"
do
perl -p -i -e "s/-1/ /g" chr$i.fm &
done

Step 2. Calculate the Regionss mean methylation level and Summary SD

perl MethSearchRoadmap.pl target.txt ESCA-space.fm.rlt.txt 23 2 &
perl MethSearchRoadMapByLoc.pl target.txt Target-Full-Tissue-space.fm.rlt.txt 2 &

wget http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/DMRs/WGBS_DMRs_v2.tsv.gz

for i in {1..5}
do
echo "perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl 6-P-1.hapInfo.txt 6-T-1.hapInfo.txt excl.txt 6-P-1.HMH" | qsub -q hotel
echo "perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl 7-P-1.hapInfo.txt 7-T-1.hapInfo.txt excl.txt 7-P-1.HMH" | qsub -q hotel
echo "perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl PC-P-1.hapInfo.txt PC-T-1.hapInfo.txt excl.txt PC-P-1.HMH" | qsub -q hotel
done

 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 
perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl 6-P-1.hapInfo.txt 6-T-1.hapInfo.txt excl.txt 6-P-1.HMH
perl ./HMHInPlasmaTest.pl 7-P-1.hapInfo.txt 7-T-1.hapInfo.txt excl.txt 7-P-1.HMH
perl ./HMHInPlasmaTest.pl PC-P-1.hapInfo.txt PC-T-1.hapInfo.txt excl.txt PC-P-1.HMH
perl ./HMHInPlasmaTest.pl 6-P-2.hapInfo.txt 6-T-2.hapInfo.txt excl.txt 6-P-2.HMH
perl ./HMHInPlasmaTest.pl 7-P-2.hapInfo.txt 7-T-2.hapInfo.txt excl.txt 7-P-2.HMH
perl ./HMHInPlasmaTest.pl PC-P-2.hapInfo.txt PC-T-2.hapInfo.txt excl.txt PC-P-2.HMH
perl ./HMHInPlasmaTest.pl 6-P-3.hapInfo.txt 6-T-3.hapInfo.txt excl.txt 6-P-3.HMH
perl ./HMHInPlasmaTest.pl 7-P-3.hapInfo.txt 7-T-3.hapInfo.txt excl.txt 7-P-3.HMH
perl ./HMHInPlasmaTest.pl PC-P-3.hapInfo.txt PC-T-3.hapInfo.txt excl.txt PC-P-3.HMH
perl ./HMHInPlasmaTest.pl 6-P-4.hapInfo.txt 6-T-4.hapInfo.txt excl.txt 6-P-4.HMH
perl ./HMHInPlasmaTest.pl 7-P-4.hapInfo.txt 7-T-4.hapInfo.txt excl.txt 7-P-4.HMH
perl ./HMHInPlasmaTest.pl PC-P-4.hapInfo.txt PC-T-4.hapInfo.txt excl.txt PC-P-4.HMH
perl ./HMHInPlasmaTest.pl 6-P-5.hapInfo.txt 6-T-5.hapInfo.txt excl.txt 6-P-5.HMH
perl ./HMHInPlasmaTest.pl 7-P-5.hapInfo.txt 7-T-5.hapInfo.txt excl.txt 7-P-5.HMH
perl ./HMHInPlasmaTest.pl PC-P-5.hapInfo.txt PC-T-5.hapInfo.txt excl.txt PC-P-5.HMH


#!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/Batch3Plasma/hapinfo
 perl ~/bin/hapinfo2mf.pl ./ > batch3AMF.txt
 perl ~/bin/hapinfo2mhl.pl  Hapinfo_File_list > batch3MHL.txt

qsub NC-P-43.*job
qsub CRC-P-3.*job
qsub CRC-P-10.*job
qsub NC-P-47.*job



cd  /oasis/tscc/scratch/ddiep/Plasma_10ngRRBS/BAMfiles
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl > /home/shg047/oasis/monod/Batch3Plasma/saminfo.txt
cd /home/shg047/oasis/monod/Batch3Plasma/
cp *bam *bai /home/shg047/oasis/monod/Batch3Plasma/bam

perl ~/bin/bam2hapInfo2PBS.pl ../saminfo.txt submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt


7-P-1.hapInfo.txt
a='80981116-80981175'

grep $a 6-P-1.hapInfo.txt
grep $a 6-T-1.hapInfo.txt
grep $a NC-P.hapInfo.txt



cp ../hapinfo/6-P-1.hapInfo.txt ./
cp ../hapinfo/6-P-2.hapInfo.txt ./
cp ../hapinfo/6-P-3.hapInfo.txt ./
cp ../hapinfo/6-P-4.hapInfo.txt ./
cp ../hapinfo/6-P-5.hapInfo.txt ./
cp ../hapinfo/7-P-1.hapInfo.txt ./
cp ../hapinfo/7-P-2.hapInfo.txt ./
cp ../hapinfo/7-P-3.hapInfo.txt ./
cp ../hapinfo/7-P-4.hapInfo.txt ./
cp ../hapinfo/7-P-5.hapInfo.txt ./
cp ../hapinfo/PC-P-1.hapInfo.txt ./
cp ../hapinfo/PC-P-2.hapInfo.txt ./
cp ../hapinfo/PC-P-3.hapInfo.txt ./
cp ../hapinfo/PC-P-4.hapInfo.txt ./
cp ../hapinfo/PC-P-5.hapInfo.txt ./

perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-2.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-3.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1


report.CanHMH.pl

/home/shg047/bak/plink/china

cd /oasis/tscc/scratch/shg047/monod/hapinfo
head -n 2000 6-P-1.hapInfo.txt > /oasis/tscc/scratch/shg047/monod/test/6-P-1.hapInfo.txt
head -n 2000 6-P-1.hapInfo.txt > /oasis/tscc/scratch/shg047/monod/test/6-T-1.hapInfo.txt

cat NC-P-*.hapInfo.txt >> /oasis/tscc/scratch/shg047/monod/test/NC-P.hapInfo.txt
cd /oasis/tscc/scratch/shg047/monod/test/


wig2bed  --zero-indexed  < wgEncodeBroadHistoneGm12878H3k4me1StdSig.wig  > wgEncodeBroadHistoneGm12878H3k4me1StdSig.bed
bedtools cluster -d 50 -i wgEncodeBroadHistoneGm12878H3k4me1StdSig.wig > wgEncodeBroadHistoneGm12878H3k4me1StdSig.bed.bed &
./PileOMeth extract -r chr1:723205-725116 --minDepth 1 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped
./PileOMethG extract -r chr1:723205-725116 --minDepth 1 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped

SIEMENS
awk '{if ($0 ~ /^>/) { print $0; } else { printf("%s%c%s\n", substr($0, 1, 9), "X", substr($0, 11, length($0) - 10))}}' in.fa > out.fa
  
/home/shg047/oasis/db/mm9/chr10.fa
samtools tview Indx01.merged.bam.sorted.bam /home/shg047/oasis/db/mm9/chr10.fa 

chr10:100222302-100223142


samtools bedcov  /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  /home/shg047/oasis/monod/bam/RRBS2/RRBS2/RRBS-6P16.sorted.clipped.bam 


grep chr10:10003965-10004290 B6ES.hapInfo_fixed.txt

grep chr10:100222148-100223142 2A4F1_miPS.hapInfo_fixed.txt

samtools tview 

chr10:100222148-100223142

#!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /oasis/tscc/scratch/shg047/monod/test
 perl report_cancerHM.pl ../hapinfo/6-P-1.hapInfo.txt ../hapinfo/6-T-1.hapInfo.txt ../hapinfo/NC-P-1.hapInfo.txt
 
 
 cd /home/shg047/bak/plink/phase  
 phase -f0 -F0.1 -p0.7 -c test.inp test.out


 cd /oasis/tscc/scratch/zhl002/hapinfo
 perl ~/bin/hapinfo2mhl.pl Hapinfo_File_List.txt > ~/MHL.output.12april.txt

 cd /oasis/tscc/scratch/zhl002/hapinfo

 grep chr10:100017509-100017813 B6ES.hapInfo_fixed.txt
 

2016-04-11

grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-T-2*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-T-3*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-P-9*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-T-1*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/DennisLo2015/hapinfo/CTR118*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/STL003SG-01*hapInfo.txt

7-T-1
torrow:

1) transfer 6 column file to 4 column file
2) try to upload to UCSC
3, try to find what happend to the plasma data?? any thing wrong? and check the row data



 #!/usr/bin/perl
 foreach my $file (glob("*bw")){
 my ($name,undef)=split /\_/,$file;
 print "track type=bigWig name=\"$name\" description=\"$file\" bigDataUrl\=ftp\:\/\/ucsd002\:ucsd002\@132.239.189.199\/wig\/$file\n";
 }


# Visulization of MONOD dataset
awk '!/track/ {print $1,$2,$3,$4}' 7-P-1_CpG.bedGraph | sort -k1,1 -k2,2n - > 7-P-1_CpG.bedGraph.Sort.V4
bedGraphToBigWig 7-P-1_CpG.bedGraph.Sort.V4 hg19.chrom.sizes 7-P-1.bw
curl -T 7-P-1.bw ftp://132.239.189.199 --user ucsd002:ucsd002


for i in `ls *bedGraph`
do
awk '!/track/ {print $1,$2,$3,$4}' $i | sort -k1,1 -k2,2n - > $i.S4
done

for i in `ls *.SortV4`
do
bedGraphToBigWig $i hg19.chrom.sizes $i.bw
done

for i in `ls *bw`
do
curl -T $i ftp://132.239.189.199 --user ucsd002:ucsd002
done


foreach my $file (glob("*bw")){
my ($name,undef)=split /\_/,$file;
print "track type=bigWig name=\"$name\" description=\"file\" bigDataUrl=ftp://ucsd002:ucsd002@132.239.189.199/wig/"
}

# utilities tools
bedtools sort -header -i 7-P-1_CpG.bedGraph > 7-P-1_CpG.bedGraph.sort
curl -T 7-P-1.bw ftp://132.239.189.199 --user ucsd002:ucsd002
find mydir -type f -exec curl -u xxx:psw --ftp-create-dirs -T {} ftp://192.168.1.158/public/demon_test/{} \;
track type=bigWig name="My" description="Lab" bigDataUrl=ftp://ucsd002:ucsd002@132.239.189.199/wig/7-P-1.bw
scp 7-P-1.bw shicheng@meangenemachine.dynamic.ucsd.edu:
@132.239.189.199/wig/7-P-1.bw

ftp 132.239.189.199
132.239.189.199/home/ucsd002/wig
samtools tview /home/shg047/oasis/monod/bam/RRBS1/bam/6-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr17:75283831-75283978
chr17:75283831-75283978 6-P-5.hapInfo.txt
chr17:75283978

setwd("/oasis/tscc/scratch/shg047/monod/mhl")
load("MONOD-Apr6.MHL.RData")
x1<-grep("chr19:58220439-58220459",rownames(data))
x2<-grep("chr19:58220481-58220515",rownames(data))
x3<-grep("chr19:58220626-58220668",rownames(data))
data[x1,]
data[x2,]
data[x3,]

[File:Pancreatic.cancer.plasma.vs.normal.plasma.significant.txt]
[[File:Lung.cancer.plasma.vs.normal.plasma.significant.txt]]
[[File:Colon.cancer.plasma.vs.normal.plasma.significant.txt]]

[[File:13E8.tm-colon.png]]
[[File:146E.tm-lung.png]]
[[File:1508.tm-pancreatic.png]]

SRX381713_normal_lung

less -S monod.mhl.List1.march30.txt

chr10:100027918-100027944

# 1
grep chr10:101302427-101302727 ../hapinfo/6-P-1.hapInfo.txt
grep chr10:100227698-100227747 ../hapinfo/SRX381713_normal_lung.
chr10:100095550-100096429
1E12P20_miPS 0.555555

/home/shg047/oasis/db/mm9/chr10.fa
samtools tview /home/shg047/oasis/db/mm9/chr10.fa

samtools bedcov  /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  /home/shg047/oasis/monod/bam/RRBS2/RRBS2/RRBS-6P16.sorted.clipped.bam 

echo "samtools bedcov /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles/RRBS-6P17.sorted.clipped.bam" | qsub -q glean -N test
How to install Perl package
perl -MCPAN -e "install Getopt::Long"
perl -MCPAN -e "use Getopt::Long"

perl GeneSearch.pl -i pancreatic.ncbi.txt -g GeneSymbolList.txt
/home/shg047/perl5/perlbrew/Getopt-Long-2.48/lib/Getopt/Long.pm
 
bedtools intersect -wao -a colon.cancer.plasma.vs.normal.plasma.significant.txt -b hg19_refGene.bed > colon.cancer.sig.mhl.sorted.Annotated.txt
awk '{print $1}' most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.txt > most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed
perl -p -i -e "s/[:-]/\t/g" most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed 
bedtools intersect -wao -a most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed -b hg19_refGene.bed > most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed.txt.Annotated.txt
awk '{print $8}' most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed.txt.Annotated.txt > most.signficant.Colon.Gene.list.txt

[[758D.tm.png|400px]]
[[1170.tm-lung-cancer-plasma.png|400px]]
[[3D4E.tm.-pancreatic-cancer-plasma.png|400px]]

chr10:123-456
samtools tview 
/home/shg047/oasis/DennisLo2015/sortbam/BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.PileOMeth.MF.txt

file="lung.cancer.plasma.vs.normal.plasma.significant.txt"
awk '{print $1}' $file > $file.bed
perl -p -i -e "s/[:-]/\t/g" $file.bed
bedtools intersect -wao -a $file.bed -b hg19_refGene.bed > $file.Annotated.txt
awk '{print $8}' $file.Annotated.txt | sort -u > $file.Gene.list.txt

d1<-read.table("xx.lung.target.txt",sep="\t")
d2<-read.table("lung.cancer.plasma.vs.normal.plasma.significant.txt.Gene.list.txt",sep="\t")
d2[na.omit(match(d1[,1],d2[,1])),]

file="pancreatic.cancer.plasma.vs.normal.plasma.significant.txt"
awk 'NR !=1 {print $1}' $file > $file.bed
perl -p -i -e "s/[:-]/\t/g" $file.bed
bedtools intersect -wao -a $file.bed -b hg19_refGene.bed > $file.Annotated.txt
awk '{print $8}' $file.Annotated.txt | sort -u > $file.Gene.list.txt

cd /home/shg047/oasis/DennisLo2015/mhl
awk '{print $1}' dennis.hapinfo2mhl.march28.txt > Dennis.bed
perl -p -i -e "s/[:-]/\t/g" Dennis.bed
/home/shg047/oasis/DennisLo2015/mhl/Dennis.bed

bedtools intersect -wb 

#!/usr/bin/perl
use Cwd;
use strict;
my $file="/home/shg047/oasis/DennisLo2015/mhl/Dennis.bed";
open F,$file;
while(<F>){
chomp;
system("PileOMeth extract -p 5 -q 10 --minDepth 1 -r $_ /home/shg047/oasis/db/hg19/hg19.fa BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam -o Tmp");
my $amf=qx/grep -v track Tmp_CpG.bedGraph | awk '{cov+=$5+\$6;nC+=\$5}END{print nC\/cov}'/;
print "$_\t$amf";
}

less -S ../mhl/dennis.hapinfo2mhl.march28.txt
PileOMeth extract -p 5 -q 10 --minDepth 1 -r chr10:120967293-120967420 /home/shg047/oasis/db/hg19/hg19.fa BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam -o Tmp2
grep -v track Tmp_CpG.bedGraph | awk '{cov+=$5+$6;nC+=$5}END{print nC/cov}' 

PileOMeth extract -p 5 -q 10 --minDepth 1 -r chr10:100069336-100069468 /home/shg047/oasis/db/hg19/hg19.fa BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam -o Tmp
grep -v track Tmp_CpG.bedGraph | awk '{cov+=$5+$6;nC+=$5}END{print nC/cov}' 


chr10:3480-6970
chr10:12345-12347

PileOMeth extract -r chr10:123-456 genome.fa alignments.bam
PileOMeth extract -l a.txt genome.fa alignments.bam  >  c.txt


cd /home/shg047/oasis/monod/bam/RRBS1
perl PileOMethPBS.pl
cd /home/shg047/oasis/monod/bam/RRBS2
perl PileOMethPBS.pl

PileOMeth extract --minDepth 10 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped
PileOMeth extract --counts -p 5 -q 10 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
gunzip chr.fa.gz
samtools view -bh RRBS-6P27.sorted.clipped.bam chr1:723205-726280 > Hanne.test.bam
samtools sort -o Hanne.test.sort.bam Hanne.test.bam
samtools index Hanne.test.sort.bam
PileOMeth extract -r chr1:723205-725116  chr1.fa  RRBS-6P27.sorted.clipped.bam  -o y1.ouput.txt
PileOMeth extract -r chr1:725240-726280  chr1.fa  Hanne.test.sort.bam  -o y2.ouput.txt
PileOMeth extract -l interval.input.txt chr1.fa  Hanne.test.sort.bam  -o y1y2.ouput.txt



grep -v track PileOMeth.output | awk '{cov+=$5+$6;nC+=$5}END{print nC/cov}'

samtools tview /home/shg047/oasis/monod/bam/RRBS1/bam/NC-P-12.sorted.clipped.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr1:725240-725280

NC-P-12.sorqted.clipped_CpG.bedGraph

chr1    725240  725241  0       0       35
chr1    725241  725242  100     1       0
chr1    725255  725256  0       0       36
chr1    725256  725257  100     1       0
chr1    725269  725270  0       0       45
chr1    725270  725271  100     1       0
chr1    725274  725275  4       2       47
chr1    725275  725276  100     1       0
chr1    725279  725280  5       1       16

Rscript --vanilla ~/bin/matrixPlot.R -i "monond-colon-plasma.txt" -o "monond-colon-plasma"
bed="/home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed"
fa="RRBS-6P11 /home/shg047/oasis/db/hg19/hg19.fa"
PileOMeth extract --fraction -q 10 -p 5 -l $bed -o  RRBS-6P11 $fa RRBS-6P11.sorted.clipped.bam  
PileOMeth extract -l /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  RRBS-6P11 /home/shg047/oasis/db/hg19/hg19.fa RRBS-7P30.sorted.clipped.bam
PileOMeth extract -l /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /home/shg047/oasis/db/hg19/hg19.fa RRBS-7P30.sorted.clipped.bam

 
head /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed > xx.bed



PileOMeth extract /home/shg047/oasis/db/hg19/hg19.fa RRBS-7P30.sorted.clipped.bam


cp /home/shg047/oasis/Tumor-WGBS/HCT116-SRX669642/mergeHapinfo/* ~/oasis/monod/hapinfo
cp /home/shg047/oasis/Tumor-WGBS/HCC-SRX332736/mergeHapinfo/* ~/oasis/monod/hapinfo
cd ~/oasis/monod/hapinfo


perl ~/bin/samInfoPrep4Bam2Hapinfo.pl  .sorted.clipped.bam > ../Saminfo4bam2hapinfo.txt
perl ~/bin/bam2hapInfo2PBS.pl ../Saminfo4bam2hapinfo.txt submit nonbismark



 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=4:00:00
 #PBS -o hapinfo2mhl.o
 #PBS -e hapinfo2mhl.e
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/SALK/mergeHapinfo
 perl ~/bin/hapinfo2mhl.pl Salk.HapinfoList > Salk.hapinfo2mhl.march29.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Salk.hapinfo2mhl.march28.txt" -o "Salk.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > Salk.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Salk.hapinfo2mf.march28.txt" -o "Salk.hapinfo2mf.march28"
 
  Rscript --vanilla ~/bin/matrixPlot.R -i "monod.mhl.List1.march29.txt" -o "monod.mhl.List1.march29"

 

 
 
head -n 50 /home/shg047/oasis/AgeBlood/hapinfo/WB-centenarian.hapInfo.txt
head -n 50 /home/shg047/oasis/monod/hapinfo/WB-centenarian.hapInfo.txt

#!/usr/bin/perl
use strict;
use Cwd;

my $dir1="/home/shg047/oasis/Estellar2016/mergeHapinfo/";
chdir $dir1;
my @hap1=glob("$dir1/*lung*hapInfo.txt");
my @hap2=glob("$dir1/*colon*hapInfo.txt");

my $dir2="/home/shg047/oasis/DennisLo2015/hapinfo"
chdir $dir2;
my @hap3=glob("CTR*hapInfo.txt");
my @hap4=glob("Pregnancy*hapInfo.txt");

my $dir3="/home/shg047/oasis/monod/hapinfo";
chdir $dir3;
my @hap5=glob("*hapInfo.txt");

my @file=(@hap1,@hap2,@hap3,@hap4,@hap5);
foreach my $file(@file){
print "$file\n";
}

 
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mhl.march28.txt" -o "Estellar2016.hapinfo2mhl.march28"

 



LC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt
LC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt
LC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt
LN	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381713_normal_lung.hapInfo.txt
NC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381553_normal_colon.hapInfo.txt
CC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt
CC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381585_metastasis_colon.hapInfo.txt


/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/



cp /oasis/tscc/scratch/ddiep/Ziller_BAMfiles/Colon_Tumor_Primary* ./

#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -o wget.log
#PBS -e wget.err
#PBS -M diep.hue.dinh@gmail.com
#PBS -m abe
cd /home/ddiep/dinh_working/Tumor_WGBS
#colon tumor primary tissue
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949210/SRR949210.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949211/SRR949211.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949212/SRR949212.sra
#hct116
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR153/SRR1536575/SRR1536575.sra

GSE46644	GSM1204465	SRX332736	Colon primary tumor
GSE60106	GSM1465024	SRX669642	HCT116

GSE16256
GSE17312
GSE31971
GSE30340


Ziller et al paper.CD14 CD56 CD19

#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -q hotel
#PBS -l walltime=12:00:00
#PBS -o wget.log
#PBS -e wget.err
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
# Colon_Tumor_Primary:HCC-SRX332736
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949210/SRR949210_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949210/SRR949210_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949211/SRR949211_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949211/SRR949211_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949212/SRR949212_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949212/SRR949212_2.fastq.gz &

# HCT116-SRX669642
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR1536575/SRR1536575_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR1536575/SRR1536575_2.fastq.gz &

cp /home/shg047/oasis/SALK/mergeHapinfo/STL*.txt ./
wc -l STL003AO-01.hapInfo.txt

# 2016-03-28
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279516/GSM1279516_CpGcontext.Brain_W.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279517/GSM1279517_CpGcontext.Breast.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279518/GSM1279518_CpGcontext.CD19.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279519/GSM1279519_CpGcontext.Colon.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279520/GSM1279520_CpGcontext.Colon_M.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279521/GSM1279521_CpGcontext.Colon_P.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279522/GSM1279522_CpGcontext.H1437.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279523/GSM1279523_CpGcontext.H157.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279524/GSM1279524_CpGcontext.H1672.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279527/GSM1279527_CpGcontext.Lung.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279532/GSM1279532_CpGcontext.U87MG.txt.gz


gzip -frtv9 *

/home/shg047/oasis/Estellar2016/mergeBam


GSE52271


library("GEOquery")
GEOSet <- getGEO("GSE56851")
data <- as.data.frame(exprs(GEOSet[[1]]))
phen <- pData(phenoData(GEOSet[[1]]))

 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=4:00:00
 #PBS -o hapinfo2mhl.o
 #PBS -e hapinfo2mhl.e
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/mhl
 Rscript --vanilla ~/bin/matrixPlot.R -i "monod.mhl.march29-List1.txt" -o "monod.mhl.march29-List1"
 cd /home/shg047/oasis/Estellar2016/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > Estellar2016.hapinfo2mhl.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mhl.march28.txt" -o "Estellar2016.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > Estellar2016.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mf.march28.txt" -o "Estellar2016.hapinfo2mf.march28"

 
  
 cd /home/shg047/oasis/monod/mhl
 perl ~/bin/hapinfo2mhl.pl hapinfo2mhl.sampleList > monod.mhl.march29.txt
 
 vim hapinfo2mhl-march29.job
 
 
 cd /home/shg047/oasis/Estellar2016/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > Estellar2016.hapinfo2mhl.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mhl.march28.txt" -o "Estellar2016.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > Estellar2016.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mf.march28.txt" -o "Estellar2016.hapinfo2mf.march28"

 
 
 
cd /home/shg047/oasis/monod/hapinfo/hapinfo
cd /home/shg047/oasis/monod/hapinfo/WGBS





[[File:66F0.tm-mhl-dennislo.png]]
[[File:6681.tm.png]]

wc -l dennis.hapinfo2mf.march28.txt
wc -l dennis.hapinfo2mhl.march28.txt

 d1<-read.table("dennis.hapinfo2mf.march28.txt",head=T,row.names=1)
 d2<-read.table("dennis.hapinfo2mhl.march28.txt",head=T,row.names=1)
 library("ggplot2")
 library(RColorBrewer)
 require(KernSmooth)
 pdf("comp-1.smoothscatter.relationship.pdf")
 par(mfrow=c(5,5),mar=c(1,1,1,1))
 for(i in 1:25){
 g = 11
 n=length(d1[,i])
 my.cols <- rev(brewer.pal(g, "RdYlBu"))
 smoothScatter(d1[,i], d2[,i], nrpoints=0.05*n, colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1",xlab="",ylab="")
 }
 dev.off()
 pdf("comp-2.smoothscatter.relationship.pdf")
 par(mfrow=c(2,2),mar=c(2,2,2,2))
 hist(d1[,1],breaks=10,main="MF",col="blue")
 hist(d2[,1],breaks=10,main="MHL",col="blue")
 plot(density(d1[,1],na.rm=T),main="MF",col="blue",lwd=2)
 plot(density(d2[,1],na.rm=T),main="MHL",col="blue",lwd=2)
 dev.off()




# 2016-03-27
vim ~/bin/samInfoPrep4Bam2Hapinfo.pl 
Data<-cbind(GSE53045NormalPBMC,GSE35069NormalPBMC)
cd /home/shg047/oasis/DennisLo2015/sortbam
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ./ > ../Saminfo4bam2hapinfo.txt
perl ~/bin/bam2hapInfo2PBS.pl ../Saminfo4bam2hapinfo.txt submit bismark

vim ~/bin/samInfoPrep4Bam2Hapinfo.pl
cd /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ./ > /home/shg047/oasis/monod/hapinfo/phase2/Saminfo4bam2hapinfo.txt
less -S /home/shg047/oasis/monod/hapinfo/phase2/Saminfo4bam2hapinfo.txt

cd /home/shg047/oasis/monod/hapinfo/phase2/
perl ~/bin/bam2hapInfo2PBS.pl Saminfo4bam2hapinfo.txt submit nonbismark

vim ~/bin/bam2hapInfo2PBS.pl



# 2016-03-26


mkdir test2
samtools view -h /home/shg047/oasis/DennisLo2015/sortbam/T21.1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam chr10:31891611-31891781 > x.sam
samtools view -bh x.sam > x.bam
samtools sort -o x.sort.bam x.bam
samtools index x.sort.bam

samtools tview x.sort.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr10:31891611-31891781
/home/shg047/bin/mergedBam2hapInfo.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /oasis/tscc/scratch/shg047/DennisLo2015/sortbam/test2x.sort.bam bismark > ../hapinfo/x.sort.bam.hapInfo.txt

samtools sort -n -o x.nsort.bam x.bam 
bismark_methylation_extractor --bedGraph --zero_based --comprehensive --cutoff 1  --mbias_off --paired-end x.nsort.bam

PileOMeth extract -q 10 -p 5 --mergeContext ~/oasis/db/hg19/hg19.fa x.sort.bam 


/home/shg047/bin/mergedBam2hapInfo.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /home/shg047/oasis/DennisLo2015/bam/test/test.bam bismark > test.read1.hapInfo.txt
1, Duplicated fragments/reads removing(BS-Seq suitable. RRBS not suitable)
2, Low sequencing quality reads removing (Phred >= 5)
3, Low mapping quality reads removing (MAPQ >= 10)

cd /home/shg047/software/IGV_2.3.68

samtools view -h T21.1.read1_val_1.fq.gz_bismark_bt2_pe.bam | head -n 30001 > x.sam

cd /home/shg047/oasis/DennisLo2015/bam/test
bismark_methylation_extractor --bedGraph --zero_based --comprehensive --cutoff 1  --mbias_off --paired-end test.bam

samtools sort -o test.sort.bam test.bam 
samtools index test.sort.bam
PileOMeth extract --fraction --mergeContext ~/oasis/db/hg19/hg19.fa test.sort.bam 

head test.sort_CpG.bedGraph
head test.bedGraph.gz.bismark.zero.cov

bedtools sort -i test.sort_CpG.bedGraph > test.sort_CpG.sort.bedGraph
bedtools sort -i test.bedGraph.gz.bismark.zero.cov > test.sort.bedGraph.gz.bismark.zero.cov

wc -l test.sort_CpG.sort.bedGraph
wc -l test.sort.bedGraph.gz.bismark.zero.cov

head test.sort_CpG.sort.bedGraph
head test.sort.bedGraph.gz.bismark.zero.cov

bedtools intersect -wo -a test.sort_CpG.sort.bedGraph -b test.sort.bedGraph.gz.bismark.zero.cov > merge.bed

data<-read.table("merge.bed")
pdf("comp.smoothscatter.relationship.pdf")
par(mfrow=c(3,3))
for(i in 1:1){
g = 11
n=length(data[,i])
my.cols <- rev(brewer.pal(g, "RdYlBu"))
smoothScatter(data[,4], data[,10], nrpoints=0.05*n, colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1",xlab="Methylation frequency",ylab="MHL")
}
dev.off()

 chr1:121485049-121485051
 
 
 samtools tview test.sort.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr1:121485049-121485051



Indx16_S12.sorted.clipped.bam  chr10:102443911-102443400


#!/usr/bin/perl

my $a="0400";
my $a2=0400;
my $b=16;
my $c=$a-$b;
my $c2=$a2-$b;
my $d=0x1a;

if($a){
	print "$a\t$b\t$c\t$a2\t$c2\t$d\n";
}



# 2016-03-24
1, transfer file to new hard-disk:  cp -r /oasis/* /media/NAS2_volume1/shg047
2, fastx_trimmer -Q33 -f 8 -l 43 -i Temp/P608-Tumor.file1.fastq -o Temp/P608-Tumor.file1.trimmed.fastq

qsub -q hotel SRR1035882_1_val_1.fq.gz.job
qsub -q hotel SRR1035895_1_val_1.fq.gz.job

/media/LTS_60T/Dinh

find /tmp -type f -size +50000k -delete



20G     ./bam
1.0K    ./biomarker
9.3G    ./tcga
166G    ./db
33G     ./meth450
1.9G    ./GEO
630G    ./Estellar2016
2.0G    ./mice
398M    ./wbc
911K    ./alice
733M    ./blue
25M     ./haib
2.9G    ./reprogramming
55M     ./ssc
438M    ./bak
422G    ./monod
44G     ./song2
1015M   ./twin
1.3T    ./

rm ../bam/*bam &



# 2016-03-23
cd 

BMT2.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err
CTR101_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR103_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR104_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR85_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
HOT240_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err


bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov.sort
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.sort.bedGraph.gz.bismark.zero.cov 
awk 'NR>8695750 && NR<8695760' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx.bed
4898107
awk '{print $1,$2,$3}' OFS="\t" CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > yy.bed
bedtools sort -i yy.bed > yy.sort.bed
awk '{a[$1];}' xx1.bed | head
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx.bed
awk '{print FILENAME, NR, FNR, $0}' file1 file2





awk 'NR<4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx1.bed &
awk 'NR>=4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx2.bed &
bedtools sort -i xx1.bed > x1.sort.bed &
bedtools sort -i xx2.bed > x2.sort.bed &
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > x3.sort.bed 

awk 'NR<8695756' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx1.bed 
bedtools sort -i xx1.bed > x1.sort.bed 
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > x3.sort.bed 



awk 'NR<4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx1.bed &
awk 'NR>4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx2.bed &

data<-read.table("CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov")


9898107
wc -l CTR150_trimmed.fq.gz_bismark_bt2.sort.bedGraph.gz.bismark.zero.cov 



BMT2.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err
CTR101_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR103_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR104_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR85_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
HOT240_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
Pregnancy.1.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err
T21.2.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err



# 2016-03-22
samtools view -bh ../../sortbam/Pregnancy.7.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam chr10:2000000-2030000 -o Pregnancy.7.bam
samtools view -bh ../../bam/Pregnancy.7.read1_val_1.fq.gz_bismark_bt2_pe.bam | head -n 3000 > Pregnancy.7.bismark.bam
bismark_methylation_extractor --bedGraph --zero_based --comprehensive Pregnancy.7.csort.bam
bismark_methylation_extractor --bedGraph --zero_based --comprehensive Pregnancy.7.nsort.bam
bismark_methylation_extractor --bedGraph --zero_based --comprehensive Pregnancy.7.bismark.bam


T21.5.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam

HWI-ST1049:8:1101:1199:2048#0/1

#!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test.o
 #PBS -e test.e
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/DennisLo2015/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > dennis.hapinfo2mhl.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "dennis.hapinfo2mhl.march28.txt" -o "dennis.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > dennis.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "dennis.hapinfo2mf.march28.txt" -o "dennis.hapinfo2mf.march28"

 
 cd /home/shg047/oasis/monod/hapinfo/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > RRBS.Plasma.Phase2.MHL.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "RRBS.Plasma.Phase2.MHL.txt" -o "RRBS-Plamsa-Phase2"

 
 
 cd /home/shg047/oasis/DennisLo2015/bam/test
 perl ~/bin/bam2hapInfo2PBS.pl test.sam submit bismark > T21.1.read1_val_1.fq.gz_bismark_bt2_pe.mhl

 
 perl ~/bin/hapinfo2mhl.pl ./ submit bismark > T21.1.read1_val_1.fq.gz_bismark_bt2_pe.mhl
 
 
 
 cd /home/shg047/oasis/DennisLo2015/hapinfo
 perl ~/bin/hapinfo2mf.pl ./ > ../dennis.hapinfo2mf.march24.txt
 
 perl ~/bin/hapinfo2mhl.pl ./ > ../dennis.mhl.march24.txt

 cd /home/shg047/oasis/DennisLo2015/bam/test
bismark_methylation_extractor --bedGraph --zero_based --comprehensive CTR98.bam


# 2016-03-22
cd /home/shg047/oasis/DennisLo2015/bam
qsub -q hotel Pregnancy.4.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel Pregnancy.6.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel Pregnancy.7.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.1.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.2.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.3.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.4.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job

/home/shg047/oasis/AgeBlood/mhl/WB.mhl.march22.txt

ls -lart *bam2sortbam*err
ls -lart *bam2sortbam*err | wc -l 

cd /home/shg047/oasis/Holger2016/fastq_trim/
qsub -q hotel SRR1035734_1_val_1.fq.gz.job 
qsub -q hotel SRR1035723_1_val_1.fq.gz.job
qsub -q hotel SRR1035725_1_val_1.fq.gz.job
qsub -q hotel SRR1035741_1_val_1.fq.gz.job
qsub -q hotel SRR1035738_1_val_1.fq.gz.job
qsub -q hotel SRR1035740_1_val_1.fq.gz.job

cd /home/shg047/oasis/Holger2016/fastq_trim/
qsub -q hotel SRR1035856_1_val_1.fq.gz.job
qsub -q hotel SRR1035829_1_val_1.fq.gz.job
qsub -q hotel SRR1035727_1_val_1.fq.gz.job
qsub -q hotel SRR1035855_1_val_1.fq.gz.job
qsub -q hotel SRR1035799_1_val_1.fq.gz.job

cd /home/shg047/oasis/Holger2016/fastq_trim/
qsub -q hotel SRR1035726_1_val_1.fq.gz.job
qsub -q hotel SRR1035739_1_val_1.fq.gz.job
qsub -q hotel SRR1035735_1_val_1.fq.gz.job
qsub -q hotel SRR1035736_1_val_1.fq.gz.job
qsub -q hotel SRR1035798_1_val_1.fq.gz.job
qsub -q hotel SRR1035804_1_val_1.fq.gz.job
qsub -q hotel SRR1035827_1_val_1.fq.gz.job
qsub -q hotel SRR1035800_1_val_1.fq.gz.job

# 2016-03-21
echo "bsrate -c ~/oasis/db/hg19/hg19.fa -o xx.bsrate xx.rdup" | qsub -q glean
/oasis/tscc/scratch/ddiep/BAMfiles
Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.mhl.march19.txt" -o "Esterllar2016"
PileOMeth extract -l  ~/oasis/db/hg19/hg19.fa





/home/shg047/oasis/DennisLo2015/mr

LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o xx.sort xx.mr
duplicate-remover -S xx_stat.txt -o xx.rdup xx.sort
bsrate -c hg19 -o xx.bsrate xx.rup
grep ˆchrM xx.mr > xx.mr.chrM
bsrate -N -c chrM.fa -o xx.chrM.bsrate xx.mr.chrM
methcounts -n -c hg19 -o xx.meth xx.mr

LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o xx.sort.end xx.mr


for clinical research. clear data would be very important.
for research the machnism should be great
 
# 2016-03-18
rmapbs -c hg19 -o ../methBam/ SRR949193_1.fastq.gz SRR949193_2.fastq.gz
perl -lane 'print "chr@F[7]\t@F[8]\t@F[8]\t@F[3]\t@F[4]\t@F[6]\t@F[10]\t@F[-3]\t@F[-2]\t@F[-1]" if ! /chrNA/' TCGA.Meth450.esca.ChAMP.DMS.txt  > DMS.bed
grep -v chrNA DMS.bed > DMS2.bed
bedtools intersect -wo -a weimarch20.txt -b DMS2.bed > batch3.region.txt
 
bzip2 -d filename.bz2
This is cutadapt 1.9.1 with Python 2.6.6
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC SRR949197_1.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
gzip: SRR949197_1.fastq.gz: unexpected end of file
cutadapt: error: In read named 'SRR949197.215195112 D1JR8ACXX130107:2:2203:10919:71639 length=99': length of quality sequence (17) and length of read (99) do not match
Cutadapt terminated with exit signal: '256'.
Terminating Trim Galore run, please check error message(s) to get an idea what went wrong...

wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX332%2FSRX332731/SRR949196/SRR949196.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX332%2FSRX332731/SRR949197/SRR949197.sra
wget -r ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/ 

#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

df = read.table(opt$file, header=TRUE)
num_vars = which(sapply(df, class)=="numeric")
df_out = df[ ,num_vars]
write.table(df_out, file=opt$out, row.names=FALSE)

Rscript --vanilla yasrs.R

wget https://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz
wget https://cran.rstudio.com/src/contrib/optparse_1.3.2.tar.gz
install.packages("getopt_1.20.0.tar.gz")
install.packages("optparse_1.3.2.tar.gz")


 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 # cd /home/shg047/oasis/monod/hapinfo
 # perl ~/bin/hapinfo2mhl.pl ./ > ../../mhl.txt
 # wget -r ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/
 cd /home/shg047/oasis/Estellar2016/mergeHapinfo
 perl ~/bin/hapinfo2mf.pl ./ > Estellar2016.mf.march21.txt
 
 cd /home/shg047/oasis/DennisLo2015/mr
 bsrate -c ~/oasis/db/hg19/hg19.fa -o xx.bsrate xx.rup

 
 
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.mhl.march19.txt" -o "Esterllar2016"

 cd /home/shg047/oasis/Ziller2013/fastq
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949196_1.fastq.bz2 &
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949196_2.fastq.bz2 &
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949197_1.fastq.bz2 &
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949197_2.fastq.bz2 &

 

 perl ~/bin/hapinfo2mhl.pl ./ > ../MHL.OUTPUT.txt


/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles

RRBS-phase2： /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles


ln -s /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles/ bam




cpg_list="/media/Ext12T/DD_Ext12T/BisRef/bisHg19_plusLambda_BWA/hg19_lambda.cpg.positions.txt"
samtools mpileup -BA -f /home/shg047/oasis/db/hg19/hg19.fa $bam_file > $pileup_file
/home/dinh/scripts/BisReadMapper/src/extractMethyl.pl $cpg_list 33 < $pileup_file > $methylFreq_file



# 2016-03-17
ln -s /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles  bam
ln -s /home/k4zhang/my_oasis_tscc/MONOD/Ecker_Tissue_WGBS/BAMfiles bam


/oasis/tscc/scratch/ddiep/Working


SRR949193.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949194.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949195.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949196.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949197.fastq.download.job
-rw------- 1 shg047 k4zhang-group 537 Mar 17 10:56 SRR949198.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949199.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949200.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949201.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949202.fastq.download.job

SRR949193 SRR949202

for i in `seq 194 202`
do 
qsub SRR949$i\fastq.download.job
done






When performing an alignment one must discriminate between different types of bisulfite-treated DNA libraries. In the first, termed directional libraries, adapters are attached to the DNA fragments such that only the original
top or bottom strands will be sequenced. Alternatively, all four DNA strands that arise through bisulfite treatment and subsequent
PCR amplification can be sequenced with the same frequency in nondirectional libraries.


# 2016-03-16
-rw-r--r-- 1 shg047 k4zhang-group          62 Jan 19 12:22 Lymphoma.run3.bam
-rw-r--r-- 1 shg047 k4zhang-group          62 Jan 15 21:56 Lymphoma.run4.bam
-rw-r--r-- 1 shg047 k4zhang-group          62 Feb 10 09:51 Pregnancy.11.bam


qdel 4537516.tscc-mgr.local  shg047      hotel    Lymphoma.run3.re    --      1     16    --   72:00:00 Q       --
qdel 4537527.tscc-mgr.local  shg047      hotel    Lymphoma.run4.re    --      1     16    --   72:00:00 Q       --
qdel 4537636.tscc-mgr.local  shg047      hotel    Pregnancy.11.rea    --      1     16    --   72:00:00 Q       --

chr10:100027918-100027944	TTTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TTTC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TTCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TTCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCTC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTTC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CCTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CCCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CCCC	1	100,027,918,100,027,000,000,000,000,000,000,000

~



qsub Pregnancy.11.read1.fq.gz.job
qsub Lymphoma.run3.read1.fq.gz.job
qsub Lymphoma.run4.read1.fq.gz.job

-rw-r--r-- 1 shg047 k4zhang-group	   62 Jan 19 12:22 Lymphoma.run3.read1_val_1.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group	   62 Jan 15 21:56 Lymphoma.run4.read1_val_1.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group	   62 Feb 10 09:51 Pregnancy.11.read1_val_1.fq.gz_bismark_bt2_pe.bam



bedtools intersect -wo -a znf154.bed -b WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
chr19   58220404	 58220671	 chr19   58220439	 58220459	 chr19:58218498-58222001,B040,4  20
chr19   58220404	 58220671	 chr19   58220481	 58220515	 chr19:58218498-58222001,B042,5  34
chr19   58220404	 58220671	 chr19   58220626	 58220668	 chr19:58218498-58222001,B047,4  42


#!/usr/bin/perl
use strict;
use Cwd;
my $bam_dir=shift @ARGV;
my $bed_dir="/home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed";

chdir $bam_dir;
my @file=glob("*.bam");
foreach my $file(@file){
my ($sample,undef)=split /\./,$file;
print "$sample\t$bam_dir$file\t$bed_dir\n";
}


# 2016-03-14
/home/shg047/monod/predict/phase2
# submit 5 each time. 
qsub SRR949198fastq.download.job
qsub SRR949199fastq.download.job
qsub SRR949200fastq.download.job
qsub SRR949201fastq.download.job
qsub SRR949202fastq.download.job

qsub SRR949203fastq.download.job
qsub SRR949204fastq.download.job
qsub SRR949205fastq.download.job
qsub SRR949206fastq.download.job
qsub SRR949207fastq.download.job

qsub SRR949208fastq.download.job
qsub SRR949209fastq.download.job
qsub SRR949210fastq.download.job
qsub SRR949211fastq.download.job
qsub SRR949212fastq.download.job

qsub SRR949213fastq.download.job
qsub SRR949214fastq.download.job
qsub SRR949215fastq.download.job


# creat bisulfite treated human reference: hg19
cp hg19.fa hg19.c2t.fa
perl -p -i -e 's/CG/M/ig' hg19.c2t.fa
perl -p -i -e 's/C/T/ig' hg19.c2t.fa
perl -p -i -e 's/M/CG/ig' hg19.c2t.fa

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/ 
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:100027865
 

chr10:101089382-101089519



Sept 2014: these batch of libraries were generated with the dRRBS protocol (MSP I & Taq I digestion)
Bam folder: /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles
Mapple_bin_hapInfo folder:
Mld_block_hapInfo folder: /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo


cd /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo
head NC-P-23.mld_blocks_r2-0.5.hapInfo.txt

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles

# index: http://genome-tech.ucsd.edu/LabNotes/index.php/Dinh/Dinh_2014/NOTES/2014-9-22 
samtools tview Indx16_S12.sorted.clipped.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:102440557-102440826
samtools view Indx16_S12.sorted.clipped.bam  chr10:102443911-102443400


# /home/shg047/oasis/DennisLo2015/bam
cd /home/shg047/oasis/DennisLo2015/fastq/
qsub Lymphoma.run2.read1.fq.gz.job
qsub CTR147.fq.gz.job

cd /home/shg047/oasis/Ziller2013/fastq
qsub SRR949197*job
qsub SRR949196*job
qsub
qsub




# 2016-03-14
File:199.tm.png
File:DF.tm.png

scp /home/shg047/oasis/haib/mhl
cp /home/shg047/monod/hapinfo/STL*.hapInfo.txt ./

setwd("")
/home/shg047/monod/mixHap/hapinfo/mMHL.whole.txt


# 2016-03-11
http://genome-tech.ucsd.edu/LabNotes/index.php/Kun:LabNotes/MONOD/2015-7-6
cd /home/shg047/oasis/monod/mixHap
cut -f 1 colon.data.plsma-hypoall.txt > colon.hypermhl.plasma.txt
cut -f 1 lung.data.plsma.hypoall.txt > lung.hypermhl.plasma.txt
cut -f 1 pancreatic.data.plsma.hypoall.txt > pancreatic.hypermhl.plasma.txt

File:Colon.hypermhl.plasma.txt
File:Lung.hypermhl.plasma.txt
File:Pancreatic.hypermhl.plasma.txt

# 2016-03-10
cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/mld_blocks_stringent_hapInfo
grep 8054556 7-P-2.mld_blocks_r2-0.5.hapInfo.txt
chr17:8054556-8054584   CCCCC   1	8054556,8054568,8054579,8054581,8054584

cd /home/shg047/monod/hapinfo
grep 8054556 7-P-14.hapInfo.txt
chr17:8054556-8054584   CCCCC   1	8054556,8054568,8054579,8054581,8054584

grep 8054556 7-P-22.hapInfo.txt
chr17:8054556-8054579   CCC     1	8054556,8054568,8054579

grep 95947328 7-P-11.hapInfo.txt
chr9:95947328-95947356  CCCCC   2	95947328,95947337,95947340,95947345,95947356
chr9:95947328-95947356  TTTTT   1	95947328,95947337,95947340,95947345,95947356

grep 100204208 PC-T-2.hapInfo.txt
chr14:100204208-100204257	CCT     1	100204208,100204221,100204232

grep 100204208 PC-P-2.hapInfo.txt
chr14:100204208-100204257	CT      1	100204250,100204257





cd /home/shg047/monod/methyblock/HM450K
bedtools intersect -wa -u -a mhb450bed.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 1551
bedtools intersect -wa -u -a mhb450bed.bed -b mhbbed.bed > mhb450GWBS.bed	      # 1258
bedtools intersect -wa -u -a mhb450GWBS.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l  # 1045

cd /home/shg047/monod/methyblock/encode_rrbs_bed
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 1551
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b /home/shg047/monod/methyblock/HM450K/mhbbed.bed | wc -l    # 8920  
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b /home/shg047/monod/methyblock/HM450K/mhbbed.bed > RRBS_GWBS.bed  
bedtools intersect -wa -u -a RRBS_GWBS.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l  # 7968

cd /home/shg047/monod/methyblock/HM450K
bedtools intersect -wa -u -a mhbbed.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 79704
bedtools intersect -wa -u -a mhbbed.bed -b ~/oasis/db/hg19/CpG.Shore.hg19.bed  | wc -l  # 26103
bedtools intersect -wa -u -a mhbbed.bed -b ~/oasis/db/hg19/CpG. | wc -l   # 3246
bedtools intersect -wa -v -a mhbbed.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 3246

# 2016-03-09
/home/shg047/oasis/monod/hapinfo/WGBS



CRC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt /home/shg047/monod/heatmap
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt /home/shg047/monod/heatmap
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt /home/shg047/monod/heatmap
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt /home/shg047/monod/heatmap

2289

CRC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt

cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt ./
cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt ./
cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt ./
cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt ./

scp Heatmap.MHL.txt shg047@genome-miner.ucsd.edu:/home/shg047/monod/heatmap/
SRX381569-tumor-colon	CCT
SRX381716-adenocarcinoma-lung	LC
SRX381719-squamous-cell-tumor-lung	LC
SRX381722-small-cell-tumor-lung	LC





# 2016-03-05
mv  *mhl.in.plsma.bed /media/LTS_33T/SG_LTS33T/monod/mixHap/
source("http://bioconductor.org/biocLite.R")
biocLite("DMRcate")
sudo apt-get install iotop
sudo iotop
sudo scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/hapinfo/SRA/* /home/ucsd002/monod/hapinfo

#!/usr/bin/perl
use strict;
my @file=glob("6-*.hapInfo.txt");
my $i;
foreach my $file(@file){
	$i++;
	my ($cancer,$type,$id)=split /[.-]/g,$file;
	my $id = sprintf("%03d",$i);
	print "$zero_num\n";
	# system("cp $file CRC.$type.$id");
	print "$id\t$i\n";
}


# 2016-03-05
cd /home/shg047/monod/mixHap/mhl
N37-Lung	 Lung
STL001LG-01     Lung
STL002LG-01     Lung
STL002PA-01     Pancreas
N37-Pancreas    Pancreas
STL003PA-01     Pancreas
N37-Colon	Colon
STL001SG-01     Colon
STL003SG-01     Colon


scp shg047@genome-miner.edu:/home/shg047/monod/rrbs_kun/*hapInfo.txt ./
cd /home/shg047/oasis/monod/mhb/hapinfo
cp Colon_primary_tumor.all_chrs.hapInfo.txt HCT116.all_chrs.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/WGBS
cd  /home/shg047/oasis/monod/hapinfo/WGBS 
cd /home/shg047/oasis/monod/haplo/wbc/mergeHapinfo
mv *hapInfo.txt /home/shg047/oasis/monod/haplo/Merge
cd /home/shg047/oasis/monod/haplo/n37/mergeHapinfo
cd /home/shg047/oasis/monod/haplo/salk/mergeHapinfo
cd /home/shg047/oasis/monod/haplo/hesc/mergeHapinfo

11092 Heyn2016.R0.3.methyblock.bed
9296 Heyn2016.R0.4.methyblock.bed
6338 Heyn2016.R0.5.methyblock.bed
3082 Heyn2016.R0.6.methyblock.bed
1366 Heyn2016.R0.7.methyblock.bed
 
 # 2016-03-05
 
 #!/bin/csh
 #PBS -N hapinfo2mhl
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 Rscript MH450DMRSmoking.R
 
 
 perl ~/bin/hapinfo2mhl.pl


# 2016-03-04
Rscript PancancerGSI.R "PancancerMethMatrixaa" "header.txt" "PancancerMethSaminfo_March2016.txt" "PancancerMethMatrixaa.txt"

cd /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles
cd /oasis/tscc/scratch/ddiep/BAMfiles
perl ~/bin/SaminfoPre4hapinfo.pl > ~/oasis/Estellar2016/SaminfoPre4hapinfo.txt
cd /home/shg047/oasis/Estellar2016/hapinfo
perl ~/bin/bam2hapInfo2PBS.pl ../SaminfoPre4hapinfo.txt

#
for i in `ls *bam`
do 
ls *bam | awk -F_ '{print $1}' | sort -u
done

cd /oasis/tscc/scratch/ddiep/BAMfiles
find ./ -name '*.bam' | {
    read firstbam
    samtools view -h "$firstbam"
    while read bam; do
	 samtools view "$bam"
    done
} | samtools view -ubS - | samtools sort -o merged - 
samtools index merged.bam
ls -l merged.bam merged.bam.bai

-rwx------ 1 shg047 k4zhang-group 1631 Mar  4 09:07 SRX381553_normal_colon.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1601 Mar  4 09:07 SRX381569_tumor_colon.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1751 Mar  4 09:07 SRX381585_metastasis_colon.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1691 Mar  4 09:07 SRX381601_tumor_prostate.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1781 Mar  4 09:07 SRX381611_metastasis_breast.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1631 Mar  4 09:07 SRX381621_tumor_breast.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1661 Mar  4 09:07 SRX381631_normal_breast.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group  510 Mar  4 09:07 SRX381646_normal_prostate.bam.merge.sh
ls -larth SRX381646_*
scp shg047@genome-miner.ucsd.edu:/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles/*bam ./
scp shg047@genome-miner.ucsd.edu:/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles/*tumor*bam ./
SRX381646_normal_prostate.bam.merge.sh



cd /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles

find ./ -name 'SRX381601_tumor_prostate*.bam' | {
read firstbam
echo "$firstbam" 
} 

 samtools view -h ./SRX381601_tumor_prostate.chr14.sorted.clipped.bam

 find .|grep "FooBar"|xargs -I{} cp "{}" ~/foo/bar

 
  ls *bam | grep -v chrLambdaNEB |xargs -I{} cp "{}" /home/shg047/oasis/Estellar2016/bam

  scp shg047@genome-miner.ucsd.edu:/home/shg047/oasis/Estellar2016/bam2/*pl ./
  
 
cd /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles
find ./ -name 'SRX381725_normal_CD19*.bam' | {
read firstbam
echo "$firstbam"
while read bam; do 
echo "$bam"
done
}

samtools view ./SRX381725_normal_CD19.chr19.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chrX.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr14.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr20.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chrY.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr15.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr7.sorted.clipped.bam | head 
samtools view ./SRX381725_normal_CD19.chrM.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr11.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chrLambdaNEB.sorted.clipped.bam  | head
samtools view ./SRX381725_normal_CD19.chr10.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr9.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr17.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr3.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr21.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr12.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr22.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr18.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr4.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr5.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr13.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr16.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr8.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr2.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr6.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr1.sorted.clipped.bam | head




find ./ -name 'SRX381601_tumor_prostate*.bam' | {
read firstbam
echo "$firstbam"
samtools view -h "$firstbam"
while read bam; do
echo "$bam"
samtools view "$bam"
done
 } | samtools view -ubS - | samtools sort - /home/shg047/oasis/Estellar2016/bam/merged
samtools index /home/shg047/oasis/Estellar2016/bam/merged.bam





shg047@tscc-login.sdsc.edu


scp * shg047@tscc-login.sdsc.edu:/home/shg047/oasis/TCGA/Meth/Data


wget --ftp-user='zhang' --ftp-password='HDIIWpP' ftp://169.228.63.66/
wget -r --user='zhang' --password='HDIIWpP' ftp://169.228.63.66/

bedtools intersect -a /home/shg047/db/hg19/encode/encode.*.hg19.bed -b en.mhb.hypo.bed

rm endo.mhb.hypo.tf
rm meso.mhb.hypo.tf
rm ecto.mhb.hypo.tf
for i in `ls /home/shg047/db/hg19/encode/encode.*.hg19.bed`
do
bedtools window -w 100 -a endo.mhb.hypo.bed -b $i >> endo.mhb.hypo.tf
bedtools window -w 100 -a meso.mhb.hypo.bed -b $i >> meso.mhb.hypo.tf
bedtools window -w 100 -a ecto.mhb.hypo.bed -b $i >> ecto.mhb.hypo.tf
done
cat endo.mhb.hypo.tf | awk '{print $7}' | sort -u > endo.mhb.hypo.uni.tf
cat meso.mhb.hypo.tf | awk '{print $7}' | sort -u > meso.mhb.hypo.uni.tf
cat ecto.mhb.hypo.tf | awk '{print $7}' | sort -u > ecto.mhb.hypo.uni.tf


for i in `ls /home/shg047/db/hg19/encode/encode.*.hg19.bed`
do
bedtools intersect -a $i -b ecto.mhb.hypo.bed 
done


/home/shg047/oasis/monod/bin/bam2hapinfo.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.chr6.mld_blocks_r2-0.5.bed  /home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/Colon_primary_tumor.chr7.sorted.clipped.bam > Colon_primary_tumor.chr7.hapInfo.txt

/home/oasis/db/
cd salk
perl ../bam2hapinfoPBS.pl ../Ecker_Tissue_WGBS_sampleInfo_WGBS-pooled-mld-blocks.txt
cd ../tumor_wgbs/
perl ../bam2hapinfoPBS.pl ../tumor_WGBS_sample_info_WGBS-pooled-mld-blocks.txt
cd ../wbc/
perl ../bam2hapinfoPBS.pl ../TSCC_whole_blood_WGBS_sampleInfo_WGBS-pooled-mld-blocks.txt
cd ../n37/
perl ../bam2hapinfoPBS.pl ../N37_10_tissue_sampleInfo_WGBS-pooled-mld-blocks.txt
cd ../hesc/
perl ../bam2hapinfoPBS.pl ../H1ESC_sampleInfo_WGBS-pooled-mld-blocks.txt

HWI-D00506:1:1102:13764:14978#0/3:F     16      chr5    177653119	14      84M     *	0	0	CAAACTAAAATACAATAACGCGATCTCAACTCACTACAACCTCTACCTCTCAAACTCAAACAATTCTCCTACTTCAACCTCCCA    bbbeeeeeggfggiiihhhiifcgfiiihihiiihhi
HWI-D00506:1:1106:21115:64108#0/1_HWI-D00506:1:1106:21115:64108#0/3_121:R	16      chr5    177653119	14      121M    *	0	0	CAAACTAAAATACAATAACTCTATCTCAACTCACTACAACCTCTACCTCTCAAACTCAAACAATTCTCCTACTTCAACCTCCCAA
HWI-D00506:1:1102:13764:14978#0/1:R     16      chr5    177653121	14      84M     *	0	0	AACTAAAATACAATAACGCGATCTCAACTCACTACAACCTCTACCTCTCAAACTCAAACAATTCTCCTACTTCAACCTCCCAAA    bbbeeeeefcggghhgfghiihhiihfhihhigfhii

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/ 
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/annotation/hg19.c2t.fa -p chr5:177653138
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/db/hg19/hg19.fa -p chr5:177653138

# I found MHL at chr5:177653138-177653229 in 6-P-1 is missing, however, there are 3 reads aligned to this region in the bam file. 
# I need check what's wrong with it? Anything should be de-bug for bam2hapinfo.pl? 
mkdir 
samtools view -b /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam chr5:177653138-177653229 -o 6-P-1.test.bam 
perl ../../bin/bam2hapinfo.pl target.bed 6-P-1.test.bam
grep 177653138 6-P-1*
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/db/hg19/hg19.fa -p chr5:177653138
samtools view 6-P-1.test.bam 
samtools view -q 6-P-1.test.bam 


Illumicode  MouseWG-6 v2.0



 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o bam2mhb.log
 #PBS -e bam2mhb.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/hapinfo
 perl ../bin/hapinfo2mhl.pl ./




cd /home/shg047/oasis/monod/haplo/n37
perl hapMergeByChrosome.pl
mv N37-CRBL.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Colon.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-FL.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Heart.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Liver.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Lung.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Pancreas.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-SI.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-SM.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Stomach.hapInfo.txt /home/shg047/oasis/monod/hapinfo/

cd /home/shg047/oasis/monod/haplo/salk
perl ../n37/hapMergeByChrosome.pl ./
mv STL001BL-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001FT-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001GA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001LG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001LV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001PO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001RV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001SB-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001SG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001SX-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001TH-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002AD-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002AO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002EG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002FT-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002GA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002LG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002OV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002PA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002PO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002SB-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002SX-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003AD-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003AO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003EG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003FT-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003GA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003LV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003PA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003PO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003RA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003RV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003SB-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003SG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003SX-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL011LI-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/

cd /home/shg047/oasis/monod/haplo/wbc
perl ../n37/hapMergeByChrosome.pl ./
mv centenarian.hapInfo.txt   /home/shg047/oasis/monod/hapinfo/
mv middle-age.hapInfo.txt   /home/shg047/oasis/monod/hapinfo/
mv new-born.hapInfo.txt   /home/shg047/oasis/monod/hapinfo/

cd /home/shg047/oasis/monod/haplo/hesc
perl ../n37/hapMergeByChrosome.pl ./
mv methylC-seq_h1+bmp4_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1+bmp4_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-msc_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-msc_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-npc_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-npc_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_mesendoderm_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_mesendoderm_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/





cd /home/shg047/monod/methyblock/HM450K/mhb450bed.bed 
bedtools intersect -wa -u -a mhb450bed.bed -b ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  # 1258
bedtools shuffle -i /home/shg047/monod/methyblock/HM450K/mhb450bed.bed -incl /home/shg047/oasis/db/hg19/CRGmapability.hg19.exclude.bed -g ~/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l  # random 48, observe: 
bedtools intersect -wa -u -a /home/shg047/monod/methyblock/HM450K/mhb450bed.bed -b /home/shg047/oasis/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l  # random 48, observe: 1258


cd /home/shg047/monod/methyblock/encode_rrbs_bed
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l   # 8920
bedtools shuffle -i /home/shg047/monod/methyblock/encode_rrbs_bed/encode.rrbs.mhb.txt -incl /home/shg047/oasis/db/hg19/CRGmapability.hg19.exclude.bed -g ~/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l  # random 337, observe: 




bedtools shuffle -i /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -incl /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -g /home/shg047/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l





/home/shg047/monod/rrbs_kun/*hapInfo.txt

load("Encode.RRBS.MHB.RData")
cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}
rrbsmhb<-cor2bed(rownames(subset1))
write.table(rrbsmhb,file="encode.rrbs.mhb.txt",col.names=F,row.names=F,quote=F,sep="\t")
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
bedtools intersect -wa -u -b encode.rrbs.mhb.txt -a ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l

../../mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_targets_LOD_CC.txt &
../../mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > RRBS_targets_LOD_LC.txt &
../../mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_targets_LOD_PC.txt &

cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/expanded_analysis_31Aug2014

perl /home/kunzhang/CpgMIP/MONOD/Data/mixMethHapAnalysisStringent_27Aug14.pl NC-P-ALL.hapInfo.txt 6-T-1.hapInfo.txt 6-P-1.hapInfo.txt  /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_RRBS_targets.BED.txt > /home/shg047/monod/mixHap/zhang/6-P-1_tumor_hap.txt

grep 105309703-105309944 NC-P-ALL.hapInfo.txt 
grep 105309703-105309944 6-T-ALL.hapInfo.txt 
grep 105309703-105309944 6-P-ALL.hapInfo.txt




mixMethHapAnalysisStringent_27Aug14

8920/
diff mixMethHapAnalysisStringent_27Aug14.pl mixMethHapAnalysisStringent_02Nov14.pl



my $hap_file_NCP = $ARGV[0];
my $hap_file_Tumor = $ARGV[1];
my $hap_file_PP = $ARGV[2];
my $mC_bed_file_NCP = $ARGV[3];
my $informative_probe_file = $ARGV[4];


my $hap_file_NCP = "NC-P-ALL.hapInfo.txt";
my $hap_file_Tumor = "6-T-ALL.hapInfo.txt";
my $hap_file_PP = "6-P-ALL.hapInfo.txt";
my $mC_bed_file_NCP = "/home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt";
my $informative_probe_file = "";





cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/expanded_analysis_31Aug2014

perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/kunzhang/CpgMIP/MONOD/DataRRBS_targets_LOD_CC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_LC.txt &
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_PC.txt &


# 2016-02-26
wget https://cran.r-project.org/src/contrib/mhsmm_0.4.14.tar.gz
install.packages("mhsmm_0.4.14.tar.gz")
source("http://bioconductor.org/biocLite.R")
biocLite("MethylSeekR")
Hidden Markov Model to identify PMDs, UMRs and LMRs

/home/shg047/oasis/monod/bin/LDR2extent.pl chr5:112039140-112047142 APC.chr5.sqr > APC.full.region.sqr
head -n 100 /home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt > head.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  APC.ext chr5:112043054-112197528 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl /home/shg047/bin/LDR2extent.pl chr5:112043054-112197528 APC.chr5.rsq /home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt > APC.chr5.rsq.ext
my $cpg_position_file="/home/shg047/oasis/db/hg19/meth/bismark/HsGenome19.CpG.positions.txt";
# APC regions
cd /home/shg047/oasis/monod/haplo/
cat /home/shg047/oasis/monod/haplo/n37/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/wbc/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/hesc/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/tumor_wgbs/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/salk/*hapInfo.txt >> hapinfo.txt

cd /home/shg047/oasis/monod/haplo/
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  APC chr5:111986595-112228720 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  ZNF154 chr19:58206356-58222715 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  SHOX2 chr3:157808200-157829413 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  DIRAS3 chr1:68511645-68516481 < /home/shg047/oasis/monod/haplo/hapinfo.txt


awk 'NR==9467146 {print}' hapinfo.txt
awk 'NR==8693459 {print}' hapinfo.txt
awk 'NR==8693437 {print}' hapinfo.txt

chr5:111986595-112228720
.
# 2016-02-24
SRR1035775
SRR1035774
SRR1035786
SRR1035787
SRR1035788
SRR1035818
SRR1035778
SRR1035893
SRR1035815
SRR1035776
SRR1035893
SRR1035814
SRR1035774
SRR1035777
SRR1035813
SRR1035777
SRR1035814
SRR1035773
SRR1035788
SRR1035787
SRR1035775
SRR1035773
SRR1035777
SRR1035786
SRR1035786
SRR1035776
SRR1035778
SRR1035822
SRR1035815

􀀛􀀚􀀘􀀃􀀫􀁒􀁗􀁈􀁏􀀃􀀦􀁌􀁕􀁆􀁏􀁈􀀃􀀶􀁒􀁘􀁗􀁋􀀏􀀃􀀰􀁌􀁖􀁖􀁌􀁒􀁑􀀃􀀹􀁄􀁏􀁏􀁈􀁜􀀏􀀃􀀶􀁄􀁑􀀃􀀧􀁌􀁈􀁊􀁒


bismark_methylation_extractor -p --no_overlap --ignore_r2 2 --comprehensive --gzip --report --multicore 4 --bedGraph --ample_memory --cytosine_report --CX --genome_folder ~/genomes/Heinz_pennellii_organelles_F1_genome/ --split_by_chromosome -o methyl_extraction/ P1.deduplicated.bam

# cytosineReport to MethylKit
awk  '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1,$2,$3,$4/($4+$5),$4+$5;}' cytosineReport.txt > outputForMethylKit

fastq-dump --split-files --skip- --gzip $id 

for i in {4391177..4391187}; do qdel $i.tscc-mgr.local; done

# SRR390728 is pair-end
fastq-dump -X 20 -Z SRR390728
fastq-dump --split-files  -X 20 -Z SRR390728 
fastq-dump --split-3  -X 20 -Z SRR390728 

fastq-dump -X 20 SRR390728
fastq-dump --split-files  -X 20 SRR390728 
fastq-dump --split-3  -X 20 SRR390728 
fastq-dump -X 5 -Z SRR390728

# SRR2542443 is single-end
fastq-dump --split-files -X 20 SRR2542443 

cd /home/shg047/oasis/Holger2016
screen -S SRA112056.Download 
wget -r -c -l 2 -nH --cut-dirs=4 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056/
qwget -r -c -l 2  -nH  ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056/
wget -r -l 2 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056

# 2016-02-21
scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt ./
cd /home/kunzhang/CpgMIP/MONOD/Data/
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_CC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_LC.txt &
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_PC.txt &

/home/kunzhang/CpgMIP/MONOD/Data/mixMethHapAnalysis_27Oct14.pl
/home/kunzhang/CpgMIP/MONOD/Data/mixMethHapAnalysisStringent_02Nov14.pl

cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/expanded_analysis

scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/MONOD/Data/*pl ./
scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/*pl ./

cd /home/kunzhang/CpgMIP


mixMethHapAnalysisStringent_27Aug14.pl

cd /home/shg047/monod/hapinfo
cat /home/shg047/monod/rrbs_kun/6-P-*.hapInfo.txt > ../mixHap/6-P-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/7-P-*.hapInfo.txt > ../mixHap/7-P-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/PC-P-*.hapInfo.txt > ../mixHap/PC-P-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/NC-P-*.hapInfo.txt > ../mixHap/NC-P-ALL.hapInfo.txt

cat /home/shg047/monod/rrbs_kun/6-T-*.hapInfo.txt > ../mixHap/6-T-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/7-T-*.hapInfo.txt > ../mixHap/7-T-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/PC-T-*.hapInfo.txt > ../mixHap/PC-T-ALL.hapInfo.txt

cd /home/shg047/monod/mixHap
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_target_LOD_CC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > RRBS_target_LOD_LC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_target_LOD_PC.txt 
./find_NC-P_high-LOD_targets.pl > NC-P_high-LOD_targets.txt
  
  
find /home/kunzhang/CpgMIP/MONOD/Data/ -name *pl | grep find_NC-P_high
cp /home/kunzhang/CpgMIP/MONOD/Data/*.sh ./


# creat bisulfite treated human reference: hg19
cp hg19.fa hg19.c2t.fa
perl -p -i -e 's/CG/M/ig' hg19.c2t.fa
perl -p -i -e 's/C/T/ig' hg19.c2t.fa
perl -p -i -e 's/M/CG/ig' hg19.c2t.fa

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/ 
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:100027865

chr10:101089382-101089519
Sept 2014: these batch of libraries were generated with the dRRBS protocol (MSP I & Taq I digestion)
Bam folder: /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles
Mapple_bin_hapInfo folder:
Mld_block_hapInfo folder: /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo


cd /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo
head NC-P-23.mld_blocks_r2-0.5.hapInfo.txt

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles

# index: http://genome-tech.ucsd.edu/LabNotes/index.php/Dinh/Dinh_2014/NOTES/2014-9-22 
samtools tview Indx16_S12.sorted.clipped.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:102440557-102440826
samtools view Indx16_S12.sorted.clipped.bam  chr10:102443911-102443400


chr10:101281221-101281239
chr10:101281274-101281291
chr10:101281879-101281896 
# 2016-02-19
head /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed 

/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles

cpg_list="/media/Ext12T/DD_Ext12T/BisRef/bisHg19_plusLambda_BWA/hg19_lambda.cpg.positions.txt"
samtools mpileup -BA -f hg19.fa $bam_file > $pileup_file
/home/dinh/scripts/BisReadMapper/src/extractMethyl.pl $cpg_list 33 < $pileup_file > $methylFreq_file

awk '{print $3-$2}' /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed

awk '{print $3-$2}' Haib.merge_RD10_80up.mld_blocks_r2-0.1.bed | sort -u
awk '{print $3-$2}' /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  |sort -n -u
awk '{print $3-$2}' /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed  | sort -n -u

scp shg047@genome-miner.ucsd.edu:/home/shg047/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed ./
scp shg047@genome-miner.ucsd.edu:/home/shg047/monod/methyblock/*.bed ./

bedtools intersect -wa -u -a  /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  | wc -l

# test 1
wc -l /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed
bedtools intersect -wa -u -a /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l 
wc -l /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed
bedtools shuffle -i /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -incl /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -g /home/shg047/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l

# test 2
wc -l /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed
bedtools intersect -wa -u -a /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l 
wc -l /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed
bedtools shuffle -i /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -excl -g /home/shg047/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l


# Cancer specific methylation haplotype
http://genome-tech.ucsd.edu/LabNotes/index.php/Cancer_specific_methylation_haplotype
http://genome-tech.ucsd.edu/LabNotes/index.php/Shicheng:Calendar/NOTES/2016-2-1
http://genome-tech.ucsd.edu/LabNotes/index.php/Kun:LabNotes/MONOD/2014-10-8_WGBS_RRBS
http://genome-tech.ucsd.edu/LabNotes/index.php/Kun:LabNotes/MONOD/2014-10-8_WGBS_RRBS


less /home/shg047/oasis/Holger2016/sra/SRR1035896.trim.err
 scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/MONOD/Data/methHapClassfier.pl ./
 methHapClassfier.pl  N37-Cerebellum.chr1.hapInfo.txt > N37-Cerebellum.methHapCounts.txt
 cat PC-T-*.methHapCounts.txt | ../scripts/get_tumor_specific_HMH_regions.pl NC-P-ALL_RRBS-dRRBS.methHapCounts.txt > PC-T_NC-plasma_specific_HMH_regions.bed
 cat 7-T-*.methHapCounts.txt | ../scripts/get_tumor_specific_HMH_regions.pl NC-P-ALL_RRBS-dRRBS.methHapCounts.txt > LC-T_NC-plasma_specific_HMH_regions.bed
 cat 6-T-*.methHapCounts.txt CTT-frozen-100ng.methHapCounts.txt  | ../scripts/get_tumor_specific_HMH_regions.pl NC-P-ALL_RRBS-dRRBS.methHapCounts.txt > CRC-T_NC-plasma_specific_HMH_regions.bed

Question 1: Cannot locate ... in @INC - Perl Maven
Question 2: How to install the module
Question 3: Where I have installed my module 
Question 4: how to load module

Answer 1:
module path is not in the @INC. You need add the path to @INC

Answer 2:  
cpan
install Sort::Array

Answer 3:  
perldoc -l XML::Simple
perldoc -l Sort::Array

Answer 4:
export PERL5LIB=$PERL5LIB:/home/shg047/perl5/perlbrew/perls/perl-5.22.0/lib/site_perl/5.22.0/Sort/
export PERLLIB=$PERLLIB:/home/shg047/perl5/perlbrew/perls/perl-5.22.0/lib/site_perl/5.22.0/Sort/
source ~/.bashrc 


qsub SRR1035860.job
qsub SRR1035849.job
qsub SRR1035856.job
qsub SRR1035860.job
qsub SRR1035863.job
qsub SRR1035869.job
qsub SRR1035879.job
qsub SRR1035888.job
qsub SRR1035895.job

 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -o bam2mhb.log
 #PBS -e bam2mhb.err
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/Haib/sortBam
 ../mergedBam2hapInfo.pl ./haib.RD10_80up.genomecov.bed haib.merge.sort.bam > Haib.merge.RD10_80up.hapinfo.txt  # get hapinfo
 ../hapInfo2mld_block.pl ./Haib.merge.RD10_80up.hapinfo.txt 0.5 >  Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed   # get methylation block


 #!/bin/csh
 #PBS -n CTR97_trimmed.fq.gz_bismark_bt2.sort.bam
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=16
 #PBS -l walltime=72:00:00
 #PBS -o CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log
 #PBS -e CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/Haib/sortBam
 # samtools cat -h header.sam -o haib.merge.bam *sort.bam
 samtools sort -@ 16 haib.encode.merge.bam -o haib.merge.sort.bam
 samtools index haib.merge.sort.bam
 bedtools genomecov -bg -split -ibam haib.merge.sort.bam >   haib.merge.bam.pool.bed
 awk '$4>9 { print $1"\t"$2"\t"$3}'  haib.merge.bam.pool.bed | bedtools merge -d 10 -i - > haib.RD10.genomecov.bed
 awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' haib.RD10.genomecov.bed > haib.RD10_80up.genomecov.bed

 for i in {1..22,"X","Y"} 
 do 
 grep chr1 haib.RD10_80up.genomecov.bed |awk '{sum+=($3-$2)} END { print "chr1 Sum = ", sum, " Average = ",sum/NR, "N = ", NR}' >> N37_WGBS_tumor_seqCap_RRBS_tumor_NC_RD10_80up.mld_blocks.summary.txt
 done
 
  for i in `seq 1 12`;
  do
  grep chr$1 haib.RD10_80up.genomecov.bed |awk '{sum+=($3-$2)} END { print "chr1 Sum = ", sum, " Average = ",sum/NR, "N = ", NR}' >> Haib_RD10_80up.mld_blocks.summary.txt
  done 
 
 
 mkdir /home/shg047/oasis/Haib/mhb
 mv 
 genome.cov.sh					 

mv haib.encode.merge.bam ../mhb			 
mv CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err  ../mhb   
mv haib.merge.sort.bam	    ../mhb		     
mv haib.merge.sort.bam.bai      ../mhb		      
mv haib.merge.bam.pool.bed	 ../mhb		    
mv haib.RD10.genomecov.bed	    ../mhb		 
mv haib.RD10_80up.genomecov.bed	 ../mhb	      
mv CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log     ../mhb

 
# 2016-2-16
Dinh methylation collection page @ genome-miner
/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge
DNA Methylation by Reduced Representation Bisulfite Seq from ENCODE/HudsonAlpha 

cd /home/kunzhang/CpgMIP/MONOD/Data
wget -r -l 2 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056


# 2016-2-13
(head -n 1 $file && tail -n +2 $file | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $sample.sort.bedGraph
bedGraphToBigWig NC-P-2.bedGraph_CpG.sort.bedGraph hg19.chrom.sizes  NC-P-2.bw
fetchChromSizes hg19 > hg19.chrom.sizes
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes

echo {0..9} | xargs -n 2

shicheng@meangenemachine2
ifconfig: 132.239.189.199
sudo passwd: fudan1108
sudo apt-get update
sudo apt-get install vim
sudo apt-get install vsftpd
sudo vim /etc/vsftpd.conf
# Listen=YES
anonymous_enables=NO
local_enable=YES
write_enable=YES
xferlog_file=/var/log/vsftpd.log
ftpd_banner= Welcome to my new FTP Serve at UCSD
save and exit
sudo cat /var/log/vsftpd.log
sudo service vsftpd restart
sudo adduser ucsd002  passwd: ucsd002

 
sudo telnet localhost 21
ps -aux | grep vsftpd
sudo service vsftpd restart
sudo netstat -ntaulp | grep vsftpd
sudo adduser genemean001  passwd: genemean001
sudo adduser ucsd002  passwd: ucsd002


 
# 2016-2-12
cd ~/bin/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigEncode
chmod 700 * 
cd -

wigCorrelate

#!/bin/csh
#PBS -n CTR97_trimmed.fq.gz_bismark_bt2.sort.bam
#PBS -q glean
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log
#PBS -e CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group

cd /home/shg047/oasis/biomark
PileOMeth extract  -q 1 -p 1 --minDepth 1 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o 7-P-3.PileOMeth.bedGraph
samtools depth -b ~/oasis/db/hg19/meth/bismark/HsGenome19.CpG.positions.txt /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam > 7-P-3.coverage
samtools tview -p chr1:10923 7-P-3.sorted.clipped.bam ~/oasis/db/hg19/meth/bismark/hg19.fa

5697312 7-P-3.coverage

37076

samtools tview -p chr1:10923 7-P-3.sorted.clipped.bam ~/oasis/db/hg19/meth/bismark/hg19.fa


[shg047@tscc-login1 biomark]$ wc -l 7-P-3.PileOMeth.bedGraph_CpG.bedGraph
37076 7-P-3.PileOMeth.bedGraph_CpG.bedGraph

/home/shg047/oasis/monod/bam
less /home/shg047/oasis/monod/bam/*coverage

chr1    10854   15
chr1    10857   16
chr1    10860   23
chr1    10866   23
chr1    10884   24
chr1    10886   24
chr1    10902   24
chr1    10907   24
chr1    10913   24
chr1    10923   23
chr1    10928   23
chr1    10930   23
chr1    10933   23
chr1    10936   1
chr1:10884
cd /home/shg047/oasis/biomark
less 7-P-3.PileOMeth.bedGraph_CpG.bedGraph
chr1    10857   10858   100     11      0
chr1    10860   10861   81      9	2
chr1    10866   10867   100     11      0
chr1    10884   10885   54      6	5
chr1    10886   10887   36      4	7
chr1    10902   10903   100     11      0
chr1    10907   10908   100     11      0
chr1    10913   10914   90      10      1
chr1    10923   10924   90      9	1
chr1    10928   10929   90      9	1
chr1    10930   10931   80      8	2

cd /home/shg047/oasis/monod/bam

samtools view -h -b 7-P-3.sorted.clipped.bam chr1:10854-10933 -o hg19.bsseq.chr1.10854-10933.bam
samtools tview -p chr1:10854 hg19.bsseq.chr1.10854-10933.bam ~/oasis/db/hg19/meth/bismark/hg19.fa
samtools mpileup -A -Q 1 --reference ~/oasis/db/hg19/meth/bismark/hg19.fa hg19.bsseq.chr1.10854-10933.bam
PileOMeth extract  -q 0 -p 1 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa hg19.bsseq.chr1.10854-10933.bam  -o hg19.bsseq.chr1.10854-10933.bam.PileOMeth.bedGraph

PileOMeth extract --mergeContext --keepSingleton --keepDiscordant /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa hg19.bsseq.chr1.10854-10933.bam  -o hg19.bsseq.chr1.10854-10933.bam.PileOMeth.bedGraph



# 2016-2-10
Getopt::Std
  getopt  	
     getopt ('lw'); $opt_l and $opt_w to accept the value from terminal 
	 getopts ('abl:w:'); 

Getopt::Long
   GetOptions ('a|all' => \$all, 'l|length=i' => \$length,'w|width=i' => \$width);
     * for the option words, a double dash is required: ‘--length 24’ is acceptible
     * reference of $varabile should be taken as destination
     * You do not need to specified the option destination. If no destination is specified, GetOptions will define variables $opt_xxx where xxx is the name of the option, just like getopt and getopts. GetOptions will also accept a reference to a hash as its first argument and deliver the option values there, again just like getopt and getopts.
   GetOptions ('foo=i' => \@values);
     * Calling this program with arguments ‘-foo 1 -foo 2 -foo 3’ will result in @values having the value (1,2,3) provided it was initially empty.
     * Also, the option destination can be a reference to a hash. In this case, option values can have the form ‘key=value’. The value will be stored in the hash with the given key.


# 2016-2-10
PileOMeth extract  -q 10 -p 5 --minDepth 5 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o 7-P-3.bedGraph 
PileOMeth extract  -q 10 -p 5 --minDepth 5 -l biomark.list.bed /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o  a

perl pipeline.pl biomark.list.bed /oasis/tscc/scratch/shg047/monod/bam/


bowtie2 -a -x ~/db/aligndb/hg19/bismark/Bisulfite_Genome/CT_conversion/BS_CT  ZNF154.fastq
bowtie2 -a -x ~/db/aligndb/hg19/bismark/Bisulfite_Genome/GA_conversion/BS_GA  ZNF154.fastq



BS_CT.rev.1.bt2
 perl -p -i -e 's/--paired-end/--single-end/g' *job

 
/oasis/tscc/scratch/ddiep/Working/Rerun_rrbs/BAMfiles/*


for i in `ls *methylFreq`
do
perl methylFreq2wig.pl $i 5 < $i > $i.wig 
gzip $i.wig
done


# 2016-2-5
@ZNF154-FP
GGTTTTTATTTTAGGTTTGA
+
aaaaaaaaaaaaaaaaaaaa
@ZNF154-RP
AAATCTATAAAAACTACATTACCTAAAATACTCTA
+
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
cd /home/shg047/work/biomarker/ZNF154
bismark --bowtie2 --non_directional --ambiguous --phred64-quals --fastq -L 10 -N 1 -s 0 --multicore 1 /home/shg047/db/aligndb/hg19/bismark ZNF154.
# 2016-2-5
bwa mem -O 0 -R '@RG    ID:ON1_S1	SM:ON1  PL:ILLUMINA     LB:On_S1' /home/zhl002/ZL_LTS33/RNA_seq/keiOTS_01262016/CAGOn-Cas_S1.fa 1-CAGOn-Cas_S1_L001_R1_001.fastq > CAGOnS1.fastq.sam
/home/kunzhang/softwares/samtools-latest/samtools view -b -t /home/zhl002/ZL_LTS33/RNA_seq/keiOTS_01262016/CAGOn-Cas_S1.fa.fai -h CAGOnS1.fastq.sam > CAGOnS1.fastq.bam
java -Xmx4g -jar /home/shg047/software/picard-tools-1.113/SortSam.jar TMP_DIR=./tmp/ INPUT=CAGOnS1.fastq.bam OUTPUT=CAGOnS1.sorted.bam QUIET=True SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
/home/kunzhang/softwares/samtools-latest/samtools index CAGOnS1.sorted.bam
java -Djava.io.tmpdir=./tmp -Xmx4g -jar /home/kunzhang/softwares/GenomeAnalysisTK-3.3/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/zhl002/ZL_LTS33/RNA_seq/keiOTS_01262016/CAGOn-Cas_S1.fa -I CAGOnS1.sorted.bam -glm INDEL -stand_call_
CAGOnS1.sh (END)
find -type f -iname '*.gz' -exec gunzip -t {}
/bsmap-2.74/bsmap -u -s 12 -v 0.04 -p 4 -a downsampled.SRR1232303_1.fastq -d hg19.fa -o downsampled.SRR1232303_1.fastq.sam
bsmap -u -s 12 -v 0.04 -p 6 -a ../fastq_trim/SRR1232304_1_trimmed.fq.gz -b ../fastq_trim/SRR1232304_2_trimmed.fq.gz -d  /home/shg047/oasis/Xliu2014/fastq -o SRR1232304_2_bam 
# 2016-2-4
cd /home/shg047/monod/hapinfo 
perl cgLD_Analysis_haploInfo_permuteRsq_v2_makeLD-plot.pl result.rsq chr10:101089204-101089305 < RRBS-7P30.hapInfo.txt
perl cgLD_Analysis_haploInfo_permuteRsq_v2_makeLD-plot.pl result.rsq chr10:101294771-101294802  < RRBS-7P30.hapInfo.txt
perl cgLD_Analysis_haploInfo_permuteRsq_v2_makeLD-plot.pl result.rsq chr10:103113834-103114495  < RRBS-7P30.hapInfo.txt
chr10:100995942-100996028

SRR1232310_1_val_1.fq.gz.temp.2_bismark_bt2_PE_report.txt
SRR1232313_1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
# 2016-1-27
Pregnancy.10.read1.fq.gz.job:samtools sort ../bam/Pregnancy.10.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.10.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.11.read1.fq.gz.job:samtools sort ../bam/Pregnancy.11.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.11.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.1.read1.fq.gz.job:samtools sort ../bam/Pregnancy.1.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.2.read1.fq.gz.job:samtools sort ../bam/Pregnancy.2.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.2.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.3.read1.fq.gz.job:samtools sort ../bam/Pregnancy.3.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.3.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.4.read1.fq.gz.job:samtools sort ../bam/Pregnancy.4.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.4.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.5.read1.fq.gz.job:samtools sort ../bam/Pregnancy.5.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.5.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.6.read1.fq.gz.job:samtools sort ../bam/Pregnancy.6.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.6.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
T21.3.read1.fq.gz.job:samtools sort ../bam/T21.3.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/T21.3.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
T21.4.read1.fq.gz.job:samtools sort ../bam/T21.4.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/T21.4.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
# 2016-1-27

/home/kunzhang/CpgMIP/MONOD/Public_data/find_LMS.pl > whole_blood_LMS.txt
./extract_clusters.pl whole_blood_LMS.txt > whole_blood_LMS_clusters.txt

./report_tumor_HMH_regions.pl UCSD-004-05_sampleInfo.txt  | sort -k 3nr > UCSD-004-05_HMH_regions.txt
./report_tumor_HMH_regionsreport_tumor_HMH_regions.pl UCSD-004-07_sampleInfo.txt  | sort -k 3nr > UCSD-004-07_HMH_regions.txt
  
cat s_1_1_ILMN_Indx04.hapInfo.txt s_1_1_ILMN_Indx05.hapInfo.txt s_1_2_ILMN_Indx04.hapInfo.txt s_1_2_ILMN_Indx05.hapInfo.txt > PC-T-2.pooled.hapInfo.txt     
cat s_2_1_ILMN_Indx02.hapInfo.txt s_2_2_ILMN_Indx02.hapInfo.txt > PC-P-2.pooled.hapInfo.txt      
cat s_2_1_ILMN_Indx14.hapInfo.txt s_2_2_ILMN_Indx14.hapInfo.txt s_2_1_ILMN_Indx15.hapInfo.txt s_2_2_ILMN_Indx15.hapInfo.txt s_2_1_ILMN_Indx16.hapInfo.txt s_2_2_ILMN_Indx16.hapInfo.txt s_2_1_ILMN_Indx27.hapInfo.txt s_2_2_ILMN_Indx27.hapInfo.txt NP-RRBS-NC-P-1ng-p2-Jun20-A013.hapInfo.txt > NC-plasma.pooled.hapInfo.txt
 
../prepare_plotting_files.pl NC-plasma.pooled.hapInfo.txt PC-T-2.pooled.hapInfo.txt PC-P-2.pooled.hapInfo.txt chr6:152128536-152129155
../prepare_plotting_files.pl NC-plasma.pooled.hapInfo.txt PC-T-2.pooled.hapInfo.txt PC-P-2.pooled.hapInfo.txt chr5:176543913-176544076

# 2016-1-20

Cython: http://cython.org/
tar xzvf Cython-0.23.4.tar.gz
cd ./Cython-0.23.4/
sudo python setup.py install

wheel: https://pypi.python.org/pypi/wheel
tar xzvf wheel-0.26.0.tar.gz
cd ./wheel-0.26.0/
sudo python setup.py install

setuptools: https://packaging.python.org/en/latest/projects/#setuptools
tar xzvf setuptools-19.4.zip
cd ./setuptools-19.4
sudo python setup.py install

NumPy: http://www.scipy.org/scipylib/download.html
git clone git://github.com/numpy/numpy.git numpy
cd ./numpy/
sudo python setup.py install

MethylPurify: https://pypi.python.org/pypi/MethylPurify
tar xzvf MethylPurify-2.0-20141116.tar.gz
cd MethylPurify-2.0-20141116/
sudo python setup.py install --user
MethylPurify 


/bsmap-2.74/bsmap -u -s 12 -v 0.04 -p 4 -a downsampled.SRR1232303_1.fastq -d hg19.fa -o downsampled.SRR1232303_1.fastq.sam
python2.7 /home/ddiep/Downloads/MethylPurify-2.0-20140819/methylpurify/bin/MethylPurify -f downsampled.SRR1232303_1.fastq.bam 
      -g hg19.fa -i /home/ddiep/Downloads/MethylPurify-2.0/methylpurify/db/CGI_hg19_slop1000.bed -c 10 -s 50 -b 300
 
 
compare the difference between new and old perl script of GATK
old：/home/ajgore/AG_Ext12T/GATK_01022012/variantCallerBwaGATK-latest
new：/home/kunzhang/bin/variantCallerBwaGATK_05012015.pl


# 2016-1-19
bismark.zero.cov

ENCFF000LUN_trimmed

chr8:22249850-22249898

samtools tview -p chr8:22249850 /home/shg047/oasis/Haib/bam/ENCFF000LUN_trimmed.fq.gz_bismark_bt2.sort.bam /home/shg047/db/hg19.fa

qsub Lymphoma.run1.*job
qsub Pregnancy.14.*job
qsub T21.2.read1.*job

-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 07:17 Lymphoma.run1.read1_val_1.fq.gz.temp.4_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  453 Jan 19 07:19 Lymphoma.run1.read1_val_1.fq.gz.temp.5_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group 3.7K Jan 19 07:19 Lymphoma.run1.read1_val_1.fq.gz.temp.4_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  453 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.5_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  23K Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.6_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.1_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.3_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.2_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 09:20 Pregnancy.14.read1_val_1.fq.gz.temp.5_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 09:20 Pregnancy.14.read1_val_1.fq.gz.temp.3_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 09:21 Pregnancy.14.read1_val_1.fq.gz.temp.4_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  451 Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group 146K Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.3_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.6_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.4_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  25K Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.5_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:18 T21.2.read1_val_1.fq.gz.temp.4_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:18 T21.2.read1_val_1.fq.gz.temp.3_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:18 T21.2.read1_val_1.fq.gz.temp.2_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.1_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  437 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.5_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.5_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  22K Jan 19 11:21 T21.2.read1_val_1.fq.gz.temp.2_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.2G Jan 19 17:39 T21.2.read1_val_1.fq.gz.temp.1_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.2G Jan 19 17:40 T21.2.read1_val_1.fq.gz.temp.4_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 2.3G Jan 19 17:40 T21.2.read1_val_1.fq.gz.temp.3_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 2.3G Jan 19 17:40 T21.2.read1_val_1.fq.gz.temp.6_bismark_bt2_pe.bam


# 2015-12-27
qsub Lymphoma.run1.read1.fq.gz.job
qsub Lymphoma.run2.read1.fq.gz.job
qsub Lymphoma.run3.read1.fq.gz.job
qsub Pregnancy.14.read1.fq.gz.job
qsub Pregnancy.9.run2.read1.fq.gz.job
qsub T21.2.read1.fq.gz.job

rm Pregnancy.13.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.12.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.15.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.7.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.8.run2.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.9.run2.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.8.run1.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.9.run1.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.14.read1_trimmed.fq.gz_bismark_bt2_*
rm T21.2.read1_trimmed.fq.gz_bismark_bt2_*
rm T21.5.read1_trimmed.fq.gz_bismark_bt2_*
rm T21.1.read1_trimmed.fq.gz_bismark_bt2_*


ls Lymphoma.run3* | grep -v val 

ls *pe 

qsub HL_05.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_06.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_07.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_08.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_09.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_14.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_15.sorted.clipped.fastq.gz.bismark.pbs

ls *LTP* | grep -v val

rm LTP2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
rm LTP2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
rm LTP3.read1_trimmed.fq.gz_bismark_bt2_pe.bam
rm LTP3.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
rm LTP4.read1_trimmed.fq.gz_bismark_bt2_pe.bam
rm LTP4.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt

ls -larth *Pregnancy* | grep -v val

-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 00:30 Pregnancy.5.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  21G Jan 15 00:30 Pregnancy.5.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 02:03 Pregnancy.4.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  23G Jan 15 02:03 Pregnancy.4.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  17G Jan 15 02:45 Pregnancy.10.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 02:45 Pregnancy.10.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 04:04 Pregnancy.6.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  22G Jan 15 04:04 Pregnancy.6.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 04:36 Pregnancy.2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  26G Jan 15 04:36 Pregnancy.2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 07:23 Pregnancy.3.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  30G Jan 15 07:23 Pregnancy.3.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 09:03 Pregnancy.11.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  24G Jan 15 09:03 Pregnancy.11.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 11:15 Pregnancy.1.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  31G Jan 15 11:15 Pregnancy.1.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:42 Pregnancy.13.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  707 Jan 15 13:42 Pregnancy.13.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:42 Pregnancy.12.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  707 Jan 15 13:42 Pregnancy.12.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:44 Pregnancy.15.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  707 Jan 15 13:44 Pregnancy.15.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:44 Pregnancy.7.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  706 Jan 15 13:44 Pregnancy.7.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.8.run2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.8.run2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.9.run2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.9.run2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.8.run1.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.8.run1.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.9.run1.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.9.run1.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:45 Pregnancy.14.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  706 Jan 15 13:45 Pregnancy.14.read1_trimmed.fq.gz_bismark_bt2_pe.bam



# back-up data
# WGBS MHL matrix(update by Shicheng) were save in (UTH):
/home/sguo/monod/hap/wgbs/All_chromosomes_combined/wgbs.mhl.txt
# RRBS(Batch 1 and 2) MHL matrix(update by Shicheng) were save in (UTH):
/home/sguo/monod/rrbs_kun/rrbs.kun.mhl.txt
/home/sguo/monod/rrbs_encode
cd /oasis/tscc/scratch/k4zhang/MONOD
cp /oasis/tscc/scratch/k4zhang/MONOD/Ecker_Tissue_WGBS/*pl ./
cp /oasis/tscc/scratch/k4zhang/MONOD/tumor_WGBS/*pl ./
TSCC_whole_blood_WGBS_sampleInfo_WGBS-pooled-mld-blocks.txt
cp /oasis/tscc/scratch/k4zhang/MONOD/whole_blood_WGBS/batch_bam2hapInfo2.pl ./
# 2015-12-20
wc -l /home/shg047/db/hg19/hg19_refGene.bed
head /home/shg047/db/hg19/hg19_refGene.Fantom.bed
cd ~/bioin/annotation/
perl /home/shg047/bioin/bin/trim.pl
perl /home/shg047/bioin/bin/select.pl
perl /home/shg047/bioin/bin/unique.pl
less /home/shg047/bioin/bin/unique.pl
# build new hg19_reference, replace enhancer with Fantom Enhancer
cd ~/db/hg19
grep -v Enhancer hg19_refGene.bed > hg19_refGene.2.bed
awk '{print $1,$2,$3,$4,$5,$6,"Enhancer",$8,$9}' OFS="\t" Enhancers.Fantom.hg19.bed > Enhancer.2.Fantom.hg19.bed9
cat Enhancer.2.Fantom.hg19.bed9 hg19_refGene.2.bed > hg19_refGene.Fantom.bed
rm Enhancer.2.Fantom.hg19.bed9
rm hg19_refGene.2.bed

cd ~/monod/methyblock
# loop source files and get all the overlap(A,B) and non-overlap(A)regions
java -jar ~/bin/M2G.jar -s ~/db/hg19/hg19_refGene.Fantom.bed -l WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed -o rlt
# select overlap regions with source by 9th column
perl /home/shg047/bioin/bin/trim.pl rlt 9
perl /home/shg047/bioin/bin/select.pl rlt.true 6
perl /home/shg047/bioin/bin/unique.pl rlt.true

# 2015-12-18
cat GSM1051152_colon.CpG.calls.txt.bedgraph | awk 'NR<2{print $0;next}{print $0| "sort -k1,1 -k2,2n"}' >  GSM1051152_colon.CpG.calls.txt.sort.bedgraph &
samtools view -h -o human.aligned.cleaned.sam song.aligned.cleaned.bam
samtools view -h -o mouse.aligned.cleaned.sam mouse.aligned.cleaned.sam

# 2015-12-17
grep -v enhancer enhancer.12.encode.cell.type.txt | grep -v file | grep -v Chromosome > enhancer.12.encode.cell.types.txt
perl -lane "next if ! /^chr/; print" enhancer.12.encode.cell.types.txt > enhancer.12.encode.cell.types.true.txt
dos2unix enhancer.12.encode.cell.types.true.txt 
perl -lane "s/\s/\t/g" enhancer.12.encode.cell.types.true.txt 
awk '{print $1, $2-1000, $2+1000,$3}' OFS="\t" enhancer.12.encode.cell.types.true.txt > enhancer.12.encode.cell.types.hg18.txt 
perl -lane "next if @F[0]<0, print" enhancer.12.encode.cell.types.hg18.txt > enhancer.12.encode.cell.types.hg18.nonneg.txt 
wc -l enhancer.12.encode.cell.types.hg18.txt
wc -l enhancer.12.encode.cell.types.hg18.nonneg.txt
bedtools sort -i enhancer.12.encode.cell.types.hg18.nonneg.txt > enhancer.12.encode.cell.types.hg18.sort.txt 
perl -lane "print if @F[-1]>0.95" enhancer.12.encode.cell.types.hg18.sort.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
bedtools merge -i enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
lifeOver

# 2015-12-10
/home/shg047/monod/methyblock/wgbs
wc -l /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/All_WGBS_pooled/mappable_bins_hapInfo/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
wc -l wgbs.bed
bedtools intersect -wao -a wgbs.bed -b ~/db/hg19/CpGI.hg19.bed | perl -lane '{print if @F[8]==0}' | sort -u | wc -l 
cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/All_WGBS_pooled/mappable_bins_hapInfo/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
bedtools intersect -v -a WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed -b ~/db/hg19/CpGI.hg19.bed | wc -l

WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed

cd 
/home/shg047/software/broadinstitute-picard-2a49ee2/src/java/picard/sam
/home/songchen/RNASeq/151211/fastq
scp /media/LTS_33T/SeqStore2016/151211_MiSeq/Fastq
git agordon/libgtextutils
/home/shg047/software/bin/fastx_trimmer

scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/STAR_Index_Gencode_75bp  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/hg38.fa  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/gencode.v23.annotation.gtf  shg047@genome-miner.ucsd.edu:/home/shg047/db/song

scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/STAR_Index_Mouse_Gencode  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/GRCm38.fa  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/gencode.vM7.annotation.gtf  shg047@genome-miner.ucsd.edu:/home/shg047/db/song

mouseGenomeIndex=/oasis/tscc/scratch/yaw004/Genome_Indexes/STAR_Index_Mouse_Gencode
mouseReferenceFasta=
mouseReferenceGTF=


cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/All_WGBS_pooled/mappable_bins_hapInfo/no_tumor/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed  /home/shg047/monod/methyblock

source("/home/shg047/monod/phase3/bedEnrichmentTest.R")



/home/shg047/db/tmp/liftOver

perl -lane "next if ! /^chr/; print" enhancer.12.encode.cell.types.txt > enhancer.12.encode.cell.types.true.txt
dos2unix enhancer.12.encode.cell.types.true.txt 
perl -lane "s/\s/\t/g" enhancer.12.encode.cell.types.true.txt 
awk '{print $1, $2-1000, $2+1000,$3}' OFS="\t" enhancer.12.encode.cell.types.true.txt > enhancer.12.encode.cell.types.hg18.txt 
perl -lane "next if @F[0]<0, print" enhancer.12.encode.cell.types.hg18.txt > enhancer.12.encode.cell.types.hg18.nonneg.txt 
wc -l enhancer.12.encode.cell.types.hg18.txt
wc -l enhancer.12.encode.cell.types.hg18.nonneg.txt
bedtools sort -i enhancer.12.encode.cell.types.hg18.nonneg.txt > enhancer.12.encode.cell.types.hg18.sort.txt 
liftOver enhancer.12.encode.cell.types.hg18.sort.txt  ~/bin/hg18ToHg19.over.chain.gz enhancer.12ct.all.hg19.txt  tmp
cp enhancer.12ct.all.hg19.txt ~/db/hg19/enhancer.12ct.encode.all.hg19.txt

perl -lane "print if @F[-1]>0.95" enhancer.12.encode.cell.types.hg18.sort.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
bedtools merge -i enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
liftOver /home/shg047/monod/RRBS_Plasma_Batch2/enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt ~/bin/hg18ToHg19.over.chain.gz enhancer.12ct.encode.hg19.bed tmp

bedtools window -u -w 0 -b enhancer.12ct.all.hg19.txt -a ~/db/hg19/Enhancers.Fantom.hg19.bed | wc -l 



# 2015-12-10
bedtools shuffle -i /home/shg047/monod/methyblock/wgbs.bed -incl ~/db/hg19/hg19_refGene.bed -g ~/db/hg19/hg19.chrom.sizes 
bedtools intersect -wo -a /home/shg047/monod/methyblock/wgbs.bed -b ~/db/hg19/hg19_refGene.bed 

perl bedAnnoEnrich.pl -i /home/shg047/monod/methyblock/wgbs.bed -r ~/db/hg19/hg19_refGene.bed -g ~/db/hg19/hg19.chrom.sizes -n 50

# compute mhl from haploinfo files
perl /oasis/tscc/scratch/shg047/RRBS/src/haplo2hml.pl /oasis/tscc/scratch/shg047/RRBS/batch2/hapinfo > ../rrbs.batch2.mhl.txt

find ./ -type f -name *pl -exec chmod 700


# copy mhl file from tscc to genome-miner
scp rrbs.batch2.mhl.txt shg047@genome-miner.ucsd.edu:/home/shg047/monod/RRBS_Plasma_Batch2

# R 




cd /oasis/tscc/scratch/shg047/RRBS/
cd /home/shg047/db/hg19
cd /oasis/tscc.old/scratch/k4zhang/MONOD
cd /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles
cd /home/ddiep/Plasma_RRBS/BAMfiles


cd /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles
scp shicheng@meangenemachine.dynamic.ucsd.edu:/home/ddiep/Plasma_RRBS/BAMfiles/RRBS-6P24.sorted.clipped.bam ./
scp shicheng@meangenemachine.dynamic.ucsd.edu:/home/ddiep/Plasma_RRBS/BAMfiles/RRBS-6P24.sorted.clipped.bam.bai ./



awk '{print $1}' combined_list_files | sort -u 
 
less /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed

cp /oasis/tscc.old/scratch/k4zhang/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed ./

scp shicheng@meangenemachine.dynamic.ucsd.edu:



/home/shg047/db/hg19/HsGenome19.CpG.positions.txt

scp shg047@genome-miner.ucsd.edu:/home/shg047/db/hg19/HsGenome19.CpG.positions.txt ./



/home/shg047/db/hg19/genomeDistribution
wget -qO- ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gff3.gz \
    | gunzip --stdout - \
    | gff2bed - \
    > annotations.bed

# 2015-12-3

screen -ls | grep Detached | cut -d. -f1 | awk '{print $1}' | xargs kill
screen -ls | grep pts | cut -d. -f1 | awk '{print $1}' | xargs kill


# 2015-12-1


for i in `ls  methylC*.txt.pairwiseR2`
do
cp pairwise.R $i.pairwise.R
R CMD BATCH "--args $i" $i.pairwise.R & 
done

for i in `ls  N37*C*.txt.pairwiseR2`
do
cp pairwise.R $i.pairwise.R
R CMD BATCH "--args $i" $i.pairwise.R & 
done


R CMD BATCH "--args WB_new-born.all_chrs.hapInfo.txt.pairwiseR2" WB_new-born.all_chrs.hapInfo.txt.pairwiseR2.pairwise.R & 





R CMD BATCH "--args N37-Liver.all_chrs.hapInfo.txt.pairwiseR2" $i.pairwise.R & 
R CMD BATCH "--args $i" $i.pairwise.R & 
R CMD BATCH "--args $i" $i.pairwise.R & 

for i in `ls *hapInfo.txt.pairwiseR2`
do
cp pairwise.R $i.pairwise.R
R CMD BATCH "--args $i" $i.pairwise.R &
done



for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a en.bed -b $i >> en.encode.share.bed
done



for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a en.bed -b $i >> en.encode.share.bed
done

for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a ec.bed -b $i >> ec.encode.share.bed
done

for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a me.bed -b $i >> me.encode.share.bed
done


awk '{print $7}' en.encode.share.bed | sort -u > en.encode.tf.txt
awk '{print $7}' ec.encode.share.bed | sort -u > ec.encode.tf.txt
awk '{print $7}' me.encode.share.bed | sort -u > me.encode.tf.txt





cd /oasis/tscc/scratch/shg047/N37
samtools merge all.bam `find /basedir/ -name "*myfiles*.bam"`

samtools view -H  Indx01.chr6.rmdup.bam > /oasis/tscc/scratch/shg047/N37/merge.sam.header
 cd /home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles
samtools merge -h /oasis/tscc/scratch/shg047/N37/merge.sam.header /oasis/tscc/scratch/shg047/N37/Indx01.bam Indx01.chr1.rmdup.bam Indx01.chr10.rmdup.bam 

samtools view -H Indx01.chr1.rmdup.bam
samtools view -H Indx01.chr2.rmdup.bam


bam2haploinfo.pl
64 to phred

WGBS_pooled_mappable_bins.chr1.mld_blocks_r2-0.5.bed

1, The refGene.txt file is a database file, and consequently is based on the internal representation.
2, Our internal database representations of coordinates always have a zero-based start and a one-based end. We add 1 to the start before displaying coordinates in the Genome Browser. 
3, If you submit data to the browser in position format (chr#:##-##), the browser assumes this information is 1-based. 
4, If you submit data in any other format (BED (chr# ## ##) or otherwise), the browser will assume it is 0-based
5, Similarly, any data returned by the browser in position format is 1-based, while data returned in BED format is 0-based.
6, Some file formats are 1-based (GFF, SAM, VCF) and others are 0-based (BED, BAM)


cd /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/
samtools view -q 20 methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam chr1:10003741-10003773 > methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.sam

/home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam
/home/shg047/monod/haplo/mergedBam2hapInfo_WGBS_25Jun15.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.chr1.mld_blocks_r2-0.5.bed  /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam > methylC-seq_h1-npc_r1.chr1.10003741.10003773.hapInfo.txt

/home/kunzhang/CpgMIP/MONOD/Data/141216_HiSeqRapidRunSeqCap/mappable_bin_hapInfo
less 141216_SeqCap_tumor_RD10_80up.chr21.hapInfo.txt

 head /home/k4zhang/human_g1k_v37/HsGenome19.CpG.positions.txt
 

grep /home/k4zhang/human_g1k_v37/HsGenome19.CpG.positions.txt chr1\s


# 2015-11-30
/home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined
cd /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/
samtools view -b /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.rmdup.bam chr1:10003741-10003773 > methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam
samtools index methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam

cd 
/home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles
mv methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam /home/shg047/monod/



# 2015-11-20

1, download the data
2, checkphred
3, trim_galore
4, compress to gz
5, fastqc
/home/k4zhang/softwares/STAR_2.3.0e.Linux_x86_64/STAR --readFilesCommand zcat --runThreadN 8  --outSAMstrandField intronMotif   --genomeDir /home/k4zhang/my_oasis_tscc/StarIndex/hg19 --readFilesIn Sample_AL-PGP1ipstoCM-102315-3_S3/AL-PGP
# db for bismarker alignment

cd ~/db/aligndb/hg19/bismark
scp hg047@genome-miner.ucsd.edu:/home/shg047/db/hg19/hg19.fa /home/shg047/db/hg19/meth/bismark/Bisulfite_Genome


bismark --bowtie2 --phred64-quals --fastq  --non_directional
bismark -q --phred33-quals --non_directional --bowtie2 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq

cd /oasis/tscc/scratch/shg047/DennisLo2015/fastq

bismark --bowtie2 --phred64-quals --fastq -L 30 -N 1 /home/shg047/db/hg19/meth/bismark -1 Lymphoma.run2.read1.fq.gz -2 Lymphoma.run2.read2.fq.gz -o ../bam




# fastq dataset
cd /oasis/tscc/scratch/shg047/DennisLo2015/fastq

cd /media/TmpStore1/DennisLo2015/


# FastQC for large numbers of samples
http://www.cureffi.org/2012/08/27/fastqc-for-large-numbers-of-samples/?utm_source=delicious+via+twitterfeed&utm_medium=twitter

# 2015-11-17
backtick · bash cmd ·

set filename but leave passwd empty and then 
ssh-keygen
ssh-copy-id shg047@genome-miner.ucsd.edu
Install perl module

wget https://s3.amazonaws.com/deqiangsun/software/moabs/moabs-v1.3.2.src.x86_64_Linux.data.tar.gz
tar xzvf moabs-v1.3.2.src.x86_64_Linux.data.tar.gz
cd moabs-v1.3.2.src.x86_64_Linux.data
cpan App::cpanminus
cpanm Config::Simple
make
pwd
export PATH=/home/shg047/software/moabs-v1.3.2.src.x86_64_Linux.data/bin/:$PATH

moabs --cf mytestrun.cfg
cd data/
moabs -i wt_r1.fq -i wt_r2.fq -i ko_r1.fq -i ko_r2.fq
mcomp -r wt_r1.bam.G.bed,wt_r2.bam.G.bed -r ko_r1.bam.G.bed,ko_r2.bam.G.bed -m wildtype -m knockout -c comp.wiVar.txt --withVariance 1 -p 4 



cd /oasis/tscc/scratch/shg047/DennisLo2015
Ricky Martin

#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
my @file=glob("*.phred64.fq");
foreach my $file(@file
my ($sample,undef)=split /\.phred64.fq/,$file;
rename $file "$sample.fq";
system("gzip $sample.fq");
print "$file\t$phred\n";
}



for i in `*fq`
do
reformat.sh in=$i out=$i.phred64.fq qin=33 qout=64
done



# 2015-11-18
cp /home/kunzhang/bin/trimFastq.pl ~/bin/
vim  ~/bin/trimFastq.pl
vim  /home/kunzhang/bin/trimFastq.pl
cd /home/zhl002/ZL_LTS33/RNA_seq/111215_PGP1iPstoCM/data
perl ~/bin/trimFastq.pl AL-PGP1ipstoCM-102315-10_S10_L001_R1_001.fastq
perl /home/kunzhang/bin/trimFastq.pl < AL-PGP1ipstoCM-102315-10_S10_L001_R1_001.fastq 0 20


ssh-keygen -t rsa
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDQig1C3Oz7sn8Z7dgYWjJL8oyZH71ul+EqmggMC6CWWRZNmXSx73kovKIDMUhrXTCiAXcQ+pHVr9N7PKewZ0Zd2IMfYWcr4tFC5orfOHL40pXyecm1tMwhFbv6U54iV3aPK5PGgGO8zy/n0Ckbj+gKWD95ySQPeptM4lFJ4dqp3z/TudIKAZJqOWj6VzG6MgUZTNpCXO76Lr8A4NXWmMumbcdK0LylSTBlCBjLfR61znP9sZpXEOZwv/OnNmr0pm7D4IM3ny+6sNXBtnexIVn3C15vRgA3OvyVHCZ89EsOLwzMxRcUU3Xp4KQVPNRmbQJ3Utfy9/uh6TE2EgFEZ+gf shg047@genome-miner2
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCtMUHApnb9L6GvhSgaQpbogwioZXjTpi3gg3yJg4YTrnvnfSyVoRz1ZBHkKYVW/HzGDOfLejVe2Pzbkyi06LoO05TiSTm6/ytV47jOMTMnuaJj/pDXMiGashYfy/Gm8TguYGFQeGtj5f+3V6hIGnIHcyF1YFf8qClDwgRQDv4DilM9FpTGU+HcuOIid9NuWcPEv96vwm4Z3mxjXLqAZS0X6aygU4ge6yGPl0zpl1aWX0lzcl0Q9c7mUZ8Z5elgUimX/Ogl+HhJapnhSNWudAV2UI7QXRUphcXdVjaVFRVAKSX+tC6bQ4vU2pjLOxHanadayM3mbpz+IIYaJFPw1B/n shg047@genome-miner2

vim ~/.ssh/id_rsa.pub

scp shg047@genome-miner.ucsd.edu:/media/TmpStore1/DennisLo2015/HOT227.fq.gz  /oasis/tscc/scratch/shg047/DennisLo2015/fastq/
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC2G7U7lO5uw0a7ft9h1FnmU8ITeP/69A49QpHX+OJay+EE+yL6PGovOG4qbq/BnUXSh78u9FLnZxP/kb2ThAdwVz2Ihent2EhpzL4ka2tUirDe4mePUlqMCnyw9O+huwIkiUKUp+C/Wo7IaFs/bltmLpnz3cANIrrXvZr7r2+6nDtai2eZpCYDI1l2R6+cEGuZBA5fwJuXkX47tYHcD1g/ptFrKJzH8WJ9EurRzM8S8F7xeRJhAl3o57+NXzrGcrepEUHX7xfs+ofM3IHPjgVxpgA9875SjsSF4wvwh8XRzqo4FmF0LAv05LLlSSR9JZ4vwsNrM2KXwRw7CReyb+mr shg047@tscc-login1.sdsc.edu

scp shg047@genome-miner.ucsd.edu:/media/TmpStore1/DennisLo2015/HOT227.fq.gz  /oasis/tscc/scratch/shg047/DennisLo2015/fastq/

scp /media/TmpStore1/DennisLo2015/HOT227.fq.gz shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/DennisLo2015/fastq/
scp *.fq.gz shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/DennisLo2015/fastq/




# 2015-11-17
for i in `ls *fq.gz` 
do 
perl ~/bin/checkphred.pl $i
done

#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
my @file=glob("*fq.gz");
foreach my $file(@file){
my $phred=`system(perl ~/bin/checkphred $file)`
print "$file\t$phred\n";
}


#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
my $file=shift @ARGV;
open F, $file;





#!/usr/bin/bash
for i in `ls *fq.gz`
do
j=`perl ~/bin/checkphred.pl $i`
echo "$i $j"
done



grep ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCG AL-LAMPCRGFP-111015-N2indx16_S40_L001_R1_001.fastq | wc -l
grep CAAGTTCAGCGTGTCCG AL-LAMPCRGFP-111015-N2indx16_S40_L001_R1_001.fastq | wc -l


grep ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCG AL-LAMPCRGFP-111015-N2indx16_S40_L001_R1_001.fastq | cut -c 1-60 > keep.fa

perl ~/bin/fa2fastq.pl keep.fa > keep.fastq




# 2015-11-16
@HISEQ:565:C6NCYACXX:5:1101:1220:1868 1:N:0:
NGTTAAATGTTAGGTGAATTTTAGTTTTATTGTTGAGTGAGATTTTTTGCGTTTAGTGTTGTTTTTATATTAGTT
+
#1=DDD>DF?FFDFFFFFIIFIIIIICACFHHHHBDDEFGF>GFIFIFE907BFFFFG=FFIIIIFFEEFEE@=;
@HISEQ:565:C6NCYACXX:5:1101:1220:1868 2:N:0:
NGTTAAATGTTAGGTGAATTTTAGTTTTATTGTTGAGTGAGATTTTTTGCGTTTAGTGTTGTTTTTATATTAGTT
+
#1=DDD>DF?FFDFFFFFIIFIIIIICACFHHHHBDDEFGF>GFIFIFE907BFFFFG=FFIIIIFFEEFEE@=;

killall screen
screen -wipe
pkill screen




prepare fa and fastq for illumina adapter and primers
cd ~/software/FastQC/Configuration
perl txt2fastq.pl adapter_list.txt > adapter_primer.fastq
perl txt2fastq.pl contaminant_list.txt >> adapter_primer.fastq
perl txt2fasta.pl adapter_list.txt > adapter_primer.fa
perl txt2fasta.pl contaminant_list.txt >> adapter_primer.fa


cat read.list |awk '{print "trim_galore -q 20 --phred33 -a AGATCGGAAGAGC "$1}' >  Run.trimgalore.single.sh


# 2015-11-05
TopHat remove the part of read names which distinguish read1 and read2, thus, when using paired-end data as single-end with TopHat.

module load cuda
module load biotools

grep AGATCGGAAGAGC

gunzip -c T21.5.read1.fq.gz | head -n 400000 > T21.5.read1.fq.head
gunzip -c T21.5.read2.fq.gz | head -n 400000 > T21.5.read2.fq.head 

trim_galore --paired -a GATCGGAAGAGCA -a2 GCTCTTCCGATCT --retain_unpaired  --trim1  T21.5.read1.fq.head T21.5.read2.fq.head  

grep AGATCGGAAGAGC T21.5.read1.fq.head_unpaired_1.fq
grep GCTCTTCCGATCT T21.5.read2.fq.head_val_2.fq

grep AGATCGGAAGAGC T21.5.read1.fq.head
grep GCTCTTCCGATCT T21.5.read1.fq.head 
grep AGATCGGAAGAGC T21.5.read2.fq.head
grep GCTCTTCCGATCT T21.5.read2.fq.head 


gunzip -c T21.5.read1.fq.gz | grep AGATCGGAAGAGC | head
gunzip -c T21.5.read2.fq.gz | grep GCTCTTCCGATCT 

gunzip -c T21.5.read1.fq.gz | grep AGATCGGAAGAGC
gunzip -c T21.5.read2.fq.gz | grep TCTAGCCTTCTCG


trim_galore support gzip compressed FastQ files
bismark also support gzip compressed FastQ files (-input files can be uncompressed or gzip-compressed (ending in .gz))

gzip vs. bzip2
1, gzip is faster than bzip2 and the compress ratio is about 35% to 40%
2, bzip2 is slower than gzip however the compress ratio is about 34% to 35%
3, gzip is more conveniment in the usage since majority software can take gzip as the input files

# 2015-11-04


# 2015-11-03
cd /home/zhl002/ZL_LTS33/mouse_WGBS/PCR_validation
bismark -q --phred33-quals --non_directional --bowtie2 -n 1 -l 5 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq
bismark2report ./
deduplicate_bismark
coverage2cytosine
bismark2bedGraph
bismark_methylation_extractor


cd /media/TmpStore1/DennisLo2015
try to analysis the RRBS data from encode project and find the difference between encode project and our own data
1, coverage region
2, coverage

# encode project
cd /media/LTS_60T/shg047/encode
bismark -q --phred64-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark wgEncodeHaibMethylRrbsAg09309UwRawDataRep1.fastq &
bismark -q --phred64-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark wgEncodeHaibMethylRrbsAg09309UwRawDataRep2.fastq &
bismark -q --phred64-quals --bowtie2 -n 1 -l 10 /home/shg047/db/aligndb/hg19/bismark a.fastq 

grep TGTGAGATAGGTTTTAGGTGAAG *.fastq
grep GTTGTAAGAGTTAATGGTTTATTATAAT *.fastq
grep TTTGTTTTTTATTAAGTTAGATGTGT  *.fastq
grep GGTGTGAGTTTGTTTTTTTTAGTTTTTA  *.fastq
grep TTAGGAGGGTAAGGAATATTTTAGG *.fastq
grep TGTGAGTTGAAGTAGGAAGGTTTTT *.fastq
grep TTTGTTTGTTGTTTTTGGAGGATA *.fastq
grep GTTGATGTATGTTTTATTAGTAAATAA *.fastq
grep AGAATGTAGGTTTTTGTTTTGTAGA *.fastq
grep ACCCTAAATTATAATCTCAAACACC  *.fastq
grep AATATCTCTAACTCCTTAACCCACC  *.fastq
grep TCCAAATAAAACAATATTTAAAATCTCC   *.fastq
grep AAACAATTTACTTCAACACAACCAA   *.fastq

GCGAGAAGGCTAGTTTGTTTGTTGTTTTTGGAGGATAGAAGAGCGGTTCAGCAGGAATGCCGAGCAACCCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATTAATAACACAAATCAGAAACCCAAGTGACTGACAGAGAAACCC

/home/shg047/db/aligndb/mm9


cd /home/zhl002/ZL_LTS33/mouse_WGBS/PCR_validation
perl checkscore.pl  ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq ZL-MeValidate-enPCR-14-Oct21_S5_L001_R1_001.fastq ZL-MeValidate-enPCR-12-Oct21_S3_L001_R1_001.fastq

grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
grep AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT *fastq

TGGCTATTTTAGCTGTACCTCAGTATGCCTTAGGGTCATCTTGGGTAATGGGTCTAGAGGACACTTTTAGTGCTAAAACCTGAAAGGGATCCAGTGGTGCCAATTTAGCAACCATGCAGGACAGTAGCCCTGCGGGGAATGACTTCT


grep GTAGGGTGGGAAGTGGTTTTTAT *fastq
grep GGTGTGAGTTTGTTTTTTTTAGTTTTTA *fastq
grep TTTGTTTGTTGTTTTTGGAGGATA *fastq
grep 

CGCAGACATTTCTTCCTGCTTCTCTGCTGTTCTCCGTTTTGTTGTCCTCTTCTTAGGCAGCAC
GCGAGAAGGCTAGTTTGTTTGTTGTTTTTGGAGGATAGAAGAGCGGTTCAGCAGGAATGCCGAGCAACCCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATTAATAACACAAATCAGAAACCCAAGTGACTGACAGAGAAACCC

TTGTTTTTAGTATAATTTTTTTATAGTGTGGTTTGGGTGGGGTTTTTAGAAGGTTTTTTTTGTTTTTAGAGTAAGGTTGATTTTTTTTAAGGTATTGTTTGTTTTTGGTTGTTTTTTATTGTTTTGTTTTTTATTAAGTTAGATGTGTGTTGTAGATATTTTTTTTTGTTTTTTTGTTGTTTTTTGTTTTGTTGTTTTTTTTTTAGGTAGTATGTAGGTTTTATTTGAATTGTTTTAGTATTGTTATGTTTTGGAGATTTTAAATATTGTTTTATTTGGAGTTTTAAATGGTTAAATTTTATGTGGTTTGAGGGTTGGTAGAAGTTATGTTTTTTTTTGGTTTTTTTTATTTTTTTTATGGTT
AATGAGTAATTTTAAGATAGTTTTTTAAAATTAAGATTATTTAGTGGTGGTTTTGTTTTGTTGTTTATTTTATGTTTGTGGGTGTTTTTTTTGTATGTATATATTTTTATTATGTATTTGTTTGTTGTTTTTGGAGGATAGAAGAGGGTATTGGAATTTTTTGGATTGGAGTTGTAGATGGTTGTGAGTTATTATATGGGTGTTTGGAGATGAGTTTAGGTTATTTTGAAGAGTAGTTATTGTATTTAATTGTTGAGTTATTTTTTGAGTTTTGGGTGTTTTGAGATAGGTAAGGGTATAGATGAAGATAATTTTTAAAGGTTTTTTGAGTTTTTTGGATTGTTATATGTTTTTT


fastqc -f fastq -t 3 T21.5.read2.fq.gz  # fastqc detect the file format automatically(fastq,bam,sam)
fastqc -t 6 T21.5.read2.fq.gz


trim_galore sequence_file.fastq
trim_galore --phred33 --fastqc --illumina --non_directional --rrbs --paired --dont_gzip -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
trim_galore --phred33 --fastqc --non_directional --rrbs --dont_gzip ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
trim_galore --phred33 --fastqc --non_directional --rrbs --dont_gzip -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq

bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq

samtools tview -p chr17:76728261 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr8:12808189 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr9:92831023 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr5:27124338 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr3:7632484 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa





cd /home/zhl002/ZL_LTS33/mouse_WGBS/PCR_validation
bismark -q --phred33-quals --non_directional --bowtie2 -n 1 -l 5 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq -o ../


samtools tview -p chr9:52407014 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm10/mm10.fa

samtools tview -p chr11:13303365 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm10/mm10.fa


samtools tview -p chrX:115198611 ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
chr11:13303365
chrX:115198611 
AGATCGGAAGAGC


fastq-mcf ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
fastq-stats ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq

fastq-clipper
fastq-join
fastq-multx


trim_galore ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S13_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S11_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S18_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S22_L001_R1_001.fastq

bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq

grep GCGAGAAGGCTAG *fastq


for i in `ls *fastq`
do
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-10-Oct21_S1_L001_R1_001.fastq
done



bismark 

bismark -q --phred64-quals -n 1 -l 40 --non_directional /data/genomes/homo_sapiens/GRCh37/s_1_sequence.txt


# 2015-10-26
ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/ ./

/home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/
/home/k4zhang/my_oasis_tscc/MONOD/Ecker_Tissue_WGBS/BAMfiles
/home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles
/home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/
/home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/
/home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015


# 2015-10-22
# chr17:16955467-16955609 in Colon_primary_tumor.chr17.sorted.clipped.bam
samtools tview -p chr17:16955467 /home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/Colon_primary_tumor.chr17.sorted.clipped.bam  /home/shg047/db/hg19.fa
# chr11:110910907-110910951 in middle-age.chr11.rmdup.bam
samtools tview -p chr11:110910907 /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam /home/shg047/db/hg19.fa
# grep chr11:110910907-110910951 in hapInfo files of middle-age.chr11.hapInfo.txt
grep chr11:110910907-110910951 /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/by_chrosomes/middle-age.chr11.hapInfo.txt
grep chr11:110910907-110910951  /home/shg047/monod/hap/wgbs/All_chromosomes_combined/WB_middle-age.all_chrs.hapInfo.txt

# RRBS
samtools tview -p chr17:47030368-47030542 /home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015/6-P-10.merged.bam /home/shg047/db/hg19.fa
samtools tview -p chr22:37908213-37908261 /home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015/6-P-10.merged.bam /home/shg047/db/hg19.fa
samtools tview -p chr20:44441518-44441605 /home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015/6-P-10.merged.bam /home/shg047/db/hg19.fa

# BSPP
samtools tview -p chr1:117909025 /home/k4zhang/my_oasis_tscc/MONOD/WGBS_BSPP/BAMfiles/6P-10.sorted.clipped.rmdup.bam /home/shg047/db/hg19.fa

samtools tview -p chr1:10003741 /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.rmdup.bam /home/shg047/db/hg19.fa


# achieve regions from methylation mhl files
cd /home/shg047/monod/dec
awk -F'[:-\t]' 'NR>1 {print $1,$2,$3}' OFS="\t" WGBS_methHap_load_matrix_20Oct2015.txt > WGBS_methHap_load_matrix_20Oct2015.bed
# transfer from genome-miner to tscc to extract aml
cd /home/shg047/monod/aml
scp shg047@genome-miner.ucsd.edu:/home/shg047/monod/dec/WGBS_methHap_load_matrix_20Oct2015.bed ./

PileOMeth extract  -q 10 -p 5 -l WGBS_methHap_load_matrix_20Oct2015.bed /home/shg047/db/hg19.fa  /home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/Colon_primary_tumor.chr17.sorted.clipped.bam  -o Colon_primary_tumor.chr17.sorted.clipped
PileOMeth extract  -q 10 -p 5 -l WGBS_methHap_load_matrix_20Oct2015.bed /home/shg047/db/hg19.fa  /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
PileOMeth extract  -q 10 -p 5 -r chr11:110910907-110910908 /home/shg047/db/hg19.fa /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
PileOMeth extract  -q 10 -p 5 -r chr11:72295556-72295573 /home/shg047/db/hg19.fa /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
PileOMeth extract  -q 10 -p 5 -l WGBS_methHap_load_matrix_20Oct2015.bed.head /home/shg047/db/hg19.fa  /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
chr11:72295556-72295573

my $ref="/home/shg047/db/hg19.fa";
my $region="chr22:16050002-17050002";
my $bamfile="methylC-seq_h1_r2.chr22.rmdup.bam";
system("PileOMeth extract -r $region $ref $bam > aa.txt");


cd /home/shg047/monod/hap/wgbs/All_chromosomes_combined
screen -s a
perl /home/shg047/monod/hap/get_avemeth_matrix_16Oct2015.pl  > WGBS_aveBedMeth_load_matrix_20Oct2015

ftp address:  137.189.133.62
ftp 137.189.133.62
Username: plamethy
Password: de$*d@s3
prompt # to turn prompt off
mget -c -r *  ./


mv middle-age.all_chrs.hapInfo.txt WB_middle-age.all_chrs.hapInfo.txt
mv new-born.all_chrs.hapInfo.txt WB_new-born.all_chrs.hapInfo.txt
mv centenarian.all_chrs.hapInfo.txt WB_centenarian.all_chrs.hapInfo.txt

perl ../../get_avemeth_matrix_16Oct2015.pl > WGBS_aveMeth_load_matrix_20Oct2015.txt

wget -m --ftp-user=plamethy --ftp-password='de$*d@s3' ftp://137.189.133.62/




