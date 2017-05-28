# Epigenomic Analysis of Lung and Liver (DNA Methylation)

[SRP072078 (PRJNA315696)](http://www.ebi.ac.uk/ena/data/view/SRP072141&display=html), [SRP072071(PRJNA315694)](http://www.ebi.ac.uk/ena/data/view/SRP072071&display=html), [SRP072075 (PRJNA315695)](http://www.ebi.ac.uk/ena/data/view/SRP072075&display=html) and [SRP072141 (PRJNA315903)](http://www.ebi.ac.uk/ena/data/view/SRP072141&display=html). 

## Background
Dr. Jenkinson published [Potential energy landscapes identify the information-theoretic nature of the epigenome](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3811.html) at Nature Genetics, 2017. In this study, the DNA methylome of normal and cancer tissue from lung and liver were created with WGBS. 

## Data
Download samplesheet(SRR and SRX) from [ebi](http://www.ebi.ac.uk/ena/data/view/SRP072078&display=html)

Each SRX samples only have one SRR files, therefore, hapinfomergebysrx.pl is not nessessary. 

```
Liver(3) and Lung(3)                        : SRP072078.txt
Brain (2)                                   : SRP072071.txt
lymphocytes (6*2 CD4) and fibroblasts (5*1) : SRP072075.txt
Stem cell (1)                               : SRP072141.txt
* In this world. Only Andrew P Feinberg have money to build WGBS like this way: each sample have biological/technical replicate. 
```

```bash
# smartbismark pipeline from fastq-dump to haplotype analysis
for i in `ls SRP*txt`; 
do
perl ~/bin/smartbismark.pl --input $i --genome hg19 --server TSCC --queue hotel --submit submit
done
# bam2hapinfo
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/oasis/db/hg19/hg19.cut10k.bed > saminfo.txt
perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt submit bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt
```

```bash
# after the hapinfo file build, transfer to genome-niner since the memory is limited in tscc
mkdir /media/Home_Raid1/shg047/NAS1/Jenkinson2017NG/hapinfo
mkdir /media/Home_Raid1/shg047/NAS1/Jenkinson2017NG/samplesheet
scp *hapInfo.txt shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS1/Jenkinson2017NG/hapinfo
scp ~/bin/
scp ~

```

## Supplementary
```
#Another method is download with wget
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP072/SRP072078
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP072/SRP072071
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP072/SRP072075
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP072/SRP072141
```

## Necessary Tools
[Cutadapt](https://github.com/marcelm/cutadapt)
```bash
# Install Trim Galore!
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
mv TrimGalore-0.4.3/trim_galore /usr/bin
```
```
awk -F_ '{print "fastq-dump",$1}' xx | sort -u
```
