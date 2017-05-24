## Charting a dynamic DNA methylation landscape of the human genome

### Data Collection and Analysis

PRJNA201480 ebi table download table: [click here](http://www.ebi.ac.uk/ena/data/view/PRJNA201480)

NCBI sra download table: [click here](http://www.ncbi.nlm.nih.gov/sra/?term=SRP028600)

```bash
cd /oasis/tscc/scratch/shg047/Ziller2013/fastq
perl ~/bin/smartbismark.pl --input saminfo.txt --genome hg19 --server TSCC --submit no --queue hotel
# Primary Colon Cancer: SRR949210 SRR949211 SRR949212 
# Adjacent normal tissue: SRR949213 SRR949214 SRR949215
for i in SRR949210 SRR949211 SRR949212 SRR949213 SRR949214 SRR949215
do
qsub $i.pbs
done
```
Recover bigwig files from Supplementary bedgraph from GEO database
### 
```{bash}
cd /home/shg047/oasis/Ziller2013/bw
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46644/suppl/GSE46644_bedFiles_set1.tar.gz  
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46644/suppl/GSE46644_bedFiles_set2.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46644/suppl/GSE46644_bedFiles_set3.tar.gz
tar xzvf GSE46644_bedFiles_set1.tar.gz  
tar xzvf GSE46644_bedFiles_set2.tar.gz  
tar xzvf GSE46644_bedFiles_set3.tar.gz  

# not same as GSM files, therefore download GSM again
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1204nnn/GSM1204465/suppl/GSM1204465_BiSeq_cpgMethylation_BioSam_1120_Colon_Primary_Tumor.BiSeq.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1204nnn/GSM1204466/suppl/GSM1204466_BiSeq_cpgMethylation_BioSam_1121_Colon_Adjacent_Normal.BiSeq.bed.gz

# change Ziller's bed file to bedgraph and bigwig
for i in `ls *bed`
do
awk -F"\t|\'|\/" '{print $1,$2,$3,$5/$6,$5,$6}' OFS="\t" $i > $i.bedgraph
bedGraphToBigWig $i.bedgraph ../../db/hg19/hg19.chrom.sizes $i.bw
done



```

### Description

Paired tumor and normal colon whole-genome bisulfite sequencing (WGBS) data. More information about this data is available on GEO ([GSE46644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46644)). 

### Citation

Ziller MJ, Gu H, Müller F, Donaghey J et al. Charting a dynamic DNA methylation landscape of the human genome. Nature 2013 Aug 22;500(7463):477-81. PMID: [23925113](https://www.ncbi.nlm.nih.gov/pubmed/23925113)

