Readme:

miRNA (hg38) mature region (seed) download from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3

241 miRNA-SNP (commonSNP150) were identified in the seed region of 2883 miRNA records.
```
cd ~/hpc/rheumatology/RA/miRNASNP
bedtools intersect -wo -a ~/hpc/db/hg38/miRNA.mature.seed.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > miRNA.seed.commonSNP150.hg38.bed
```
9 miRNA, 17 miRNA-SNPs and 19 RA-GWAS-SNP were collected in the study.
```
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $8}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b RA.GWAS.SNP.hg38.commonSNP.uni.sort.bed | awk '{print $13}' | sort -u | wc -l
```

82 miRNA, 93 miRNA-SNPs and 100 immune-disease-GWAS-SNP were collected in the study.
```
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $8}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $13}' | sort -u | wc -l
```
CHB and CHS performance for these 93 SNPs
```

```

