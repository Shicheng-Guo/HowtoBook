## Genome-wide Association Study 

### Plan A: miRNA (mature), GWAS and commonSNP150

1000 genome dataset re-trim (from vcf to plink without any filtering)
```
cd ~/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo plink --vcf ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --make-bed --out ./plink/chr$i >> chr$i.job
qsub chr$i.job
done
```

1982 auto-immnue GWAS-SNP and 2016 immune-disease GWAS records(hg38 bed)
```
bedtools intersect -wo -a AutoImmue.GWAS.SNP.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > AutoImmue.GWAS.SNP.hg38.commonSNP.bed
awk '{print $7}' AutoImmue.GWAS.SNP.hg38.commonSNP.bed | sort -u | wc -l 
```
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

56 miRNA, 64 miRNA-SNPs and 195 immune-disease-GWAS-SNP were collected in the study.
```
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $4}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $8}' | sort -u | wc -l
bedtools window -w 500000 -a miRNA.seed.commonSNP150.hg38.bed -b AutoImmue.GWAS.SNP.hg38.commonSNP.bed | awk '{print $16}' | sort -u | wc -l
```
CHB and CHS performance for these 93 SNPs
```

```
CpG-loss-SNP occurred in RA-hypermethylation regions play predictive roles while CpG-gain-SNPs occurred in RA-hypomethylation regions play risk roles
```

```

Plan A: miRNA, allSNP150, CHB and CHS

