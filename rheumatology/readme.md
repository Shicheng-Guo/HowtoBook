Readme:

miRNA (hg38) mature region (seed) download from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3

241 miRNA-SNP were identified in the seed region of 2883 miRNA records.
```
bedtools intersect -wo -a ~/hpc/db/hg38/miRNA.mature.seed.hg38.bed -b ~/hpc/db/hg38/commonSNP150.hg38.bed > miRNA.seed.commonSNP150.hg38.bed
```
