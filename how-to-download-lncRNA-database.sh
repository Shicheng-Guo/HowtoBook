#!/usr/bin/sh
# where to download lncRNA databse
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz
tar xzvf gencode.v19.long_noncoding_RNAs.gtf.gz
awk 'NR>5 {print $1,$4,$5,$10}' gencode.v19.long_noncoding_RNAs.gtf > lncRNA.hg19.bed
perl -p -i -e "s/[\";]//g" lncRNA.hg19.bed
