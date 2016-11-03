#!/usr/bin/sh

tar xvf GSE39775_RAW.tar
tar xzvf *.tar.gz
gunzip *.gz

cat GSM1083881.*BED > GSM1083881.Placenta.Placenta3.BED
cat GSM1083880.*BED > GSM1083880.Placenta.Placenta3.BED

cat placenta_rep2_MethylC-seq_chr*.BED > GSM978964.Placenta.Placenta1-1.BED
cat placenta_rep1_MethylC-seq_chr*.BED > GSM1083881.Placenta.Placenta1-2.BED
GSM978965
GSM1083880_placenta2_MethylC-seq_chr20.BED
# hg18 to hg19

