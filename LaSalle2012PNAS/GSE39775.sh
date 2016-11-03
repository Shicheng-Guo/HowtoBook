#!/usr/bin/sh

tar xvf GSE39775_RAW.tar
tar xzvf *.tar.gz
gunzip *.gz


cat placenta_rep1_MethylC-seq_chr*.BED > GSM978964.Placenta.Placenta1-1.BED
cat placenta_rep2_MethylC-seq_chr*.BED > GSM978965.Placenta.Placenta1-2.BED
cat GSM1083880_*BED > GSM1083880.Placenta.Placenta2.BED
cat GSM1083881_*BED > GSM1083881.Placenta.Placenta3.BED
cat kidney_MethylC-seq_chr*.BED > GSM978966.Kidney.BED
cat cerebellum_MethylC-seq_chr*.BED > GSM978968.Brain.BED
cat NK_MethylC*.BED > GSM978967.NK.BED

rm *chr*

for i in `ls *BED`
do
awk '{print $1,$2,$3,$4}' $i | sort -k1,1 -k2,2n > $i.sort
done

# hg18 to hg19

