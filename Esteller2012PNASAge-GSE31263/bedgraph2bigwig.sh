#!/usr/bin/sh

for i in `ls *bedGraph`
do
sort -k1,1 -k2,2n $i > $i.sort
bedGraphToBigWig $i.sort ~/NAS3/db/hg19/hg19.chrom.sizes $i.sort.bw
echo $i
done
