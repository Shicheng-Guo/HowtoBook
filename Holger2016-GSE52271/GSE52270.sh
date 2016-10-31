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

