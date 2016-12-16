for i in `ls *5mC-P*txt`
do
sort -k1,1 -k2,2n $i >$i.sort
perl tobedgraph.pl $i.sort 
bedGraphToBigWig $i.sort.bedgraph ~/db/hg19/hg19.chrom.sizes $i.bw 
rm $i.sort 
rm $i.sort.bedgraph
done
