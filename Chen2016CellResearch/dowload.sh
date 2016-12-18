
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63183/suppl/GSE63183_RAW.tar

for i in `ls *5mC-P*` 
do
sort -k1,1n -k2,2 $i >$i.sort
perl tobedgraph.pl $i.sort 
bedGraphToBigWig $i.sort.bedgraph ~/db/hg19/hg19.chrom.sizes $i.bw &
done


