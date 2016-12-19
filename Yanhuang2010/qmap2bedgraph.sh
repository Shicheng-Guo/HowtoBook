#!/usr/bin/R
# qmap to bedgraph and then to bigwig. Yanhuang methylome is based on hg18, then transfer to hg19 and hg38
# Contact: Shicheng Guo
# Version 1.3
# Update: 12/19/2016

for i in `ls *.txt`
do
perl qmap2bedgraph.pl $i > $i.bedgraph
liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.bedgraph.hg19 tmp
liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.bedgraph.hg38 tmp
mv $i.bedgraph $i.bedgraph.hg18
bedGraphToBigWig $i.bedgraph.hg18 ~/work/db/hg18/hg18.chrom.sizes $i.hg18.bw
bedGraphToBigWig $i.bedgraph.hg19 ~/work/db/hg19/hg19.chrom.sizes $i.hg19.bw
bedGraphToBigWig $i.bedgraph.hg38 ~/work/db/hg38/hg38.chrom.sizes $i.hg38.bw
done
 
