#!/usr/bin/sh
# qmap to bedgraph and then to bigwig. Yanhuang methylome is based on hg18, then transfer to hg19 and hg38
# Contact: Shicheng Guo
# Version 1.3
# Update: 12/19/2016
# interesting: sort hg18 bedgraph will be non-sorted after liftOver to hg19 and hg38

# qmap to bedgraph, liftOver hg18 to hg19 and hg38, finally sort the bedgraph 

 for i in `ls *.txt`
 do
 echo $i
 perl qmap2bedgraph.pl $i > $i.bedgraph
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.bedgraph.hg19 tmp
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.bedgraph.hg38 tmp
 mv $i.bedgraph $i.bedgraph.hg18.sort
 sort -k1,1 -k2,2n $i.bedgraph.hg19 > $i.bedgraph.hg19.sort
 sort -k1,1 -k2,2n $i.bedgraph.hg38 > $i.bedgraph.hg38.sort
 done
 
 # merge bedgraph 

 cat *hg18.sort > GSE17972.YanHuang.hg18.bedgraph 
 cat *hg19.sort > GSE17972.YanHuang.hg19.bedgraph 
 cat *hg38.sort > GSE17972.YanHuang.hg38.bedgraph 

 # bedgraph to bigwig 
 for i in `ls *.bedgraph`
 do
 echo $i
 bedGraphToBigWig $i ~/work/db/hg18/hg18.chrom.sizes $i.bw
 bedGraphToBigWig $i ~/work/db/hg19/hg19.chrom.sizes $i.bw
 bedGraphToBigWig $i ~/work/db/hg38/hg38.chrom.sizes $i.bw
 done
