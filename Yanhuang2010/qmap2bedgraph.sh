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
 done
 
 # merge liftOver 
 for i in `ls *.bedgraph`
 do
 echo $i
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.hg19 tmp
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.hg38 tmp
 done
 
 # Sort liftOver bedgraph
 for i in `ls *.bedgraph.hg19`
 do
 echo $i
 sort -k1,1 -k2,2n $i > $i.sort
 done

 for i in `ls *.bedgraph.hg38`
 do
 echo $i
 sort -k1,1 -k2,2n $i > $i.sort
 done

 # merge bedgraph by chrosome
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
