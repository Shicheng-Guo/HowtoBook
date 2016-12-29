 #!/usr/bin/sh
 
 liftOver GSE17972.hg18.bedgraph ~/work/db/hg18/hg18ToHg19.over.chain GSE17972.hg19.bedgraph tmp
 sort -u -k1,1 -k2,2n GSE17972.hg19.bedgraph > GSE17972.hg19.bedgraph.sort
 bedGraphToBigWig GSE17972.hg19.bedgraph.sort ~/work/db/hg19/hg19.chrom.sizes GSE17972.hg19.bw
 track type=bigWig color=0,0,255 visibility=2 maxHeightPixels=128:30:11 smoothingWindow=16 windowingFunction=mean name="PBMC" description="Yanhuang-methylome" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Yanhuang2010/GSE17972.hg19.bw

 liftOver GSE17972.hg18.bedgraph ~/work/db/hg18/hg18ToHg38.over.chain GSE17972.hg38.bedgraph tmp
 sort -u -k1,1 -k2,2n GSE17972.hg38.bedgraph > GSE17972.hg38.bedgraph.sort
 bedGraphToBigWig GSE17972.hg38.bedgraph.sort /media/Home_Raid1/shg047/work/db/hg38/hg38.chrom.sizes GSE17972.hg38.bw
 track type=bigWig color=0,0,255 visibility=2 maxHeightPixels=128:30:11 smoothingWindow=16 windowingFunction=mean name="PBMC" description="Yanhuang-methylome" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Yanhuang2010/GSE17972.hg38.bw

 
 
