#!/usr/bin/bash

# Download from GEO
for i in {1..22} X Y M
do
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17972/suppl/GSE17972_HUMtg5lib.qmap.chr$i.txt.gz
done

# gunzip with -f parameter
for i in `ls *.gz`
do
gunzip -f $i
done

# download qmap2bedgraph.pl to qmap file fold
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/801acad12ff839ffe977342312a6cb46912f2551/qmap2bedgraph.pl ./

# transfer qmap to bedgraph
for i in `ls GSE17972_HUMtg5lib.qmap*.txt`
do
perl qmap2bedgraph.pl $i > $i.bedGraph &
done

# Download bigwig related softwares
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig ./
fetchChromSizes hg18 > hg18.chrom.sizes
fetchChromSizes hg19 > hg19.chrom.sizes
fetchChromSizes mm9 > mm9.chrom.sizes
fetchChromSizes mm10 > mm10.chrom.sizes

# bedGraph to bigwig
for i in `ls*bedGraph`
do
sort -k1,1 -k2,2n $i > $i.sort.bedGraph
bedGraphToBigWig $i.sort.bedGraph hg19.chrom.sizes $i.bw
rm 
done


