wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip wgEncodeRegTfbsClusteredV3.bed.gz
awk '{print $1"\t"$2"\t"$3"\t"$4}' wgEncodeRegTfbsClusteredV3.bed > wgEncodeRegTfbsClusteredV3.hg19.bed

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz
gunzip wgEncodeRegTfbsClusteredWithCellsV3.bed.gz

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tfbsConsFactors.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt.gz


TERT启动子区域（chr5:1,278,386-1,301,979）内所有转录因子进行诸葛扫描
chr5	1278386	1301979	TERT

wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make
