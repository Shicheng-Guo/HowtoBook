cd /gpfs/home/guosa/hpc/db/hg19/H3k27ac

wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k27ac/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me1/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegMarkH3k4me3/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/
wget -r -l 1 -nd -e robots=off --reject jpg,html http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/

rm index*
gunzip *.gz 

/gpfs/home/guosa/hpc/db/hg19/H3k27ac/wgEncodeRegTfbsClusteredWithCellsV3.bed
/gpfs/home/guosa/hpc/db/hg19/H3k27ac/wgEncodeRegDnaseClusteredV3.bed

wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pharmacogenomics/VIP.hg19.bid
awk '{print "chr"$1,$2,$3,$4}' OFS="\t" VIP.hg19.bid > VIP.hg19.bed
awk '{print $4}' OFS="\t" VIP.hg19.bid > VIP.genelist.txt
grep -F -f VIP.genelist.txt ~/hpc/db/hg19/ref

* TFBS cluster should be occured in >=3 cells


data<-read.table("VIP.genelist.txt",head=F)
db<-read.table("~/hpc/db/hg19/refGeneV2.hg19.bed",head=F)
newdb<-subset(db,V8=="Exon1"| V8=="Enhancer" | V8=="Promoter")
output<-newdb[newdb$V6 %in% data$V1,]
write.table(output,file="VIP.Regulatory.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)

bedtools sort -i VIP.Regulatory.hg19.bed > VIP.Regulatory.hg19.sort.bed
mv VIP.Regulatory.hg19.sort.bed VIP.Regulatory.hg19.bed

bedtools intersect -b VIP.Regulatory.hg19.bed -a wgEncodeRegTfbsClusteredWithCellsV3.bed | sort -u > VIP.TFBS.hg19.bed
bedtools intersect -b VIP.Regulatory.hg19.bed -a wgEncodeRegDnaseClusteredV3.bed | sort -u > VIP.Dnase.hg19.bed

bedtools intersect -a VIP.Dnase.hg19.bed -b VIP.TFBS.hg19.bed | sort -u > VIP.TFBS.Dnase.hg19.bed


for i in `ls *bigWig`
do 
bigWigToBedGraph $i $i.bedgraph
done


