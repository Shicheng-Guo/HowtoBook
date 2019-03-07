Genome-wide epigenetic meta-analysis identifies SHOX2 associated with multiple autoimmune disease

DRB1*15:01



 	HLA-DRB1	HDAC-INHIBITOR	NASID
RA	Yes	Yes	Yes
SLE	No	No	Yes
SSC	No	No	No
SS	No	No	Yes
AAV	Yes	No	No
BD	Yes	No	No
AS	Yes	No	Yes
PSA	Yes	No	Yes
OA	No	No	Yes
GOUT	No	No	Yes
IBD	Yes	No	No
CD	Yes	No	No


https://www.nature.com/articles/s12276-019-0215-5


```
data1<-read.table("/gpfs/home/guosa/hpc/db/hg19/H3k27ac/VIP.genelist.txt",head=F)
gene1<-unique(data[,1])
data2<-read.table("/gpfs/home/guosa/hpc/GWAS_Catalog/immnue.Gene.txt")
gene2<-names(which(table(data2[,1])>10))

db<-read.table("~/hpc/db/hg19/refGeneV2.hg19.bed",head=F)
newdb<-subset(db,V8=="Exon1"| V8=="Enhancer" | V8=="Promoter")

output1<-newdb[newdb$V6 %in% gene1,]
output2<-newdb[newdb$V6 %in% gene2,]

output<-unique(rbind(output1,output2))
unique(output$V6)
write.table(output,file="VIP.PID.Regulatory.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)

bedtools sort -i VIP.PID.Regulatory.hg19.bed > VIP.PID.Regulatory.hg19.sort.bed
mv VIP.PID.Regulatory.hg19.sort.bed VIP.PID.Regulatory.hg19.bed

bedtools intersect -b VIP.PID.Regulatory.hg19.bed -a wgEncodeRegTfbsClusteredWithCellsV3.bed | sort -u > VIP.TFBS.hg19.bed
bedtools intersect -b VIP.PID.Regulatory.hg19.bed -a wgEncodeRegDnaseClusteredV3.bed | sort -u > VIP.Dnase.hg19.bed
bedtools intersect -a VIP.Dnase.hg19.bed -b VIP.TFBS.hg19.bed | sort -u > VIP.TFBS.Dnase.hg19.bed
bedtools sort -i VIP.TFBS.Dnase.hg19.bed >  VIP.TFBS.Dnase.hg19.sort.bed
bedtools merge -d 500 -i VIP.TFBS.Dnase.hg19.sort.bed  | awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" > VIP.TFBS.Dnase.hg19.merge.bed

bedtools intersect -wo -a VIP.TFBS.Dnase.hg19.merge.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed | grep PTPN

# bw2tab
mkdir temp
for i in `ls *bigWig`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job 
echo bigWigAverageOverBed $i VIP.TFBS.Dnase.hg19.merge.bed $i.tab >> ./temp/$i.job
qsub ./temp/$i.job
done


tab2matrix.pl > Matrix.txt

data<-read.table("Matrix.txt")
data$RowMean<-rowMeans(data)
write.table(data,file="Matrix.RowMean.bed",sep="\t",quote=F,row.names=T,col.names=NA)

awk '{print $1,$23}' OFS="\t" Matrix.RowMean.bed | grep -v Encod > Matrix.RowMean.bid
perl -p -i -e 's/:/\t/' Matrix.RowMean.bid
perl -p -i -e 's/-/\t/' Matrix.RowMean.bid

bedtools intersect -wo -a Matrix.RowMean.bid -b ~/hpc/db/hg19/refGeneV2.hg19.bed > Matrix.RowMean.Symbol.bid

data<-read.table("Matrix.RowMean.Symbol.bid",head=F)
newdata<-subset(data,V11=="Exon1"| V11=="Enhancer" | V11=="Promoter")
dim(newdata)

output<-c()
for(i in unique(data$V10)){
input<-unique(subset(data,V10==i))
new<-unique(input[,1:4])
if(nrow(new)>1){
Max<-which.max(new$V4)
output<-rbind(output,input[Max,])
}else{
output<-rbind(output,input)
}
}
rlt<-unique(output[,1:4])
write.table(rlt,file="VIP.PID.Regulatory.MaxScore-Region.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(output,file="VIP.PID.Regulatory.MaxScore-Gene.hg19.bed",sep="\t",quote=F,row.names=F,col.names=F)
```
