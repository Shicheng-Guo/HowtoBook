cd /gpfs/home/guosa/nc
mkdir ./intersect_between_dmer_and_gwas
cp RA-OA.DMER.GRCH37.bed ./intersect_between_dmer_and_gwas/
cp GWAS-RA-378.GRCH37.bed ./intersect_between_dmer_and_gwas/
cd ./intersect_between_dmer_and_gwas/
mkdir ./Shuffle/
mkdir ./temp/
for i in `ls *RA-OA.DMER.GRCH37.bed`
do
for j in {1..10000}
do
echo \#PBS -N $i.$j  > $i.$j.job
echo \#PBS -o ./temp/ >>$i.$j.job
echo \#PBS -e ./temp/>> $i.$j.job
echo cd $(pwd) >> $i.$j.job
echo bedtools shuffle -i $i -g ~/hpc/db/hg19/hg19.chrom.sizes \> ./Shuffle/$i.$j.shuffle >> $i.$j.job
echo bedtools sort -i ./Shuffle/$i.$j.shuffle  \> ./Shuffle/$i.$j.shuffle.sort  >> $i.$j.job 
echo bedtools closest -a ./Shuffle/$i.$j.shuffle.sort -b GWAS-RA-378.GRCH37.bed \> ./Shuffle/$i.$j.GWAS.Cloest >> $i.$j.job
qsub $i.$j.job
done
done

## Statistic 
bedtools closest -a RA-OA.DMER.GRCH37.bed -b GWAS-RA-378.GRCH37.bed > RA-OA.DMER-GWAS.N378.bed
## R 
file="RA-OA.DMER-GWAS.N378.bed"
data<-read.table(file)
xx1<-abs(as.numeric(as.character(data[,2]))-as.numeric(as.character(data[,7])))
xx2<-abs(as.numeric(as.character(data[,3]))-as.numeric(as.character(data[,7])))
xx2[which(xx1<xx2)]=xx1[which(xx1<xx2)]
len<-c(mean(xx2,na.rm=T))
sum<-c(sum(xx2<50000,na.rm=T)/length(xx2))
# R
setwd("./")
file=list.files(pattern="*GWAS.Cloest")
Length<-c()
Sum<
for(i in 1:length(file)){
  marker<-unlist(strsplit(file[i],"[.]"))[1]
  data<-read.table(file[i])
  xx1<-abs(as.numeric(as.character(data[,2]))-as.numeric(as.character(data[,7])))
  xx2<-abs(as.numeric(as.character(data[,3]))-as.numeric(as.character(data[,7])))
  xx2[which(xx1<xx2)]=xx1[which(xx1<xx2)]
  Length<-c(Length,mean(xx2,na.rm=T))
  Sum<-c(Sum,sum(xx2<50000,na.rm=T)/length(xx2))
}
LengthS1000000<-sample(Length,100000,replace=T)
