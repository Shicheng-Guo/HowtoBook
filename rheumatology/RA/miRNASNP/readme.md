
```
setwd("/gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/All_Target_Locations.hg19.bed")
d1<-read.table("/gpfs/home/guosa/hpc/rheumatology/RA/KEGG_90_RA_GeneList.txt")
d2<-read.table("/gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/Design/Target.miRNA.mature.bed")
data<-c()
file=list.files(pattern="*.txt")
for(i in 1:length(file)){
print(i)
d3<-read.table(file[i])
d4<-d3[d3[,4]%in%d1[,1],]
d5<-d4[d4[,5]%in%d2[,8],]
data<-rbind(data,d5)
}
```
