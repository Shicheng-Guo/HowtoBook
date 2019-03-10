#####################################################################
#####################################################################
#####################################################################

setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
files=list.files(pattern="*gdc_hg38.txt$",recursive = T)
methdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=T,sep="\t",row.names = 1)
  methdata<-cbind(methdata,temp[,1])
  print(i)
}
colnames(methdata)<-files
save(methdata,file="methdata.pancancer.RData")
save.image("methdata.pancancer.env.RData")

#####################################################################
#####################################################################
#####################################################################

files=list.files(pattern="*.FPKM-UQ.txt$",recursive = T)
rnaseqdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=F,sep="\t",row.names = 1)
  rnaseqdata<-cbind(rnaseqdata,temp[,1])
  print(i)
}
colnames(rnaseqdata)<-files
save(rnaseqdata,file="rnaseqdata.pancancer.RData")
save.image("rnaseqdata.pancancer.env.RData")
