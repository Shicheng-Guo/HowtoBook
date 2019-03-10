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
