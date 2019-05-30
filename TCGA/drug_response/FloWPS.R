ENSG<-unlist(lapply(colnames(train.cv),function(x) unlist(strsplit(x,"[.]"))[1]))
ENSG2Symbol(ENSG)
library("readxl")
setwd("/home/guosa/hpc/temp/Response")
file=list.files(pattern="xlsx")

rg<-c()
for(i in 1:10){
  data=data.frame(read_xlsx(file[1],sheet=i))
  rg<-c(rg,as.character(colnames(data)))
}

for(i in 1:10){
  data=data.frame(read_xlsx(file[2],sheet=i))
  rg<-c(rg,as.character(colnames(data)))
}

write.table(unique(rg),file="drugResponseGene.txt",quote=F,row.names = F)
