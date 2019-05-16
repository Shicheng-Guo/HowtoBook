```
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/meth450Pancancer.R")
load("methdata.pancancer.nomissing.RData")
map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
methdata<-input[na.omit(match(map[which(map$V5 %in% "SOX11"),4],rownames(input))),grep("ESCA",colnames(input))]
methdata<-data.frame(phen=id2bin(colnames(methdata)),t(methdata))
methdata[methdata[,1]==11,1]=0
rlt<-Table2Generator(methdata)
write.table(rlt,file="../SOX11.TCGA.ESCA.txt",sep="\t",col.names=NA,row.names=T)
```
