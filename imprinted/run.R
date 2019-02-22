
######################
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/imprinted/167.imprint.gene.txt")
db<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/refGene.hg19.bed")
rlt<-db[db$V5 %in% data[,1],]
write.table(rlt,file="imprinted.hg19.bed",sep="\t",quote=F,col.names = F,row.names = F)
######################
