args = commandArgs(trailingOnly=TRUE)
cor2bed<-function(cor){
  cor<-as.character(cor)
  unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
 }
sam<-read.table("~/hpc/epimarker/bedgraph/SampleMatrix.txt",head=T,sep="\t")
data=read.table(args[1],head=T,sep="\t",row.names=1,check.names=F)
chr=unlist(strsplit(args[1],split=".tab.matrix.rlt"))
filename=unlist(lapply(strsplit(colnames(data),split="[.]"),function(x) x[1]))
newsam=sam[match(filename,sam[,1]),]
blood=which(newsam[,3]=="blood")
solid=which(newsam[,3]=="tissue")

Mean<-rowMeans(data[,blood],na.rm=T)

bed<-rownames(data)[which(Mean<0.3)]
BED<-c()
for(i in 1:length(bed)){
BED<-rbind(BED,c(cor2bed(bed[i]),bed[i]))
}
output=paste(args[1],"BUR.bed",sep="")
write.table(BED,file=output,sep="\t",quote=F,col.names=F,row.names=F)
