
setwd("G:\\")
luad<-read.table("/home/sguo/lung/luad/methylation-expression-regulation.luad.2.txt",head=T,sep="\t")
lusc<-read.table("/home/sguo/lung/lusc/methylation-expression-regulation.LUSC.txt",head=T,sep="\t")

lusc<-lusc[na.omit(match(luad[,2],lusc[,2])),]
luad<-luad[na.omit(match(lusc[,2],luad[,2])),]
lung<-data.frame(luad,lusc)

sig<-subset(lung,beta<0 & pvalue<exp(-6) & beta.1<0 & pvalue.1<exp(-6))
write.table(sig,file="methylation-expression-regulation.LUAD.LUSC.Sig.txt",col.names=T,row.names=F,sep="\t")

diffexp<-read.table("/home/sguo/lung/luad.lusc.paired.differential.expresison.pvalue2.txt",head=T,sep="\t",as.is=T)
sigdiffexp<-subset(diffexp,Pvalue<exp(-5) & Pvalue.1<exp(-5))
genename=sapply(sigdiffexp$GeneName,function(x) unlist(strsplit(x,"[|]"))[1])

rlt<-sig[which(sig[,1] %in% genename),]
output1<-sig[which(sig[,1] %in% names(sort(table(rlt$genename),decreasing=T)[1:10])),]
output2<-sigdiffexp[which(genename %in% names(sort(table(rlt$genename),decreasing=T)[1:10])),]

write.table(output1,"methylation.expression.sig.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(output2,"Different.expression.sig.txt",col.names=T,row.names=F,sep="\t",quote=F)

setwd("/home/sguo/tcga/lung")
target<-read.table("target.txt")


luad<-read.table("/home/sguo/tcga/lung/TCGA.RNAseqV2.diff.luad.paired.and.total.sample.txt",head=T,sep="\t",row.names=1)
lusc<-read.table("/home/sguo/tcga/lung/TCGA.RNAseqV2.diff.lusc.paired.and.total.sample.txt",head=T,sep="\t",row.names=1)
gene1<-sapply(rownames(luad),function(x) unlist(strsplit(x,"[|]"))[1])
gene2<-sapply(rownames(lusc),function(x) unlist(strsplit(x,"[|]"))[1])
order1<-match(target[,1],gene1)
order2<-match(target[,1],gene2)
rlt_luad<-luad[order1,]
rlt_lusc<-lusc[order2,]
rlt_luad
rlt_lusc
