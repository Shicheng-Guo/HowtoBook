source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/meth450Pancancer.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/meth450Pancancer.R")

setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
methdata[1:5,1:5]
phen4<-id2phen4(colnames(methdata))
phen3<-id2phen3(colnames(methdata))
bin<-id2bin(colnames(methdata))
pid<-id2pid(colnames(methdata))
phen<-data.frame(phen4=phen4,phen3=phen3,pid=pid,bin=bin)
exclude<-which(c(phen$bin !=1 & phen$bin !=11))
phen<-phen[-exclude,]
input<-methdata[,-exclude]
Seq<-paste(phen$pid,phen$bin,sep="-")
head(phen)
input[1:5,1:5]

pca <- prcomp(t(beta),center=T,scale = F)
pdf("MCRI.GBM.BUR.PCA_SDEV.pdf")
plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
dev.off()
scores <- data.frame(phen, pca$x[,1:10])
pdf("MCRI.GBM.BUR.PCA_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
pdf("MCRI.GBM.BUR.PCA_2_3.pdf")
plot(x=scores$PC2,y=scores$PC3, xlim=c(min(scores$PC2),max(scores$PC2)),ylim=c(min(scores$PC3),max(scores$PC3)),xlab="PC2",ylab="PC3",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
pdf("MCRI.GBM.BUR.PCA_2_4.pdf")
plot(x=scores$PC2,y=scores$PC4, xlim=c(min(scores$PC2),max(scores$PC2)),ylim=c(min(scores$PC4),max(scores$PC4)),xlab="PC2",ylab="PC4",pch=16,col=as.numeric(as.factor(phen$disease))+1)
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
pdf("MCRI.GBM.BUR.PCA_3_4.pdf")
plot(x=scores$PC3,y=scores$PC4, xlim=c(min(scores$PC3),max(scores$PC3)),ylim=c(min(scores$PC4),max(scores$PC4)),xlab="PC3",ylab="PC4",pch=16,col=as.numeric(as.factor(phen$disease)))
phen$col=as.numeric(as.factor(phen$disease))+1
legend("topright",legend=c("GBM","LGG","Control"),pch=16,col=2:5,bty="n",cex=1)
dev.off()
