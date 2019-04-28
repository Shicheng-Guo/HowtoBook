load("MCRI.GBM.full.RData")
beta<-GBM$beta
phen<-GBM$phen
beta[1:5,1:5]
phen[1:5,1:3]
phen$disease=as.character(phen$disease)
phen$tissue=as.character(phen$tissue)

pdf("MCRI.GBM.BUR.PCA.dataset_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$dataset))+1)
phen$col=as.numeric(as.factor(phen$dataset))+1
legend("topright",legend=c("GSE103659","GSE114534","GSE41826","GSE66351","GSE74486","GSE89707","TCGA-Brain"),pch=16,col=2:8,bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.BUR.PCA.tissue_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$tissue))+1)
phen$col=as.numeric(as.factor(phen$tissue))+1
legend("topright",legend=c("blood","brain"),pch=16,col=2:3,bty="n",cex=1)
dev.off()

pdf("MCRI.GBM.BUR.PCA.gender_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$gender))+1)
phen$col=as.numeric(as.factor(phen$gender))
legend("topright",legend=c("male","female"),pch=16,col=c(2,3),bty="n",cex=1)
dev.off()
