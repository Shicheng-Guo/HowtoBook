## Senarior-A
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/miRNA/data/TCGA-Pancancer.miRNAseq.RData")
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/miRNA/data/TCGA-Pancancer.miRNAseq.RData")
library("tsne")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
phen<-id2bin(colnames(miRNA))
colors = rainbow(length(unique(phen)))
names(colors) = unique(phen)
ecb = function(x,y){plot(x,t='n'); text(x,labels=phen, col=colors[phen]) }
pdf("pancancer.miRNA.tnse.pdf")
par(mfrow=c(2,2))
tsne_iris = tsne(t(miRNA), epoch_callback = ecb, perplexity=10)
tsne_iris = tsne(t(miRNA), epoch_callback = ecb, perplexity=30)
tsne_iris = tsne(t(miRNA), epoch_callback = ecb, perplexity=50)
tsne_iris = tsne(t(miRNA), epoch_callback = ecb, perplexity=70)
dev.off()


## Scenario-B
clic<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/clinical.txt",head=T,sep="\t")
pid<-clic[match(id2phen3(colnames(miRNA)),clic[,1]),2]
pid<-pid[grep("TCGA",pid)]
input<-miRNA[,grep("TCGA",pid)]
colors = rainbow(length(unique(pid)))
names(colors) = unique(phen)
ecb = function(x,y){plot(x,t='n'); text(x,labels=phen, col=colors[phen]) }
pdf("pancancer.miRNA.tnse.pid.pdf")
tsne_iris10 = tsne(t(miRNA), epoch_callback = ecb, perplexity=10,max_iter = 100)
dev.off()


## Scenario-C


