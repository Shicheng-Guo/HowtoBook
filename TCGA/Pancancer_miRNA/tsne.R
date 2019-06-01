load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/miRNA/data/TCGA-Pancancer.miRNAseq.RData")
library("tsne")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
phen<-id2bin(colnames(miRNA))
colors = rainbow(length(unique(phen)))
names(colors) = unique(phen)
ecb = function(x,y){plot(x,t='n'); text(x,labels=phen, col=colors[phen]) }
tsne_iris = tsne(t(miRNA), epoch_callback = ecb, perplexity=10)
