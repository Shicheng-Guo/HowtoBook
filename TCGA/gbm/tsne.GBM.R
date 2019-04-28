
load("MRCI.GBM.BUR.RData")
beta<-GBMBUR$beta
phen<-GBMBUR$phen
pca <- prcomp(t(beta),center=T,scale = F)
RG1<-rotation[head(order(rotation[,1],decreasing=T),n=10),]
RG2<-rotation[head(order(rotation[,1],decreasing=T),n=10),]
RG<-c(rownames(RG1),rownames(RG1))
heatmapinput<-beta[match(RG,rownames(beta)),]
mydata=t(heatmapinput)
phen=phen$disease
prefix="MCRI.GBM"
library("tsne")
library("ggplot2")
colors = rainbow(length(unique(phen)))
names(colors) = unique(phen)
ecb = function(x,y){ plot(x,t='n',xlab="Coordinate 1",ylab="Coordinate 2"); text(x,labels=phen, col=colors[phen]) }
tsne_rlt = data.frame(tsne(mydata, epoch_callback = ecb, perplexity=50))
colnames(tsne_rlt)<-c("xtsne","ytsne")
chart = ggplot(data.frame(tsne_rlt), aes(xtsne,ytsne))+geom_point(size=1,alpha=1,aes(colour = phen))+ggtitle("tSNE dimensions colored by digit")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file="Tsne-20PC1-20PC2.pdf")



