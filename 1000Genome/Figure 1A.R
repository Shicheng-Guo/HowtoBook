# Figure 1A
data<-read.table("DMER-GWAS-SNP-Distance.txt")
pdf("Figure.1.total.dmer.distance.distitbution.pdf")
dis<-data$V8-(data$V2+data$V3)/2
dis<-dis[abs(dis)<10000000]
hist(dis,main="DMER",breaks=100,cex.axis=0.4,col="blue")
dev.off()
