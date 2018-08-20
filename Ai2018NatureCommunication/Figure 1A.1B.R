# Figure 1A
data<-read.table("DMER-GWAS-SNP-Distance.txt")
pdf("Figure.1.total.dmer.distance.distitbution.pdf")
dis<-data$V8-(data$V2+data$V3)/2
dis<-dis[abs(dis)<10000000]
hist(dis,main="DMER",breaks=100,cex.axis=0.4,col="blue")
dev.off()
# Figure 1B
data<-read.table("saturation.txt",head=T,sep="\t")
data<-subset(data,Distance<200000)
pdf("Figure.1B.sauration.pdf")
plot(data[,2]~data[,1],col="blue",type="o",pch=21,lwd=3,cex=2,xlab="distance",ylab="percentage of total GWAS-SNP")
dev.off()
