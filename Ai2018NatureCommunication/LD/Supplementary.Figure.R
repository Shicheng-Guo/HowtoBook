######################################################################
######################################################################
pdf("../Figure 1-CEU-CHB-YRI.pdf")
par(mfrow=c(3,2))

input<-"DMER.twoside.dp.max.CEU.txt"
data<-read.table(input,fill=T,row.names=1)
p<-c()
q<-c()
for(i in seq(0,1,0.05)){
  p=c(p,sum(data[,1]>i))
  q=c(q,sum(data[,2]>i,na.rm = T))
}

plot(p~seq(0,1,0.05),col="blue",main="Max DP for DMER and Symmetry regions",xlab="D'",ylab="Number of DMERs (CEU)",type="o",pch=16,cex=0.8,cex.main=0.85)
lines(q~seq(0,1,0.05),type="o",col="red",pch=16)
legend("bottomleft",col=c("red","blue"),legend=c("DMER","Symmetry"),lty=1,pch=16,bty = "n")

plot(density(data[,2]),col="red",xlab="D'",main="Distribution of max(D')",xlim=c(0,0.98),lwd=1.5,cex=0.8,cex.main=0.85)
lines(density(data[,1]),col="blue")
legend("topleft",col=c("red","blue"),legend=c("DMER","Symmetry"),lty=1,bty = "n",lwd=2)

input<-"DMER.twoside.dp.max.CHB.txt"
data<-read.table(input,fill=T,row.names=1)
p<-c()
q<-c()
for(i in seq(0,1,0.05)){
  p=c(p,sum(data[,1]>i))
  q=c(q,sum(data[,2]>i,na.rm = T))
}

plot(p~seq(0,1,0.05),col="blue",main="Max DP for DMER and Symmetry regions",xlab="D'",ylab="Number of DMERs (CHB)",type="o",pch=16,cex=0.8,cex.main=0.85)
lines(q~seq(0,1,0.05),type="o",col="red",pch=16)
legend("bottomleft",col=c("red","blue"),legend=c("DMER","Symmetry"),lty=1,pch=16,bty = "n")

plot(density(data[,2]),col="red",xlab="D'",main="Distribution of max(D')",xlim=c(0,0.98),lwd=1.5,cex=0.8,cex.main=0.85)
lines(density(data[,1]),col="blue")
legend("topleft",col=c("red","blue"),legend=c("DMER","Symmetry"),lty=1,bty = "n",lwd=2)


input<-"DMER.twoside.dp.max.YRI.txt"
data<-read.table(input,fill=T,row.names=1)
p<-c()
q<-c()
for(i in seq(0,1,0.05)){
  p=c(p,sum(data[,1]>i))
  q=c(q,sum(data[,2]>i,na.rm = T))
}
plot(p~seq(0,1,0.05),col="blue",main="Max DP for DMER and Symmetry regions",xlab="D'",ylab="Number of DMERs (YRI)",type="o",pch=16,cex=0.8,cex.main=0.85)
lines(q~seq(0,1,0.05),type="o",col="red",pch=16)
legend("bottomleft",col=c("red","blue"),legend=c("DMER","Symmetry"),lty=1,pch=16,bty = "n")
plot(density(data[,2]),col="red",xlab="D'",main="Distribution of max(D')",xlim=c(0,0.98),lwd=1.5,cex=0.8,cex.main=0.85)
lines(density(data[,1]),col="blue")
legend("topleft",col=c("red","blue"),legend=c("DMER","Symmetry"),lty=1,bty = "n",lwd=2)
dev.off()

data.frame(seq(0,1,0.05),p,q)





