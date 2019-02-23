for(i in 1:22){
div<-read.table(paste("chr",i,".CHB.hap.txt",sep=""))
map<-read.table(paste("chr",i,".CHB.hap.map.txt",sep=""))
pdf(paste("chr",i,".CHB.pdf",sep=""))
plot(div[,1]~map[,2],cex=0.35,col="blue",ylab="number of haplotypes within 500bp",main=paste("chr",i,sep=""),xlab="Position in Chromosome")
dev.off()
}

rlt<-c()
for(i in 1:22){
  div<-read.table(paste("chr",i,".CHB.hap.txt",sep=""))
  map<-read.table(paste("chr",i,".CHB.hap.map.txt",sep=""))
  loci<-which(div[,1]>100)
  rlt<-rbind(rlt,data.frame(map[loci,],div[loci,]))
}

write.table(rlt,file="highDiversityHaplotype.txt",sep="\t",quote=F,col.names = F,row.names = F)

awk '{print "chr"$1,$2,$3,$4}' OFS="\t" highDiversityHaplotype.txt >  highDiversityHaplotype.hg19.bed
/gpfs/home/guosa/hpc/db/hg19/1000Genome/highDiversityHaplotype.hg19.bed
bedtools intersect -wao -a highDiversityHaplotype.txt  -b 68194.GWASCatlog.hg19.bed


##################################################################################

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/1000Genome")

d1<-read.table("chr1.CHB.hap.txt")
d2<-read.table("chr1.CEU.hap.txt")
d3<-read.table("chr1.YRI.hap.txt")

png("hapcount500.png")
par(mfrow=c(3,1))
plot(d1[,1],type="l",col="red",ylab="# haplotype",main="CHB-N=103",xlab="chr1-500bp-window")
plot(d2[,1],type="l",col="blue",ylab="# haplotype",main="CEU-N=99",xlab="chr1-500bp-window")
plot(d2[,1],type="n",ylab="# haplotype",,main="CHB+CEU",xlab="chr1-500bp-window")
lines(d1[,1],type="l",col="red")
lines(d2[,1],type="l",col="blue")
dev.off()

png("hapcount5000.png")
par(mfrow=c(3,1))
plot(d1[,1],type="l",col="red",ylab="# haplotype",main="CHB-N=103",xlab="chr1-10000bp-window")
plot(d2[,1],type="l",col="blue",ylab="# haplotype",main="CEU-N=99",xlab="chr1-10000bp-window")
plot(d3[,1],type="l",col="green",ylab="# haplotype",main="YRI-N=108",xlab="chr1-10000bp-window")
dev.off()



d1<-read.table("chr1.CHB.hap.txt")
d2<-read.table("chr1.CEU.hap.txt")
d3<-read.table("chr1.YRI.hap.txt")

quantile(d1[,1])

max(d1[,1])/206
max(d2[,1])/198
max(d3[,1])/216
png("hapcount5000.png")
par(mfrow=c(3,1))
plot(d1[,1]/206,type="l",col="red",ylab="# haplotype",main="CHB-N=103",xlab="chr1-5000bp-window",ylim=c(0,1))
plot(d2[,1]/198,type="l",col="blue",ylab="# haplotype",main="CEU-N=99",xlab="chr1-5000bp-window",ylim=c(0,1))
plot(d3[,1]/216,type="l",col="blue",ylab="# haplotype",main="YRI-N=108",xlab="chr1-5000bp-window",ylim=c(0,1))
dev.off()





median(d1[d1>0])
mean(d1[d1>0])
plot(density(d1[d1>0]),xlim=c(0,10))

quantile(d1[d1>0])
mean
