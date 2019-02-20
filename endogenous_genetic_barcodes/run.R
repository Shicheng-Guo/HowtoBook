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
