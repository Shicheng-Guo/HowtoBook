# 2016-12-30
cd /media/Home_Raid1/shg047/work/db/hg19
perl rmsk.pl -i rmsk.hg19 -o rmsk.hg19.bed
sort -u rmsk.hg19.bed > rmsk.hg19.bed.temp
sort -k1,1 -k2,2n rmsk.hg19.bed.temp > rmsk.hg19.bed
rm rmsk.hg19.bed.temp

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
done

perl ~/bin/bigWigAverageOverBed2Matrix.pl > Repeat.txt

data<-read.table("Repeat.txt")
info<-read.table("/media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed")
group<-unlist(lapply(info[,4],function(x) unlist(strsplit(as.character(x),":"))[3]))
input<-data.frame(data,group)

library("reshape2")
colnames(input)<-c("T","N","T","N","group")
input.long<-melt(input, id.vars=c("group"))

library(ggplot2)
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =1,outlier.colour="white")+ coord_flip()+ geom_point(position = position_jitter(width = 0.2))
dev.off()


# it is so strange,therefore, I change the methylation status of kidney from Roadmap sample(GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig.gz.bw)
cd /media/Home_Raid1/shg047/work/Roadmap/bw
for i in GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig.gz.bw GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.gz.bw
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
done
perl ~/bin/bigWigAverageOverBed2Matrix.pl > Repeat.txt

data<-read.table("Repeat.txt")
info<-read.table("/media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed")
group<-unlist(lapply(info[,4],function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(data,group)
library("reshape2")
colnames(input)<-c(c("Kidney","Esophagus"),"group")
input.long<-melt(input, id.vars=c("group"))
library(ggplot2)
pdf("rmsk.methylation.pdf")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA,outlier.colour="white")+ coord_flip()
dev.off()
input.long<-data.frame(input.long)
tapply(input.long$value,input.long$group,function(x) median(x))



