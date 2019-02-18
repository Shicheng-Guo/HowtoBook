

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/GWAS")
data1<-read.table("RA_GWASmeta_TransEthnic_v2.txt",head=T)
data2<-read.table("eyre_2012_23143596_ra_efo0000685_1_ichip.sumstats.tsv",sep="\t",head=T)
data3<-read.table("stahl_2010_20453842_ra_efo0000685_1_gwas.sumstats.tsv",head=T)


head(data1)
head(data2)
head(data3)

newdata1<-data1[order(data1[,7]),]
newdata2<-data2[order(data2[,6]),]
newdata3<-data3[order(data3[,6]),]

sig1<-newdata1[newdata1$P.val<10^-8,1]
sig2<-newdata2[newdata2$p<10^-8,3]
sig3<-newdata3[newdata3$p<10^-8,3]

head(sig1)
head(sig2)
head(sig3)
length(sig1)
length(sig2[sig2 %in% sig1])/length(sig2)
length(sig3[sig3 %in% sig1])/length(sig3)
