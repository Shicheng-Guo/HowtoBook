
methylation hm450k datase

setwd("/media/Home_Raid1/shg047")

data<-read.table("./db/hg19/Normal.PBMC.GEO.HM450K.Beta.txt",sep="\t")

input<-read.table("test.txt",sep="\t")

output<-data[match(input[,1],data[,1]),]



