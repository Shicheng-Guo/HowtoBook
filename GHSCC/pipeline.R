#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("~/hpc/project/pmrp/Exom2/annovar")
file=paste(args[1],".9.hg19_multianno.csv",sep="")
data<-read.table(file,head=T,sep=",",check.names = F)
RiskAAC<-function(data){
  risk1<-grep("D",data$SIFT_pred)
  risk2<-grep("D|P",data$Polyphen2_HDIV_pred)
  risk3<-grep("D|P",data$Polyphen2_HVAR_pred)
  risk4<-grep("D",data$LRT_pred)
  risk5<-grep("D",data$MutationTaster_pred)
  risk6<-grep("H|M",data$MutationAssessor_pred)
  risk7<-grep("D",data$FATHMM_pred)
  risk8<-grep("D",data$PROVEAN_pred)
  risk9<-grep("D",data$MetaSVM_pred)
  risk10<-grep("D",data$MetaLR_pred)
  risk11<-grep("D",data$fathmm.MKL_coding_pred)
  risk12<-grep("D",data$M.CAP_pred)
  risk13<-grep("Name",data$gwasCatalog)
  risk14<-grep("Name",data$tfbsConsSites)
  risk15<-grep("Name",data$targetScanS)
  risk16<-grep("Name",data$wgRna)
  rlt=unique(c(risk1,risk2,risk3,risk4,risk5,risk6,risk7,risk8,risk9,risk10,risk11,risk12,risk13,risk14,risk15,risk16))
  return(rlt)
}
Riskloci<-RiskAAC(data)
rlt<-data[Riskloci,]
ID=paste(rlt$Chr,":",dataset$Start,sep="")
rlt<-data.frame(rlt,ID)
write.table(rlt,file=paste(args[1],"anno.txt",sep="."),col.names=F,row.names=F,quote=F)



data<-read.table("/gpfs/home/guosa/hpc/project/pmrp/Exom2/IBD/relationship-2.txt")
ID<-c()
for(i in unique(data[,2])){
  index=which(data[,2] %in% i)
  id=index[which.max(data[index,3])]
  ID=rbind(ID,c(i,as.character(data[id,1]),data[id,3]))
}
data<-write.table(ID,file="/gpfs/home/guosa/hpc/project/pmrp/Exom2/IBD/relationship-3.txt",sep="\t",quote=F)




setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/biobank/clinicaldata")
data<-read.table("GHCC.RA.CCP.uni.R",head=T,sep="\t",row.names=1)
load("fit.RData")
input<-fit$points
input<-input[input[,1]> -5096,]
input<-input[input[,2]> -2596,]
input<-input[input[,2]< 2596,]

temp<-data[match(rownames(input),rownames(data)),]
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
plot(rep(1,50),col=(colfunc(50)), pch=19,cex=2)

temp$CCP[temp$CCP<=25]<-1
temp$CCP[temp$CCP>25]<-2
temp$CRP[temp$CRP<=10]<-1
temp$CRP[temp$CRP>10]<-2
temp$ESR[temp$ESR<=30]<-1
temp$ESR[temp$ESR>30]<-2
temp$RFIGM[temp$RFIGM<=30]<-1
temp$RFIGM[temp$RFIGM>30]<-2
temp$WBC[temp$WBC<=9.16]<-1
temp$WBC[temp$WBC>9.16]<-2
temp$HCT[temp$HCT<=median(temp$HCT,na.rm=T)]<-1
temp$HCT[temp$HCT>median(temp$HCT,na.rm=T)]<-2
temp$HGB[temp$HGB<=median(temp$HGB,na.rm=T)]<-1
temp$HGB[temp$HGB>median(temp$HGB,na.rm=T)]<-2
temp$HGB[temp$HGB<=median(temp$HGB,na.rm=T)]<-1
temp$HGB[temp$HGB>median(temp$HGB,na.rm=T)]<-2
temp$MCH[temp$MCH<=median(temp$MCH,na.rm=T)]<-1
temp$MCH[temp$MCH>median(temp$MCH,na.rm=T)]<-2
temp$MCV[temp$MCV<=median(temp$MCV,na.rm=T)]<-1
temp$MCV[temp$MCV>median(temp$MCV,na.rm=T)]<-2
temp$MOP[temp$MOP<=median(temp$MOP,na.rm=T)]<-1
temp$MOP[temp$MOP>median(temp$MOP,na.rm=T)]<-2

pdf("MDSplot-chg4.pdf") 
par(mfrow=c(3,3))
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$CCP+1,main="CCP"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$ESR+1,main="ESR"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$RFIGM+1,main="RFIGM"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$WBC+1,main="WBC"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$HCT+1,main="HCT"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$HGB+1,main="HGB"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$MCH+1,main="MCH"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$MCV+1,main="MCV"); 
plot(input, xlab="Coordinate 1", ylab="Coordinate 2", pch=16,cex=0.5,col=temp$MOP+1,main="MOP"); 
dev.off()

head(temp)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/biobank/clinicaldata")
library(Rtsne)  
library("ggplot2")  
data<-read.table("GHCC.RA.CCP.uni.R",head=T,sep="\t",row.names=1)
data=bdata(data)
d <- dist(data) # euclidean distances between the rows
rlt<-Rtsne(d, is_distance=T)
pdf("tsne_RA_clinical.pdf")
tsne = as.data.frame(rlt$Y)  
plot(tsne$Y,pch=16,cex=0.5)
dev.off()
load("tsne.RData")
cols <- rainbow(16)
label=kmeans(tsne,16)
plot(tsne, pch=label$cluster,col=label$cluster,cex=0.4)

plot(tsne, type="n")
label=kmeans(tsne,14)
text(tsne, labels=label$cluster,col=label$cluster,cex=0.8)


par(mfrow=c(3,2))
plot(tsne, pch=16,col=varible2color(data$CCP),cex=0.8,main="CCP")
plot(tsne, pch=16,col=varible2color(data$RFIGG),cex=0.8,main="RF-IGG")
plot(tsne, pch=16,col=varible2color(data$RFIGA),cex=0.8,main="RF-IGA")
plot(tsne, pch=16,col=varible2color(data$RFIGM),cex=0.8,main="RF-IGM")
plot(tsne, pch=16,col=varible2color(data$ESR),cex=0.8,main="ESR")
plot(tsne, pch=16,col=varible2color(data$CRP),cex=0.8,main="CRP")
dev.off()
  

g9=which(label$cluster==9)
c9=which(label$cluster !=9)

rlt<-apply(data,2,function(x) t.test(x[g9],x[c9],na.rm=T))
rbind(lapply(rlt,function(x) x$p.value),lapply(rlt,function(x) x$statistic))



varible2color<-function(x){
  x<-x+1
  x[is.na(x)]<-0
  x
}


bdata<-function(temp){
temp$CCP[temp$CCP<=25]<-1
temp$CCP[temp$CCP>25]<-2
temp$CRP[temp$CRP<=10]<-1
temp$CRP[temp$CRP>10]<-2
temp$ESR[temp$ESR<=30]<-1
temp$ESR[temp$ESR>30]<-2
temp$RFIGM[temp$RFIGM<=30]<-1
temp$RFIGM[temp$RFIGM>30]<-2
temp$RFIGG[temp$RFIGG<=18]<-1
temp$RFIGG[temp$RFIGG>18]<-2
temp$RFIGA[temp$RFIGA<=18]<-1
temp$RFIGA[temp$RFIGA>18]<-2
temp$WBC[temp$WBC<=9.16]<-1
temp$WBC[temp$WBC>9.16]<-2
return(temp)
}

