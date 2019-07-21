```
# CHG1

cd /home/guosa/hpc/project/TCGA

library("TCGAbiolinks")

pid<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/Pid.drugResponse.txt",head=F,sep="\t")
pid<-as.character(pid[,1])

drugResponse<-c()
for(i in pid){
website="https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/tcga.";
file=paste(website,i,".drug.txt",sep="")
idat<-read.table(file,sep="\t",head=T)
drugResponse<-rbind(drugResponse,idat)
print(file)
}
head()

drugResponse$drug_name<-tolower(drugResponse$drug_name)
drugResponse<-subset(drugResponse,therapy_types=="Chemotherapy" & measure_of_response!="")
head(drugResponse)
sort(table(drugResponse$measure_of_response))
write.table(drugResponse,file="pancancer.drugResponse.txt",sep="\t",col.names=T,row.names=F,quote=F)


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/miRNA/data")
load("x.RData")
x<-sort(table(drugName)[table(drugName)>10])
pdf("tcga.pancancer.drugname.pdf")
par(cex.lab=1,cex.axis=1,las=3,mar=c(10,5,5,5))
barplot(x,horiz=F,col=1:9,ylim=c(0,300))
dev.off()
drugName[drugName=="5-fu"]<-"5-FU"
drugName[drugName=="fluorouracil"]<-"5-FU"
drugName[drugName=="5-fluorouracil"]<-"5-FU"
drugName[drugName=="taxotere"]<-"taxol"
drugName[drugName=="temodar"]<-"temozolomide"
```
