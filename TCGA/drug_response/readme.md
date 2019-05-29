```
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
