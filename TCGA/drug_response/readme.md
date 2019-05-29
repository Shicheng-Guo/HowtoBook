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
head()

drugResponse$drug_name<-tolower(drugResponse$drug_name)
drugResponse<-subset(drugResponse,therapy_types=="Chemotherapy" & measure_of_response!="")
head(drugResponse)
sort(table(drugResponse$measure_of_response))
write.table(drugResponse,file="pancancer.drugResponse.txt",sep="\t",col.names=T,row.names=F,quote=F)
