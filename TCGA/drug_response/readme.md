```
library("TCGAbiolinks")
pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]
drugResponse<-c()
for(i in pid){
website="https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/tcga.";
file=paste(website,i,".drug.txt",sep="")
idat<-read.table(file,sep="\t",head=T)
drugResponse<-rbind(drugResponse,idat)
print(file)
}
