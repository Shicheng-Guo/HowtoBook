install.packages("openxlsx", dependencies = TRUE)
install.packages("readxl", dependencies = TRUE)
library("openxlsx")
library("readxl")
file=list.files(pattern="frontier_research_topic*")
MFR<-c()
for(i in 1:length(file)){
  data = read_excel(file[i],sheet = 1)
  data = as.data.frame(data)
  MFR<-rbind(MFR,data)
}
MFR
MFR[,2]<-"Dr"
MFR[,6]<-"Research"
MFR[,2]<-"Dr"
M1<-which(is.na(MFR[,3]))
MFR[M1,3]<-MFR[M1,5]
MFR[,4]<-""
write.table(MFR,file="../AuthorList.txt",sep="\t",quote=F,col.names = T,row.names = F)
