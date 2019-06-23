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
  rlt=c(risk1,risk2,risk3,risk4,risk5,risk6,risk7,risk8,risk9,risk10,risk11,risk12)
  return(rlt)
}


file=list.files(pattern="*.csv")
rlt<-c()
for(i in 1:length(file)){
  data<-read.csv(file[i])
  Riskloci<-names(which(table(RiskAAC(data))>=6))
  rlt<-rbind(rlt,data[Riskloci,])
  print(i)
}
dim(rlt)
write.table(rlt,file=paste("autism","anno.10748.txt",sep="."),col.names=F,row.names=F,quote=F)
