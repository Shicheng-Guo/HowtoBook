setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/Gnomad/exome/aloft-exome-rec/annovar")

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

file=list.files(pattern="multianno.csv")
LOF<-c()
for(i in 1:length(file)){
  data<-read.csv(file[i])
  num<-RiskAAC(data)
  summ<-table(num)
  Nsyn6I12<-names(summ)[which(summ>6)]
  frameshift_deletion<-grep('\\bframeshift deletion',data$ExonicFunc.refGene)
  frameshift_insertion<-grep('\\bframeshift insertion',data$ExonicFunc.refGene)
  stop<-grep('stop',data$ExonicFunc.refGene)
  lof<-sort(c(Nsyn6I12,frameshift_deletion,frameshift_insertion,stop))
  write.table(data[lof,1:2],file=paste(file[i],".lof",sep=""),sep="\t",col.names = F,row.names = F,quote=F)
  print(file[i])
}





