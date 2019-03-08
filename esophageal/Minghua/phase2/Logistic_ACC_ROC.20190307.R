if (!require("pROC")) biocLite("pROC")
if (!require("readxl")) biocLite("readxl")
library(pROC)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/esophageal_carcinoma/phase2/profile")

data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
head(rowname)

data[1,]
methdata<-data.matrix(data[,c(12:ncol(data))])
methdata[1,]

rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
head(genesymbol)

phen=rep(0,ncol(methdata))  
phen[grep("T",colnames(methdata))]<-1
phen

methdata=data.frame(phen,t(methdata))
methdata[1:5,1:5]

methdata=ColNARemove(methdata)
methdata[1:5,1:5]
genesymbol= unlist(lapply(colnames(methdata)[2:ncol(methdata)], function(x) strsplit(as.character(x),"_")[[1]][1]))
head(genesymbol)

rlt<-Table2Generator(methdata)
genesymbol= unlist(lapply(rownames(rlt), function(x) strsplit(as.character(x),"_")[[1]][1]))
rlt<-data.frame(genesymbol,rlt)

write.table(rlt,file="Phase2.ESCC.Sen.Spe.Acc.txt",col.names=NA,row.names = T,quote=F,sep="\t")
xueqing<-rlt[rlt$genesymbol %in% c("SOX11","LINE.1","SHOX2"),]
write.table(xueqing,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.txt",col.names=NA,row.names = T,quote=F,sep="\t")

rlt1<-combineAUC(methdata,c("LINE.1","SHOX2"))
rlt2<-combineAUC(methdata,c("LINE.1","SOX11"))
rlt3<-combineAUC(methdata,c("SHOX2","SOX11"))
rlt4<-combineAUC(methdata,c("LINE.1","SHOX2","SOX11"))

write.table(rlt1$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt1.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt1$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt1.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt2$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt2.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt2$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt2.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt3$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt3.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt3$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt3.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt4$matrix,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt4.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt4$model,file="Phase2.ESCC.Sen.Spe.Acc.xueqing.rlt4.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

rlt6<-combineAUC(methdata,recombination=c("LINC00466","LncRNA"))
rlt7<-combineAUC(methdata,recombination=c("LINC00466","MIR129"))
rlt8<-combineAUC(methdata,recombination=c("MIR129","LncRNA"))
rlt9<-combineAUC(methdata,recombination=c("LINC00466","MIR129","LncRNA"))

write.table(rlt6$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt6.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt6$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt6.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt7$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt7.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt7$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt7.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt8$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt8.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt9$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt8.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

write.table(rlt9$matrix,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt9.matrix.txt",col.names=NA,row.names = T,quote=F,sep="\t")
write.table(rlt9$model,file="Phase2.ESCC.Sen.Spe.Acc.LINC00466.rlt9.AUC.txt",col.names=F,row.names = T,quote=F,sep="\t")

plot(rlt6$roc,col=3)
lines(rlt7$roc,col=4)
lines(rlt8$roc,col=5)
lines(rlt9$roc,col=2)
legend("bottomright",cex=0.9,legend=c("LINC00466+TCONS_00021941","LINC00466+MIR129-2","MIR129-2+TCONS_00021941","LINC00466+MIR129-2+TCONS_00021941"),col=c(3,4,5,2),lty=1,bty="n")

rlt10<-combineAUC(methdata,recombination=c("LINC00466"))
rlt11<-combineAUC(methdata,recombination=c("SOX11"))
rlt12<-combineAUC(methdata,recombination=c("LINE.1"))
rlt13<-combineAUC(methdata,recombination=c("SHOX2"))

pdf("SOX11.LINC00466.pdf")
par(mfrow=c(2,2))
plot(rlt10$roc,col=3,main="LINC00466")
text(x=0.25,y=0.25,paste("AUC=",round(rlt10$model[3],3),sep=""))
plot(rlt11$roc,col=3,main="SOX11")
text(x=0.25,y=0.25,paste("AUC=",round(rlt11$model[3],3),sep=""))
plot(rlt12$roc,col=3,main="LINE-1")
text(x=0.25,y=0.25,paste("AUC=",round(rlt12$model[3],3),sep=""))
plot(rlt13$roc,col=3,main="SHOX2")
text(x=0.25,y=0.25,paste("AUC=",round(rlt13$model[3],3),sep=""))
dev.off()

AUC<-c()
GENE<-unlist(lapply(unique(genesymbol),function(x) unlist(strsplit(x,"[-]"))[1]))
for(i in GENE){
  print(i)
  rlt<-combineAUC(methdata,recombination=c(i))
  AUC<-rbind(AUC,rlt$model)
}
rownames(AUC)<-unique(genesymbol)
write.table(AUC,file="Phase2.ESCC.Sen.Spe.Average.AUC.txt",col.names=NA,row.names = T,quote=F,sep="\t")


combineAUC(methdata,recombination=c("SHOX2","CDO1"))
combineAUC(methdata,recombination=c("SHOX2"))

recombination=c("SHOX2","CDO1")
recombination=c("SHOX2")


combineAUC<-function(methdata,recombination="."){
  Table<-list()
  temp <- methdata[,grepl(paste(recombination, collapse="|"), colnames(methdata))]
  genesymbol= unlist(lapply(colnames(temp), function(x) strsplit(as.character(x),"_")[[1]][1]))
  temp<-t(apply(temp,1,function(x) tapply(x, genesymbol,function(x) mean(x,na.rm=T))))
  head(temp)
  if(nrow(temp)==1){
  temp<-t(temp)
  colnames(temp)<-unique(genesymbol)
  }
  head(temp)
  phen=rep(0,nrow(temp))  
  phen[grep("T",rownames(temp))]<-1
  temp=data.frame(phen,temp)
  head(temp)
  glm.fit  = glm(phen~ ., data = temp, family = "binomial")
  summary(glm.fit)
  logOR = log(exp(summary(glm.fit)$coefficients[,1]),base = 10)
  Logistic.P = summary(glm.fit)$coefficients[,4]
  CI.upper=log(exp(confint(glm.fit)[,2]),base = 10)
  CI.lower = log(exp(confint(glm.fit)[,1]),base = 10)
  #Do the analysis of the sens, spec, and AUC
  predicted.value = predict(glm.fit)
  predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
  logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
  logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
  AUC = logistic.rocobj$auc[[1]]
  #Find the best Sens and Spec
  logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
  seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
  Sens = logistic.rocdata[seq.max,1]
  Spec = logistic.rocdata[seq.max,2]
  Table$matrix = data.frame(logOR, CI.upper, CI.lower, Logistic.P)
  Table$model=c(Sen=Sens,Spec=Spec,AUC=AUC)
  Table$roc=logistic.rocobj
  return(Table)
}


Table2Generator = function(methdata){
  seq.case = which(methdata[,1] ==1)
  seq.control = which(methdata[,1] == 0)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methdata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methdata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case],na.rm=T)$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  #Logistic regression analysis
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(ncol(methdata)-1 )){
    temp = methdata[,c(1,i+1 )]
    temp[,1] = ifelse(temp[,1] ==1,1,0)
    temp[,1] = as.factor(temp[,1])
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i]=log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}


options(digits = 2)
index2type<-function(index){
  sampletype=ifelse(as.numeric(index) %% 2,"Case","Normal") # for chol project, odds is case while even is control
}
data2summary <- function(data, varname, groupnames){
  # require(plyr)
  # c(mean(x)-2*sem,mean(x)+2*sem)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem=sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])),
      iqr=as.numeric(quantile(x[[col]],na.rm=T)[4]-quantile(x[[col]],na.rm=T)[2]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


