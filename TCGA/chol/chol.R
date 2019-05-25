library("readxl")
library("arm")
library("nlme")
library("randomForest")

setwd("/mnt/bigdata/Genetic/Projects/shg047/chol/phase2")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-stDNA")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
data[1:10,1:10]
methdata<-data.matrix(data[,c(7:180)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rowname)
phen=as.numeric(colnames(methdata)) %%2   # 0 control, 1 case
methdata=data.frame(phen,t(methdata))
methdata<-t(na.omit(t(methdata)))
RF <- randomForest(as.factor(phen) ~ ., data=methdata, importance=TRUE,proximity=T)
imp<-RF$importance
imp<-imp[order(imp[,4],decreasing = T),]
topvar<-match(rownames(imp)[1:20],colnames(methdata))
methdata<-data.frame(methdata[,c(1,topvar)])
glmRlt <- bayesglm(phen ~ . ,family=binomial, data = methdata)
Stepmodel <- stepAIC(glmRlt,trace =T,direction="both",steps=2000)
summary(Stepmodel)
OR= log(exp(summary(Stepmodel)$coefficients[2,1]),base = 10)
P = summary(Stepmodel)$coefficients[2,4]


Table2Generator = function(methydata){
  seq.case = which(methydata[,1] ==0)
  seq.control = which(methydata[,1] == 1)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  McaM = apply(methydata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=T))} )
  McoM = apply(methydata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=T))} )
  Pvalue=apply(methydata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case])$p.value)})
  Pvalue=p.adjust(Pvalue,method="fdr")
  #Logistic regression analysis
  library(pROC)
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(dim(methydata)[2] -1 )){
    temp = methydata[,c(1,i+1 )]
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial")
    OR[i] = log(exp(summary(glm.fit)$coefficients[2,1]),base = 10)
    Logistic.P[i] = summary(glm.fit)$coefficients[2,4]
    CI.upper[i] = log(exp(confint(glm.fit)[2,2]),base = 10)
    CI.lower[i] = log(exp(confint(glm.fit)[2,1]),base = 10)
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit,type=c("response"))
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = logistic.rocobj$auc[[1]]
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = logistic.rocdata[seq.max,1]
    Spec[i] = logistic.rocdata[seq.max,2]
    print(paste(i,colnames(methydata)[i]))
  }
  Logistic.P = p.adjust(Logistic.P, method = "fdr")
  options(digits = 2)
  Table = data.frame(McaM, McoM, Pvalue, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  return(Table)
}

Table<-Table2Generator(methdata)
write.table(Table,file="chol.rlt.txt",sep="\t",quote=F,col.names = NA,row.names =T)
