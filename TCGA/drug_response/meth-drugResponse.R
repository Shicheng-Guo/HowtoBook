######################################################################################
######## DNA methylation based deep-learning for drug-response prediction ###########
######################################################################################
library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")
setwd("/home/guosa/hpc/project/TCGA")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/bin/id2phen4.R")

setwd("/mnt/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/meth450")
file=list.files(pattern="*mirnas.quantification.txt$",recursive = TRUE)
manifest2barcode("gdc_manifest.pancancer.miRNA.2019-05-29.txt")
barcode<-read.table("barcode.txt",sep="\t",head=T)
data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,sep="\t",as.is=F)  
  data<-cbind(data,tmp[,3])
  print(paste(i,"in",length(file),file[i],sep=" "))
  rownames(data)<-tmp[,1]
}

colnames(data)<-barcode[match(unlist(lapply(file,function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name),ncol(barcode)]
data<-data[,grep("TCGA",colnames(data))]
colnames(data)<-id2phen4(colnames(data))
data<-data[,match(unique(colnames(data)),colnames(data))]
miRNA<-data
save(miRNA,file="TCGA-Pancancer.miRNAseq.RData")

                                            
system("cp ~/hpc/methylation/Pancancer/methdata.pancancer.nomissing.RData ./")
load("methdata.pancancer.nomissing.RData")
colnames(input)<-id2phen4(colnames(input))
input<-input[,grep("-01",colnames(input))]

phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
phen$ID4<-paste(phen$bcr_patient_barcode,"-01",sep="")

input<-input[,colnames(input) %in% phen$ID4]
input<-input[,match(unique(colnames(input)),colnames(input))]
phen<-phen[na.omit(unlist(lapply(colnames(input),function(x) match(x,phen$ID)[1]))),]
dim(input)
dim(phen)

sort(table(phen$bcr_patient_barcode))
levels(phen$measure_of_response)
levels(phen$measure_of_response)<-c(0,1,1,0)

input<-data.frame(phen=phen$measure_of_response,t(input))
# input<-input[,unlist(apply(input,2,function(x) sd(x)>0))]
# library("SIS")
# x=data.matrix(input[,2:ncol(input)])
# y=as.numeric(input[,1])-1
# sisrlt<-SIS(x,y,family = c( "binomial"),penalty = c("lasso"))
# newx<-x[,sisrlt$ix]
# cpgs<-colnames(input)[sisrlt$ix]
# ncpgs<-data.frame(cpgs,cpg2symbol(cpgs))
# pdr<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/hsa01524.kegg.txt")
# write.table(ncpgs,file="pancancer.drugresponse.SIS.variables.txt",quote=F,row.names = F,col.names = F)

set.seed(49)
cv.error <- NULL
k <- 10
pbar <- create_progress_bar('text')
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  
  P=apply(train.cv[,2:ncol(train.cv)],2,function(x) summary(bayesglm(as.factor(train.cv[,1])~x,family=binomial))$coefficients[2,4])

  train.cv<-train.cv[,c(1,match(names(P[head(order(P),n=6000)]),colnames(train.cv)))]
  test.cv<-test.cv[,c(1,match(names(P[head(order(P),n=6000)]),colnames(test.cv)))]
  
  meth<-list()
  meth$train.cv=train.cv
  meth$test.cv=test.cv
  save(meth,file="meth.t6000.mxnet.RData")
  
  print(paste(ncol(train.cv),"variables passed P-value threshold and enrolled in SIS model"))
  RF <- randomForest(as.factor(phen) ~ ., data=train.cv, importance=TRUE,proximity=T)
  imp<-RF$importance
  imp<-imp[order(imp[,4],decreasing = T),]
  head(imp)
  write.table(imp,file=paste("RandomForest.VIP.Meth.",i,".txt",sep=""),sep="\t",quote=F,row.names = T,col.names = NA)
  topvar<-match(rownames(imp)[1:30],colnames(input))
  
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]
  
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,hidden=c(5,3),act.fct = "logistic",linear.output = F,threshold = 0.1)
  pr.nn <- neuralnet::compute(nn,test.cv)
  trainRlt<-data.frame(phen=train.cv[,1],pred=unlist(nn$net.result[[1]][,1]))
  testRlt<-data.frame(phen=test.cv[,1],pred=unlist(pr.nn$net.result[,1]))
  rownames(trainRlt)=row.names(train.cv)
  rownames(testRlt)=row.names(test.cv)
  rlt1<-rbind(rlt1,trainRlt)  
  rlt2<-rbind(rlt2,testRlt)
  print(i)
}
data1<-na.omit(data.frame(rlt1))
data2<-na.omit(data.frame(rlt2))
model.glm1 <- bayesglm(phen~.,data=rlt1,family=binomial(),na.action=na.omit)
model.glm2 <- bayesglm(phen~.,data=rlt2,family=binomial(),na.action=na.omit)
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))
