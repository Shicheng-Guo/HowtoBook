######################################################################################
#########  mRNA-seq based deep-learning for drug-response prediction #################
######################################################################################
install.packages("caret")

library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")
library("caret")

setwd("/home/guosa/hpc/project/TCGA")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")

setwd("~/hpc/project/TCGA/pancancer/FPKM")
file=list.files(pattern="*FPKM-UQ.txt$",recursive = TRUE)
manifest2barcode("gdc_manifest.pancancer.FPKM.2019-05-29.txt")
barcode<-read.table("barcode.txt",sep="\t",head=T)
data<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=F,sep="\t",as.is=F)  
  data<-cbind(data,tmp[,2])
  print(paste(i,"in",length(file),file[i],sep=" "))
  rownames(data)<-tmp[,1]
}

barcode$file_name<-gsub(".gz","",barcode$file_name)
colnames(data)<-id2phen4(barcode[match(unlist(lapply(file,function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name),]$cases.0.samples.0.submitter_id)
data<-data[,grep("TCGA",colnames(data))]
data<-data[,match(unique(colnames(data)),colnames(data))]
save(data,file="TCGA-Pancancer.mRNAseq.RData")

phen<-id2bin(colnames(data))
data<-data[,which(phen==1 | phen==11)]
phen<-id2bin(colnames(data))
input=data.frame(phen=phen,t(data))
input<-input[,unlist(apply(input,2,function(x) sd(x)>0))]
index <- sample(1:nrow(input),round(0.9*nrow(input)))
train.cv <- input[index,]
test.cv <- input[-index,]
mRNA<-data[,grep("-01",colnames(data))]

setwd("~/hpc/project/TCGA/pancancer/FPKM")
load("Pancancer.DrugResponse.V5292.N1462.RData")
input<-newinput

set.seed(49)
cv.error <- NULL
k <- 10
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  
  P=apply(input[,2:ncol(input)],2,function(x) summary(glm(as.factor(input[,1])~x,family=binomial))$coefficients[2,4])
  
  train.cv<-train.cv[,c(1,match(names(P[head(order(P),n=5200)]),colnames(train.cv)))]
  test.cv<-test.cv[,c(1,match(names(P[head(order(P),n=5200)]),colnames(test.cv)))]
  
  meth<-list()
  meth$train.cv=train.cv
  meth$test.cv=test.cv
  save(meth,file="mRNA.t5200.mxnet.RData")
  
  RF <- randomForest(as.factor(phen) ~ ., data=train.cv, importance=TRUE,proximity=T)
  imp<-RF$importance
  head(imp)
  imp<-imp[order(imp[,4],decreasing = T),]
  write.table(imp,file=paste("RandomForest.VIP.",i,".txt",sep=""),sep="\t",quote=F,row.names = T,col.names = NA)
  topvar<-match(rownames(imp)[1:30],colnames(input))
  
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]
  
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,hidden=c(10,3),act.fct = "logistic",linear.output = T)
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

### ROC 
pdf("mRNA.drugresponse.pdf")
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))
dev.off()

### heatmap
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
setwd("~/hpc/project/TCGA/pancancer/FPKM")
load("Pancancer.DrugResponse.V5292.N1462.RData")
input<-newinput
RF <- randomForest(as.factor(phen) ~ ., data=input, importance=TRUE,proximity=T)
imp<-RF$importance
head(imp)
imp<-imp[order(imp[,4],decreasing = T),]

newinput<-t(log(input[,match(rownames(imp)[1:50],colnames(input))]+1,2))
colnames(newinput)<-input[,1]
newinput[1:5,1:5]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
pdf("mRNA.heatmap.randomForest.n2.pdf")
HeatMap(newinput)
dev.off()
save.image("miRNAseq-N2.RF.heatmap.RData")
