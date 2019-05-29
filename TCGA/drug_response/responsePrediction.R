library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")

setwd("/home/guosa/hpc/project/TCGA")

phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
barcode<-read.table("/home/guosa/hpc/project/TCGA/pancancer/FPKM/barcode.txt",head=T,sep="\t")
load("rnaseqdata.pancancer.env.RData")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/bin/id2phen4.R")

ncn<-barcode[match(unlist(lapply(colnames(rnaseqdata),function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name<-gsub(".gz","",barcode$file_name)),]
ncol<-match("cases.0.samples.0.submitter_id",colnames(ncn))
colnames(rnaseqdata)<-ncn[,ncol]
phen$ID<-paste(phen$bcr_patient_barcode,"-01",sep="")

rnaseq<-rnaseqdata[,na.omit(match(unique(phen$ID),id2phen4(colnames(rnaseqdata))))]
colnames(rnaseq)<-id2phen4(colnames(rnaseq))
newphen<-phen[unlist(lapply(colnames(rnaseq),function(x) match(x,phen$ID)[1])),]
sort(table(newphen$bcr_patient_barcode))
table(newphen$measure_of_response)
levels(newphen$measure_of_response)<-c(0,1,0,0)
input<-data.frame(phen=newphen$measure_of_response,t(rnaseq))
input[1:10,1:10]

SD<-unlist(apply(input[,2:ncol(input)],2,function(x) sd(x)))

set.seed(49)
cv.error <- NULL
k <- 10
pbar <- create_progress_bar('text')
pbar$init(k)
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  
  RF <- randomForest(as.factor(phen) ~ ., data=train.cv, importance=TRUE,proximity=F)
  imp<-RF$importance
  head(imp)
  imp<-imp[order(imp[,4],decreasing = T),]
  topvar<-match(rownames(imp)[1:10],colnames(input))
  
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]
  
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,hidden=c(3),act.fct = "logistic",linear.output = FALSE)
  #plot(nn,lwd=0.85,cex=1)
  pr.nn <- neuralnet::compute(nn,test.cv)
  trainRlt<-data.frame(phen=train.cv[,1],pred=unlist(nn$net.result))
  testRlt<-data.frame(phen=test.cv[,1],pred=unlist(pr.nn$net.result))
  rownames(trainRlt)=row.names(train.cv)
  rownames(testRlt)=row.names(test.cv)
  rlt1<-rbind(rlt1,trainRlt)  
  rlt2<-rbind(rlt2,testRlt)
  pbar$step()
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
