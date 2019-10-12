load("~/hpc/project/TCGA/pancancer/FPKM/Pancancer.DrugResponse.V5292.N1462.RData")

load("~/hpc/project/TCGA/pancancer/FPKM/TCGA-Pancancer.mRNAseq.RData")
load("~/hpc/project/TCGA/methdata.pancancer.trim.RData")
load("~/hpc/project/TCGA/pancancer/miRNA/pancancer.miRNA.drugResponse.RData")

triple<-names(which(sort(table(c(rownames(miRNA),rownames(meth),rownames(mRNA))))==3))

miRNA<-miRNA[match(triple,rownames(miRNA)),]
mRNA<-mRNA[match(triple,rownames(mRNA)),]
meth<-meth[match(triple,rownames(meth)),]

input<-cbind(phen=miRNA[,1],miRNA[,2:ncol(miRNA)],mRNA[,2:ncol(mRNA)],meth[,2:ncol(meth)])
SD<-apply(input,2,function(x) sd(x))
input<-input[,-which(SD==0)]
input[,2:ncol(input)]<-scale(input[,2:ncol(input)])
input[1:5,1:5]
#######################################################
# index <- sample(1:nrow(input),round(0.9*nrow(input)))
# train.cv <- data.matrix(input[index,])
# test.cv <-  data.matrix(input[-index,])
# library("mxnet")
# mx.set.seed(0)
# model <- mx.mlp(train.cv[,2:ncol(train.cv)],train.cv[,1]-1, hidden_node=c(20,5), out_node=2, num.round=20,out_activation="softmax", array.batch.size=15, learning.rate=0.07, momentum=0.9,eval.metric=mx.metric.accuracy)
# graph.viz(model$symbol)
# preds = predict(model, test.cv[,2:ncol(test.cv)])
# pred.label = max.col(t(preds))-1
# pred.label
# Tab<-table(pred.label, test.cv[,1]-1)
# Tab
# acc=(Tab[1,1]+Tab[2,2])/sum(Tab)
# acc
#######################################################
set.seed(49)
cv.error <- NULL
k <- 2
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  print(i)
  index <- sample(1:nrow(input),round(0.5*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  
  P=apply(train.cv[,2:ncol(train.cv)],2,function(x) summary(glm(as.factor(train.cv[,1])~x,family=binomial))$coefficients[2,4])

  train.cv<-train.cv[,c(1,match(names(P[head(order(P),n=sum(P<0.05/length(P)))]),colnames(train.cv)))]
  test.cv<-test.cv[,c(1,match(names(P[head(order(P),n=sum(P<0.05/length(P)))]),colnames(test.cv)))]
  
  print(paste(ncol(train.cv),"variables passed P-value threshold and enrolled in SIS model"))
  RF <- randomForest(as.factor(phen) ~ ., data=train.cv, importance=TRUE,proximity=T)
  imp<-RF$importance
  imp<-imp[order(imp[,4],decreasing = T),]
  head(imp)
  write.table(imp,file=paste("RandomForest.VIP.Meth.",i,".txt",sep=""),sep="\t",quote=F,row.names = T,col.names = NA)
  topvar<-match(rownames(imp)[1:min(nrow(input)/2,nrow(imp)/2)],colnames(input))
  
  train.cv <- input[index,c(1,topvar)]
  print(train.cv[1:5,1:5])
  test.cv <- input[-index,c(1,topvar)]
  print(test.cv[1:5,1:5])
          
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,stepmax=5*10^5,hidden=c(10,3),act.fct = "logistic",linear.output = F,threshold = 0.1)
  pr.nn <- neuralnet::compute(nn,test.cv)
  trainRlt<-data.frame(phen=train.cv[,1],pred=unlist(nn$net.result[[1]][,1]))
  testRlt<-data.frame(phen=test.cv[,1],pred=unlist(pr.nn$net.result[,1]))
  rownames(trainRlt)=row.names(train.cv)
  rownames(testRlt)=row.names(test.cv)
  print(head(rownames(testRlt)))
  rlt1<-rbind(rlt1,trainRlt)  
  rlt2<-rbind(rlt2,testRlt)
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
pdf("meth.ROC.drugresponse.pdf")
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))
dev.off()

### heatmap
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
setwd("~/hpc/project/TCGA/pancancer/meth450")
input<-data.frame(input)
P=apply(input[,2:ncol(input)],2,function(x) summary(glm(as.factor(input[,1])~x,family=binomial))$coefficients[2,4])

RF <- randomForest(as.factor(phen) ~ ., data=input, importance=TRUE,proximity=T)
imp<-RF$importance
head(imp)
imp<-imp[order(imp[,4],decreasing = T),]
topvar<-match(rownames(imp)[1:200],colnames(input))
newinput <- t(input[,topvar])
colnames(newinput)<-newinput[,1]
newinput[1:5,1:5]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
pdf("triple.heatmap.randomForest.n2.pdf")
HeatMap(newinput)
dev.off()
        

               
