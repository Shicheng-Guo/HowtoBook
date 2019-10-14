# CHG1

library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")

Symbol2ENSG<-function(Symbol){
  db<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed",sep="\t")
  ENSG<-na.omit(as.character(db[match(Symbol,db$V4),8]))
  return(ENSG)
}

load("~/hpc/project/TCGA/TCGA-Pancancer.mRNAseq.trim.RData")
vip<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/pharmacogenomics/RA/VIP.Gene.hg19.bed",head=F)
ensg<-Symbol2ENSG(vip$V4)
grep(ensg,colnames(input))
xsel<-unlist(lapply(ensg,function(x) grep(x,colnames(mRNA))))

input<-data.matrix(mRNA[,xsel])
load("triple.id.RData")
load("drugResponse.phen.RData")
input<-input[match(triple,rownames(input)),]
input<-data.frame(phen,log(input+1,2))
input[1:5,1:5]

colnames(input)<-newinput[,1]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
pdf("mRNA.vip.heatmap.n2.pdf")
HeatMap(input)
dev.off()

save(input,file="vip.mRNA.drugresponse.RData")
load("vip.mRNA.drugresponse.RData")
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
