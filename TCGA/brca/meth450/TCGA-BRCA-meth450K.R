source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/meth450Pancancer.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/pancancer/methylation/meth450Pancancer.R")

id2phen4<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
}

id2phen3<-function(filename){
  library("stringr")
  as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*"))
}

id2bin<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
  as.numeric(lapply(strsplit(filename,"-"),function(x) x[4]))
}

filename=colnames(input)

id2pid<-function(filename){
  library("stringr")
  filename<-as.array(str_extract(filename,"edu_...."))
  unlist(lapply(filename,function(x) unlist(strsplit(x,"[_]"))[2]))
}

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*ncol(data)
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>=threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    output<-data[-NaRAW,]
  }else{
    output<-data;
  }
  output
}

cpg2symbol<-function(cpg){
  map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
  symbol<-map[match(cpg,map[,4]),5]
  return(symbol)
}

tsneplot<-function(mydata,phen,plot="tsne.plot.pdf"){
  library("tsne")
  data=data.frame(phen,mydata)
  pdf(plot)
  colors = rainbow(length(unique(data$phen)))
  names(colors) = unique(data$phen)
  ecb = function(x,y){plot(x,t='n'); text(x,labels=data$phen, col=colors[data$phen]) }
  tsne_iris = tsne(data, epoch_callback = ecb, perplexity=10)
  dev.off()
}
##################################################################################################### 
####################### Step 1: Read methylation 450K files ######################################### 
##################################################################################################### 
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
methdata[1:5,1:5]
phen4<-id2phen4(colnames(methdata))
phen3<-id2phen3(colnames(methdata))
bin<-id2bin(colnames(methdata))
pid<-id2pid(colnames(methdata))
##################################################################################################### 
####################### Step 2: Select Cancer Type ################################################## 
##################################################################################################### 
pid<-id2pid(colnames(methdata))
input<-methdata[,which(pid=="BRCA")]
phen4<-id2phen4(colnames(input))
phen3<-id2phen3(colnames(input))
bin<-id2bin(colnames(input))
pid<-id2pid(colnames(input))
input<-input[,c(which(bin==1),which(bin==11))]
phen4<-id2phen4(colnames(input))
phen3<-id2phen3(colnames(input))
bin<-id2bin(colnames(input))
pid<-id2pid(colnames(input))
phen<-data.frame(phen4=phen4,phen3=phen3,pid=pid,bin=bin)
head(phen)
##################################################################################################### 
####################### Step 3: Select BUR CpG probes ############################################### 
##################################################################################################### 
library("SIS")
library("arm")
library("randomForest")
map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
system("wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_mh450/Normal.PBMC.GEO.HM450K.Beta.RData")
load("Normal.PBMC.GEO.HM450K.Beta.RData")
system("wc -l ~/hpc/db/hg19/GPL13534_450K_hg19_PBMC_BUR.bed")
BUR<-read.table("~/hpc/db/hg19/GPL13534_450K_hg19_PBMC_BUR.bed")
input<-input[rownames(input) %in% BUR$V4,]
PDMR<-read.table("TCGA-Pancancer-MH450.Meta.diff.txt",head=T,row.names=1,sep="\t")
DMR<-subset(PDMR,beta>0.1 & pval<10^-5)
DMG<-na.omit(cpg2symbol(rownames(DMR)))
N<-length(unique(DMG))
input<-input[rownames(input)%in%rownames(DMR),]
bin<-id2bin(colnames(input))
input<-data.frame(phen=bin,t(input))
input<-data.frame(t(na.omit(t(input))))
input$phen[input$phen==11]<-0
save(input,file="TCGA-BRCA-BUR-PAN-BRCA-Beta.RData")

P=apply(input[,2:ncol(input)],2,function(x) summary(bayesglm(as.factor(input[,1])~x,family=binomial))$coefficients[2,4])
newinput<-input[,c(1,match(names(P[head(order(P),n=2000)]),colnames(input)))]
RF <- randomForest(as.factor(phen) ~ ., data=newinput, importance=TRUE,proximity=T)
imp<-RF$importance
imp<-imp[order(imp[,4],decreasing = T),]
imp<-imp[1:500,]
newimp<-data.frame(imp,Symbol=cpg2symbol(row.names(imp)))
newimp<-na.omit(newimp)
RLT<-data.frame(newimp,map[match(rownames(newimp),map[,4]),])
START=RLT$V2-150
END=RLT$V3+150
CHR=RLT$V1
RLT<-data.frame(CHR,START,END,RLT)
write.table(RLT,file="RF.BRCA.BUR.PAN.Top200VIP.txt",sep="\t",quote=F,col.names = NA,row.names = T)
##################################################################################################### 
####################### Step 4: Tumor Suppressor GENES ############################################### 
##################################################################################################### 
TSGTARGET<-newimp[newimp$Symbol %in% TSG[,2],]
GENE<-cpg2symbol(row.names(imp[1:500,]))
TSG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/TSGene2.0.txt",head=T,sep="\t")
GENE[GENE%in%TSG[,2]]
##################################################################################################### 
####################### Step 5: Tumor Suppressor GENES ############################################### 
##################################################################################################### 
set.seed(49)
cv.error <- NULL
k <- 10
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  P=apply(train.cv[,2:ncol(train.cv)],2,function(x) summary(bayesglm(as.factor(train.cv[,1])~x,family=binomial))$coefficients[2,4])
  train.cv<-train.cv[,c(1,match(names(P[head(order(P),n=200)]),colnames(train.cv)))]
  test.cv<-test.cv[,c(1,match(names(P[head(order(P),n=200)]),colnames(test.cv)))]
  meth<-list()
  meth$train.cv=train.cv
  meth$test.cv=test.cv
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
##################################################################################################### 
####################### Step 6: TSNE to show subtype ############################################### 
##################################################################################################### 
tsneplot(input[,match(rownames(RLT),colnames(input)),2:ncol(input)],input[,1],plot="TCGA.BRCA.BUR.TSNE.TopVariable.pdf")


