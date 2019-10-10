# CGH1
install.packages("caret")
library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")
library("caret")

setwd("/home/guosa/hpc/project/TCGA")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")

phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
barcode<-read.table("/home/guosa/hpc/project/TCGA/pancancer/FPKM/barcode.txt",head=T,sep="\t")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/bin/id2phen4.R")

load("rnaseqdata.pancancer.env.RData")

ncn<-barcode[match(unlist(lapply(colnames(rnaseqdata),function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name<-gsub(".gz","",barcode$file_name)),]
ncol<-match("cases.0.samples.0.submitter_id",colnames(ncn))
colnames(rnaseqdata)<-ncn[,ncol]
phen$ID<-paste(phen$bcr_patient_barcode,"-01",sep="")
rnaseq<-rnaseqdata[,na.omit(match(unique(phen$ID),id2phen4(colnames(rnaseqdata))))]
rnaseq<-rnaseq[which(unlist(apply(rnaseq,1,function(x) sd(x)>0))),]
colnames(rnaseq)<-id2phen4(colnames(rnaseq))
newphen<-phen[unlist(lapply(colnames(rnaseq),function(x) match(x,phen$ID)[1])),]
sort(table(newphen$bcr_patient_barcode))
levels(newphen$measure_of_response)<-c(0,1,1,0)
newinput<-data.frame(phen=newphen$measure_of_response,t(rnaseq))

######################################################################################
P=apply(newinput[,2:ncol(newinput)],2,function(x) summary(glm(as.factor(newinput[,1])~x,family=binomial))$coefficients[2,4])
pQQ(P, nlabs =length(P), conf = 0.95, mark = F) 
newinput<-newinput[,c(1,which(P<8.5*10^-7)+1)]
save(newinput,file="Pancancer.DrugResponse.V5292.N1462.RData")
load("Pancancer.DrugResponse.V5292.N1462.RData")
######################################################################################
library("SIS")
x=data.matrix(newinput[,2:ncol(newinput)])
y=as.numeric(newinput[,1])-1
sisrlt<-SIS(x,y,family = c( "binomial"),penalty = c("lasso"))
newx<-x[,sisrlt$ix]

fit<-bayesglm(phen~.,data=input,family=binomial(),na.action=na.omit)
pred <- predRisk(fit)
plotROC(data=input,cOutcome=1,predrisk=cbind(pred))

pdf("ROC.pdf")
par(cex.lab=1.5,cex.axis=1.5)
plotROC(data=input,cOutcome=1,predrisk=cbind(pred))
dev.off()

input<-data.frame(phen=y,newx)
ENSG<-unlist(lapply(colnames(newx),function(x) unlist(strsplit(x,"[.]"))[1]))
ENSG2Symbol(ENSG)
RF <- randomForest(as.factor(phen) ~ ., data=input,mtry=20,importance=TRUE,proximity=T)
RF$importance<-RF$importance[order(RF$importance[,4],decreasing = T),]
Symbol=ENSG2Symbol(unlist(lapply(rownames(RF$importance),function(x) unlist(strsplit(x,"[.]"))[1])))
RLT<-data.frame(RF$importance,Symbol)
write.table(RLT,file="SIS.RF.VIP.txt",sep="\t",quote=F,row.names = T,col.names = NA)
######################################################################################
input<-newinput
RF <- randomForest(as.factor(phen) ~ ., data=input,mtry=10,importance=TRUE,proximity=T)
RF$importance<-RF$importance[order(RF$importance[,4],decreasing = T),]
Symbol=ENSG2Symbol(unlist(lapply(rownames(RF$importance),function(x) unlist(strsplit(x,"[.]"))[1])))
RLT<-data.frame(RF$importance,Symbol)
write.table(RLT,file="RF.VIP.txt",sep="\t",quote=F,row.names = T,col.names = NA)
######################################################################################
#########  mRNA-seq based deep-learning for drug-response prediction #################
######################################################################################
library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")

setwd("/home/guosa/hpc/project/TCGA/pancancer/FPKM")

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
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/FPKM")
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
######################################################################################
############  miRNA based deep-learning for drug-response prediction #################
######################################################################################
library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")

setwd("/home/guosa/hpc/project/TCGA/pancancer/miRNA/data")
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
colnames(data)<-id2phen4(barcode[match(unlist(lapply(file,function(x) unlist(strsplit(x,"[/]"))[2])),barcode$file_name),ncol(barcode)])
data<-data[,match(unique(colnames(data)),colnames(data))]
save(data,file="TCGA-Pancancer.miRNAseq.RData")
miRNA<-data[,grep("-01",colnames(data))]

load("TCGA-Pancancer.miRNAseq.RData")

barcode$id4=id2phen4(barcode$cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id)
phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
phen$ID4<-paste(phen$bcr_patient_barcode,"-01",sep="")

miRNA<-miRNA[,colnames(miRNA) %in% phen$ID4]
phen<-phen[na.omit(unlist(lapply(colnames(miRNA),function(x) match(x,phen$ID)[1]))),]
dim(miRNA)
dim(phen)

sort(table(phen$bcr_patient_barcode))
table(levels(phen$measure_of_response))
levels(phen$measure_of_response)<-c(0,1,1,0)

input<-data.frame(phen=phen$measure_of_response,t(miRNA))
input<-input[,unlist(apply(input,2,function(x) sd(x)>0))]
miRNA<-input
save(miRNA,file="pancancer.miRNA.drugResponse.RData")

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
  train.cv<-train.cv[,c(1,which(P<0.05/length(P))+1)]
  
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
  
  nn <- neuralnet(f,data=train.cv,hidden=c(10,3),act.fct = "logistic",linear.output = F)
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
          
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
newinput<-t(input[,match(rownames(imp)[1:50],colnames(input))])
colnames(newinput)<-input[,1]
pdf("heatmap.randomForest.pdf")
HeatMap(newinput)
dev.off()
save.image("RNAseq-N2.RF.heatmap.RData")



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


map<-read.table("GPL13534_450K_hg19_V2.bed",sep="\t")
map<-map[,c(1:5,7,8)]
write.table(map,file="GPL13534_450K_hg19_V3.bed",sep="\t",quote=F,col.names = F,row.names = F)

cpg2symbol<-function(cpg){
map<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/GPL13534_450K_hg19_V3.bed")
symbol<-map[match(cpg,map[,4]),5]
return(symbol)
}


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/FPKM")
load("mRNA.t5200.mxnet.RData")
train.cv=data.matrix(meth$train.cv)
test.cv=data.matrix(meth$test.cv)

mx.set.seed(0)
model <- mx.mlp(train.cv[,2:ncol(train.cv)],train.cv[,1]-1, hidden_node=c(200,4), out_node=2, out_activation="softmax",num.round=20, array.batch.size=15, learning.rate=0.07, momentum=0.9,eval.metric=mx.metric.accuracy)
graph.viz(model$symbol)
preds = predict(model, test.cv[,2:ncol(test.cv)])
pred.label = max.col(t(preds))-1
Tab<-table(pred.label, test.cv[,1]-1)
Tab
acc=(Tab[1,1]+Tab[2,2])/sum(Tab)
acc


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/meth450")
require(mxnet)
require(mlbench)
load("meth.t6000.mxnet.RData")

load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/FPKM/mRNA.t3000.mxnet.RData")
mRNA<-rbind(mRNA$train.cv,mRNA$test.cv)
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/meth450/meth.t6000.mxnet.RData")
meth<-rbind(meth$train.cv,meth$test.cv)
load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/miRNA/pancancer.miRNA.drugResponse.RData")

triple<-names(which(sort(table(c(rownames(miRNA),rownames(meth),rownames(mRNA))))==3))

miRNA<-miRNA[match(triple,rownames(miRNA)),]
mRNA<-mRNA[match(triple,rownames(mRNA)),]
meth<-meth[match(triple,rownames(meth)),]

data<-cbind(miRNA,mRNA,meth)
input<-cbind(miRNA[,1],mRNA[,2:ncol(mRNA)],meth[,2:ncol(meth)])
input<-input[,c(1,sample(1:ncol(input),2500))]
colnames(input)[1]<-"phen"
input[1:5,1:5]

index <- sample(1:nrow(input),round(0.9*nrow(input)))
train.cv <- data.matrix(input[index,])
test.cv <-  data.matrix(input[-index,])

apply(input,2,function(x) summary(glm(as.numeric(input$phen)~input[,2]))$coefficient[2,4])

source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/PCAPlot.R")
PCAPlot(input[,2:ncol(input)],as.numeric(input[,1]),"drugResponse.pca.pdf")

library("tsne")
data=input
colors = rainbow(length(unique(data$phen)))
names(colors) = unique(data$phen)
ecb = function(x,y){plot(x,t='n'); text(x,labels=data$phen, col=colors[data$phen]) }
tsne_iris = tsne(data.matrix(data), epoch_callback = ecb, perplexity=30)

pca <- prcomp(input$data,center=T,scale = F)

input$data=input[,2:ncol(input)]
input$phen=as.numeric(input[,1])
pca <- prcomp(input$data,center=T,scale = F)  # Here, input file: row is individual and column is variable
scores <- data.frame(input$phen, pca$x[,1:3])
par(cex=1)
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
col = as.numeric(as.factor(input$phen))
for(i in 1:length(scores$PC1)){
  points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(input$phen))[i],col=col[i]+1,cex=0.8,lwd=2)
}
legend("topright",legend=names(table(input$phen)),pch=1:length(table(input$phen)),col=1:length(table(input$phen))+1,bty="n",cex=2)


library(randomForest)
input[1:5,1:5]
input<-data.matrix(input)
RF<-randomForest(phen~.,data=input)


library("mxnet")
mx.set.seed(0)
model <- mx.mlp(train.cv[,2:ncol(train.cv)],train.cv[,1]-1, hidden_node=c(20,5), out_node=2, num.round=20,out_activation="softmax", array.batch.size=15, learning.rate=0.07, momentum=0.9,eval.metric=mx.metric.accuracy)
graph.viz(model$symbol)
preds = predict(model, test.cv[,2:ncol(test.cv)])
pred.label = max.col(t(preds))-1
pred.label
Tab<-table(pred.label, test.cv[,1]-1)
Tab
acc=(Tab[1,1]+Tab[2,2])/sum(Tab)
acc



load("TCGA-Pancancer.mRNAseq.RData")
phen<-id2bin(colnames(data))
data<-data[,which(phen==1 | phen==11)]
phen<-id2bin(colnames(data))
input=data.frame(phen=phen,t(data))
input<-input[,unlist(apply(input,2,function(x) sd(x)>0))]
index <- sample(1:nrow(input),round(0.9*nrow(input)))
train.cv <- input[index,]
test.cv <- input[-index,]
mRNA<-list()
mRNA$train.cv<-train.cv
mRNA$test.cv<-test.cv

input<-input[,c(1,sample(2:ncol(input),2000))]
save(input,file="TCGA-Pancancer.mRNAseq.CancerNormal.subset.RData")


load("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/TCGA/pancancer/FPKM/data/TCGA-Pancancer.mRNAseq.CancerNormal.subset.RData")
input$phen[input$phen==11]<-0
index <- sample(1:nrow(input),round(0.9*nrow(input)))
train.cv <- data.matrix(input[index,])
test.cv <- data.matrix(input[-index,])
mx.set.seed(0)
model <- mx.mlp(train.cv[,2:ncol(train.cv)],train.cv[,1], hidden_node=c(30,4), out_node=2, num.round=10,out_activation="softmax", array.batch.size=15, learning.rate=0.07, momentum=0.9,eval.metric=mx.metric.accuracy)
graph.viz(model$symbol)
preds = predict(model, test.cv[,2:ncol(test.cv)])
pred.label = max.col(t(preds))
table(pred.label)
Tab<-table(pred.label, test.cv[,1]-1)
Tab
acc=(Tab[1,1]+Tab[2,2])/sum(Tab)
acc




