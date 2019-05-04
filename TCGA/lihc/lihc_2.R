if (!require("randomForest")) install.packages("randomForest")
if (!require("stringr")) install.packages("stringr")
if (!require("impute")) install.packages("impute")
if (!require("rpart")) install.packages("rpart")
if (!require("DAAG")) install.packages("DAAG")
if (!require("PredictABEL")) install.packages("PredictABEL")
if (!require("Deducer")) install.packages("Deducer")
if (!require("Rcmdr")) install.packages("Rcmdr")
if (!require("ROCR")) install.packages("ROCR")
if (!require("MASS")) install.packages("MASS")
if (!require("e1071")) install.packages("e1071")
if (!require("BayesTree")) install.packages("BayesTree")
if (!require("class")) install.packages("class")
if (!require("randomForest")) install.packages("randomForest")
if (!require("gbm")) install.packages("gbm")

library("BART")
library("Rtsne")  
library("ggplot2")  
library("readxl")
library("BayesTree")
library("randomForest")
library(tidyverse)
library(caret)
library(MASS)
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC")
source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/assess.R")
beta<-read.table("methyLevel_54samples.txt",head=T,sep="\t")
ROC<-read_xlsx("ROC.xlsx",sheet=1)

rownames(beta)<-paste(beta[,1],":",beta[,2],"-",beta[,3],sep="")
beta<-beta[,4:ncol(beta)]
phen<-read.table("sample_infomation.txt",sep="\t")
match(colnames(beta),phen$V1)

head(beta)
write.table(colMeans(beta,na.rm = T),file="colMeans.HCC.txt",sep="\t",quote=F)
# PCA
pca <- prcomp(t(na.omit(beta)),center=T,scale = F)
pdf("HCC.lowPass.PCA_SDEV.pdf")
par(cex.lab=1.5,cex.axis=1.5)
plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
dev.off()

scores <- data.frame(phen$V2, pca$x[,1:10])
pdf("HCC.MCRI.lowpass.PCA_1_2.pdf")
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),
     xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$V2))+1,cex=1.5)
phen$col=as.numeric(as.factor(phen$V2))+1
legend("topright",legend=c("HCC","HCC_Surgery","Normal"),pch=16,col=2:4,bty="n",cex=1.5)
dev.off()

pdf("HCC.MCRI.lowpass.Virus.PCA_1_2.pdf")
par(cex.lab=1.5,cex.axis=1.5)
plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),
     xlab="PC1",ylab="PC2",pch=16,col=as.numeric(as.factor(phen$V3))+1,cex=1.5)
phen$col=as.numeric(as.factor(phen$V3))+1
unique(data.frame(phen$V3,phen$col))
legend("bottom",legend=c("acute hepatitis","advanced HCC","after surgery","alcoholic cirrhosis","chronic hepatitis","cirrhosis","early stage HCC","health","nash-related cirrhosis"),pch=16,col=2:10,bty="n",cex=1)
dev.off()


source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/assess.R")
beta<-read.table("methyLevel_54samples.txt",head=T,sep="\t")
rownames(beta)<-paste(beta[,1],":",beta[,2],"-",beta[,3],sep="")
beta<-beta[,4:ncol(beta)]
phen<-read.table("sample_infomation.txt",sep="\t")
phen$V4=phen$V2

levels(phen$V4)=c("HCC","nonHCC","nonHCC")
input<-na.omit(data.frame(phen=phen$V4,t(beta)))

model<-randomForest(phen ~ ., data = input, ntree = 500, mtry = 5, importance = TRUE)
imp<-model$importance
vip<-imp[order(imp[,4],decreasing = T),]

case=which(phen$V4=="HCC")
con=which(phen$V4=="nonHCC")
P=apply(beta,1,function(x) wilcox.test(x[case],x[con],na.rm=T)$p.value)
newbeta<-data.matrix(beta[which(P<0.05/length(P)),])
boxplot(newbeta[1,]~phen$V4,col=2:3)

input<-data.frame(phen=phen$V4,t(na.omit(newbeta)))
glmRlt <- glm((as.numeric(as.factor(input$phen))-1) ~ . ,family=binomial, data = input)
Stepmodel <- stepAIC(glmRlt,trace =F,direction="both",steps=2000)

Stepmodel$anova
summary(Stepmodel)

beta<-data.matrix(beta)
P<-c()
for(i in 1:nrow(beta)){
  glmRlt <- summary(glm((as.numeric(as.factor(input$phen))-1) ~ beta[i,],family=binomial))
  P<-c(P,glmRlt$coefficients[2,4])
}

newbeta<-beta[which(P<0.002),]
newinput<-data.frame(phen=abs(as.numeric(as.factor(input$phen))-2),t(na.omit(beta)))
glmRlt <-glm(phen~.,family=binomial,data=newinput)
StepModel<-stepwise(glmRlt)
summary(StepModel)


# Rtsne
rlt<-Rtsne(t(na.omit(beta)))
input<-data.frame(x = rlt$Y[,1], y = rlt$Y[,2], col = as.numeric(colnames(methdata)) %% 2 +2, ID=colnames(methdata))
plot(input[,1:2],pch=16,cex=1.5,col=input$col,xlab="Dimention 1",ylab="Dimention 2")
text(x=input$x-0.25,y=input$y,labels=colnames(methdata),cex=0.85)
legend("bottomright",legend=c("CHOL","NORMAL"),pch=16,col=c(2,3),bty="n",cex=1.5)

# Five-fold cross-validation
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC")
beta<-read.table("methyLevel_54samples.txt",head=T,sep="\t")
rownames(beta)<-paste(beta[,1],":",beta[,2],"-",beta[,3],sep="")
beta<-beta[,4:ncol(beta)]
phen<-read.table("sample_infomation.txt",sep="\t")
match(colnames(beta),phen$V1)
phen$V4=phen$V2
levels(phen$V4)=c("1","0","0")
input<-data.frame(t(na.omit(beta)))
phen<-phen$V4

N=10
K=5  
rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-matrix(NA,N*K,6)
j<-1
for(n in 1:N){
  CV<-CvSampling(nrow(input),5)
  for(i in 1:K){
    print(j)
    train=input[CV[[i]]$train,]
    trainlabel<-as.factor(phen[CV[[i]]$train])
    test=input[CV[[i]]$test,]
    testlabel<-phen[CV[[i]]$test]
    ### logistic model
    
    model<-glm(trainlabel~.,family=binomial,data=train)
    predict(model,train,type="response")
    predict(model,test,type="response")
    predict(model, train)
    x<-predict(model, train)
    rlt1[j,]<-assess(model,train,test,trainlabel,testlabel)
    ## 2. KNN model
    model1<-knn(train=input[CV[[i]]$train,],test=input[CV[[i]]$train,],cl=phen[CV[[i]]$train],k=6)
    model2<-knn(train=input[CV[[i]]$train,],test=input[CV[[i]]$test,],cl=phen[CV[[i]]$train],k=6)
    rlt2[j,]<-assess2(model1,model2,trainlabel,testlabel)    
    ## 3. N*K-fold cross-validation for SVM classification    
    model<-svm(x=train,y=trainlabel)
    rlt3[j,]<-assess3(model,train,test,trainlabel,testlabel)    
    ## 4. N*K-fold cross-validation for BayesTree
    #model1 = (apply(pnorm(bart(train,trainlabel,train)$yhat.test),2,mean))
    #model2 = (apply(pnorm(bart(train,trainlabel,test)$yhat.test),2,mean))
    #rlt4[j,]<-assess4(model1,model2,trainlabel,testlabel)    
    ## 5. N*K-fold cross-validation for randomForest
    model <- randomForest(trainlabel~., data=train)
    rlt5[j,]<-assess3(model,train,test,trainlabel,testlabel)    
    ## 6. N*K-fold cross-validation for Boosting
    model <- gbm(trainlabel~.,data=train,distribution = "adaboost")
    model1 <- predict.gbm(model,train,100,type="response")
    model2 <- predict.gbm(model,test,100)
    rlt6[j,]<-assess2(model1,model2,trainlabel,testlabel)    
    j=j+1
  }
}


colnames(rlt5)<-c("Spe.train","Sen.train", "Accu.train","Spe.test","Sen.test","Acc.test")
write.table(na.omit(rlt5),file="RandomForest_5_fold.txt",sep="\t",quote=F,col.names = NA,row.names = T)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC")
source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/assess.R")
data<-read_xlsx("ROC.xlsx",sheet=1)
data<-na.omit(data.frame(data))
dim(data)

model.glm1 <- (bayesglm(Phen ~ Sex+Age,family=binomial,data=data,na.action=na.omit))
model.glm2 <- (bayesglm(Phen ~ Sex+Age+ALT+AST+Tbil,family=binomial,data=data,na.action=na.omit))
model.glm3 <- (bayesglm(Phen ~ Sex+Age+AFP,family=binomial,data=data,na.action=na.omit))
model.glm4 <- (bayesglm(Phen ~ Sex+Age+WGBS,family=binomial,data=data,na.action=na.omit))
model.glm5 <- (bayesglm(Phen ~ Sex+Age+AFP+WGBS,family=binomial,data=data,na.action=na.omit))

library("PredictABEL")
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
pred3 <- predRisk(model.glm3)
pred4 <- predRisk(model.glm4)
pred5 <- predRisk(model.glm5)
length(pred1)
length(pred2)
length(pred3)
length(pred4)
length(pred5)

par(cex.lab=1.5,cex.axis=1.5,lwd=1.5)
plotROC(data=data,cOutcome=3,predrisk=cbind(pred1,pred2,pred3,pred4,pred5))
legend("bottomright",cex=1.5,legend=c("Sex+Age","Sex+Age+ALT+AST+Tbil","Sex+Age+AFP","Sex+Age+WGBS","Sex+Age+AFP+WGBS"),lty=1:5,col=1:5,bty="n")




model.glm1 <- (bayesglm(Phen ~ AFP,family=binomial,data=data,na.action=na.omit))
model.glm2 <- (bayesglm(Phen ~ WGBS,family=binomial,data=data,na.action=na.omit))
model.glm3 <- (bayesglm(Phen ~ AFP+WGBS,family=binomial,data=data,na.action=na.omit))

library("PredictABEL")
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
pred3 <- predRisk(model.glm3)
length(pred1)
length(pred2)
length(pred3)

par(cex.lab=1.5,cex.axis=1.5,lwd=1.5)
plotROC(data=data,cOutcome=3,predrisk=cbind(pred1,pred2,pred3))
legend("bottomright",cex=1.5,legend=c("AFP","WGBS","AFP+WGBS"),lty=1:5,col=1:5,bty="n")
