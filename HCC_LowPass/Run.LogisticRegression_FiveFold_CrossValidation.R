# Five-fold cross-validation
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC")
source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/assess.R")


assess<-function(model,train,test,trainlabel,testlabel){
  t1<-table(chan(predict(model, train)), trainlabel)
  t2<-table(chan(predict(model, test)), testlabel)
  if(nrow(t1)==2 && nrow(t2)==2){
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    rlt
  }
}

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
    rlt<-apply(train,2,function(x) summary(bayesglm(trainlabel~x,family=binomial))$coefficients[2,4])
    train=train[,head(order(rlt,decreasing = F),n=nrow(input))]
    test=test[,head(order(rlt,decreasing = F),n=nrow(input))]
    model<-stepwise(bayesglm(trainlabel~.,family=binomial,data=train),criterion="AIC")
    predict(model,train,type="response")
    predict(model,test,type="response")
    predict(model, train)
    x<-predict(model, train)
    temp<-assess(model,train,test,trainlabel,testlabel)
    if(! is.null(temp)){
      rlt1[j,]<-temp
    }
    j=j+1
  }
}
colnames(rlt1)<-c("Spe.train","Sen.train", "Accu.train","Spe.test","Sen.test","Acc.test")
rlt1<-na.omit(rlt1)
write.table(na.omit(rlt2),file="HCC_LowPass_LogisticRegression_FiveFold.txt",sep="\t",quote=F,col.names = NA,row.names = T)
save.image("LogisticRegression.LowPass.RData")

