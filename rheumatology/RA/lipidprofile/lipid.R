# this code is for lipid profile analysis in RA patients. 

#######################################################################################
######################### load relevant librarsy  #########################
#######################################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("randomForest")
# biocLite("stringr")
# biocLite("impute")
# biocLite("rpart")

library("stringr")
library("randomForest")
library("impute")
library("rpart")


#######################################################################################
######################### function definition and description #########################
#######################################################################################
calcOddsRatio <- function(mymatrix,alpha=0.05,referencerow=2){
  numrow <- nrow(mymatrix)
  myrownames <- rownames(mymatrix)
  for (ii in 1:numrow){
    rowname <- myrownames[ii]
    DiseaseUnexposed <- mymatrix[referencerow,1]
    ControlUnexposed <- mymatrix[referencerow,2]
    if (ii != referencerow){ 
      DiseaseExposed <- mymatrix[ii,1]
      ControlExposed <- mymatrix[ii,2] 
      totExposed <- DiseaseExposed + ControlExposed
      totUnexposed <- DiseaseUnexposed + ControlUnexposed
      probDiseaseGivenExposed <- DiseaseExposed/totExposed
      probDiseaseGivenUnexposed <- DiseaseUnexposed/totUnexposed
      probControlGivenExposed <- ControlExposed/totExposed
      probControlGivenUnexposed <- ControlUnexposed/totUnexposed
      # calculate the odds ratio
      oddsRatio <- (probDiseaseGivenExposed*probControlGivenUnexposed)/(probControlGivenExposed*probDiseaseGivenUnexposed)
      # calculate a confidence interval
      confidenceLevel <- (1 - alpha)*100
      sigma <- sqrt((1/DiseaseExposed)+(1/ControlExposed)+(1/DiseaseUnexposed)+(1/ControlUnexposed))
      # sigma is the standard error of our estimate of the log of the odds ratio
      z <- qnorm(1-(alpha/2))
      lowervalue <- oddsRatio * exp(-z * sigma)
      uppervalue <- oddsRatio * exp( z * sigma)
      OR=paste(round(oddsRatio,2)," (",round(lowervalue,2),",",round(uppervalue,2),")",sep="")
    }
  }
  OR
}
matrixanova<-function(x){
  # x is var matrix, 0 or 1, + or - 
  # remmber not to include sample id
  rlt<-list()
  data<-data.frame(x)
  var<-unlist((lapply(data,function(x) length(names(table(x))))))
  catvar<-as.numeric(unlist(which(lapply(data,function(x) length(names(table(x))))<5)))
  nullvar1<-unlist(lapply(data[,catvar],function(x) sum(names(table(x))=="")))
  nullvar2<-var[catvar]
  nullvar2-nullvar1
  contvar<-as.numeric(unlist(which(lapply(data,function(x) length(names(table(x))))>=5)))
  pvalue<-matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
  table<-matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
  OR<-matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
  for(i in catvar){
    for (j in c(catvar,contvar)){
      if(j !=i ){
        tmp<-which(apply(cbind(data[,i],data[,j]),1,function(x) all(!is.na(x))))
        if(j %in% catvar & any(table(data[tmp,1],data[tmp,j])>5)){
          pvalue[i,j]<-chisq.test(data[tmp,i],data[tmp,j],simulate.p.value = TRUE)$p.value
          tmp<-table(data[tmp,i],data[tmp,j])
          table[i,j]<-paste(tmp[1,1],"/",tmp[1,2],":",tmp[2,1],"/",tmp[2,2],sep="")
          OR[i,j]<-calcOddsRatio(tmp,alpha=0.05,referencerow=1)
        }else if(j %in% contvar){
          pvalue[i,j]<-oneway.test(as.numeric(data[tmp,j]) ~ data[tmp,i],na.action="na.omit")$p.value
        }
      }
    }
  }
  colnames(pvalue)=colnames(data)
  rownames(pvalue)=colnames(data)
  colnames(OR)=colnames(data)
  rownames(OR)=colnames(data)
  colnames(table)=colnames(data)
  rownames(table)=colnames(data)
  rlt$pvalue=pvalue
  rlt$table=table
  rlt$OR=OR
  return(rlt)
}
OR95CIPV<-function(dat){
  out<-array()
  for(i in 1:dim(dat)[2]){
    glm<-lapply(dat,function(x) glm(abs(as.numeric(phen)-2)~gender+age+x,data,family=binomial(logit)))
    pvalue=summary(glm[[i]])$coefficients[4,4]
    beta=summary(glm[[i]])$coefficients[4,1]
    se=summary(glm[[i]])$coefficients[4,2]
    fit<-round(exp(beta),2)
    up<-round(exp(beta+1.96*se),2) 
    low<-round(exp(beta-1.96*se),2)  
    out[i]<-paste(names(glm)[i],"\t",fit," (",low,"-",up,")\t",pvalue)  
  }
  write.table(out, "marginal.pvalue_OR.txt",sep="\t",quote=F,row.names=T,col.names=NA)
}
CvSampling<- function(Nobs=1000,K=5){
  #======================
  # return sample row number: used to sampling the subset to cross-validation predition.
  # rank<-CvSampling(Nobs=dim(data)[1],K=5)
  # train=x[[1]]$train
  # test=x[[1]]$test
  #====================
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))], test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}
chan<-function(x){
  x[x<0.5]<-0; x[x>=0.5]<-1
  x
}
assess<-function(model,train,test,trainlabel,testlabel){
  t1<-table(chan(predict(model, train)), trainlabel)
  t2<-table(chan(predict(model, test)), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
assess2<-function(model1,model2,trainlabel,testlabel){
  t1<-table(model1,trainlabel)
  t2<-table(model2,testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
assess3<-function(model,train,test,trainlabel,testlabel){
  t1<-table((predict(model, train)), trainlabel)
  t2<-table((predict(model, test)), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
assess4<-function(model1,model2,trainlabel,testlabel){
  t1<-table(chan(model1), trainlabel)
  t2<-table(chan(model2), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}


#######################################################################################
################################## Data input and description #########################
#######################################################################################
setwd("/home/sguo/Dropbox/Project/rheumatology/SerumLipid/result/dat")
case<-read.table("case.txt",sep="\t",as.is=F,head=T)
cont<-read.table("control.txt",sep="\t",head=T,as.is=F)
cont2<-read.table("continfo.txt",sep="\t",head=T,as.is=F)
head(case)
head(cont)
head(cont2)
con<-merge(cont2,cont)
head(con)
write.table(con,file="con.txt",col.names=NA,row.names=T,sep="\t",quote=F)
#con<-read.table("con.txt",head=T,sep="\t",as.is=F)
head(con)
summary(con)
c1<-case[,c(1,2,3,7:10,18)]
c2<-con[,c(1,2,3,12:16)]
summary(c1)
summary(c2)
colnames(c1)<-colnames(c2)
data<-rbind(c1,c2)
head(data)


#######################################################################################
################# uni and multi variates logistic regression  #########################
#######################################################################################
OR95CIPV(data[,4:7])
summary(data)
glm<-summary(glm(abs(as.numeric(phen)-2)~age+gender+tc+tg+hdl+ldl,data=data,family=binomial(link = "logit")))
pvalue=glm$coefficients[,4]
beta=glm$coefficients[,1]
se=glm$coefficients[,2]
fit<-round(exp(beta),2)
up<-round(exp(beta+1.96*se),2) 
low<-round(exp(beta-1.96*se),2)  
ci<-paste(round(fit,5)," (",low,"-",up,")",sep="")
pvalue=format(pvalue,digits=3)
output<-cbind(pvalue,ci)
output
write.table(output, "partial.pvalue_OR.txt",sep="\t",quote=F,row.names=T,col.names=NA)

#######################################################################################
############################## matrix anovar analyis ##################################
#######################################################################################

y1<-matrixanova(data.frame(na.omit(case[,2:17])))


head(data)
y1<-matrixanova(na.omit(data[,2:8]))
dat1<-cor(na.omit(sapply(data[,2:8], as.numeric)),use="pairwise.complete.obs")
rcorr<-rcorr(na.omit(sapply(data[,2:8], as.numeric)))
rcorr$r
rcorr$P
write.table(rcorr$r,file="CorrelationMatrix.R2.total.txt",sep="\t",row.names=T,col.names=NA)
write.table(rcorr$P,file="CorrelationMatrix.P.total.txt",sep="\t",row.names=T,col.names=NA)
getwd()

rcorr<-rcorr(na.omit(sapply(c1[,2:8], as.numeric)))
rcorr$r
rcorr$P
write.table(rcorr$r,file="CorrelationMatrix.R2.case.txt",sep="\t",row.names=T,col.names=NA)
write.table(rcorr$P,file="CorrelationMatrix.P.case.txt",sep="\t",row.names=T,col.names=NA)
getwd()

rcorr<-rcorr(na.omit(sapply(c2[,2:8], as.numeric)))
head(c2)
rcorr$r
rcorr$P
write.table(rcorr$r,file="CorrelationMatrix.R2.control.txt",sep="\t",row.names=T,col.names=NA)
write.table(rcorr$P,file="CorrelationMatrix.P.control.txt",sep="\t",row.names=T,col.names=NA)
getwd()


t=as.data.frame(sapply(data[,2:8],as.numeric))
head(t)
t$phen=abs(t$phen-2)
dat<-as.data.frame(impute.knn(data.matrix(t[,1:7]))$data)
head(dat)
aip<-log(dat$tc/dat$hdl,10)
dat$aip=aip
head(dat)
table(dat$phen)
glm<-summary(glm(phen~age+gender+aip,data=dat,family=binomial(link = "logit")))
glm
pvalue=glm$coefficients[,4]
beta=glm$coefficients[,1]
se=glm$coefficients[,2]
fit<-round(exp(beta),2)
up<-round(exp(beta+1.96*se),2) 
low<-round(exp(beta-1.96*se),2)  
ci<-paste(round(fit,5)," (",low,"-",up,")",sep="")
pvalue=format(pvalue,digits=3)
output<-cbind(pvalue,ci)
output




mean <-aggregate(dat, by=list(data$gender),FUN=mean, na.rm=TRUE)
range <-aggregate(dat, by=list(data$gender),FUN=IQR, na.rm=TRUE)  # interquartile range
print(mean)
print(range)

max(data$tg)

pairs(dat[1:7], main = "lipid profile of RA and Normal",pch = 21, bg = c("red", "blue")[unclass(dat$phen)])

c11<-as.data.frame(sapply(c1[,2:7],as.numeric))
c22<-as.data.frame(sapply(c2[,2:7],as.numeric))

head(c1)
head(c11)

mean1 <-aggregate(c11, by=list(c11$gender),FUN=mean, na.rm=TRUE)
mean2 <-aggregate(c22, by=list(c22$gender),FUN=mean, na.rm=TRUE)
print(mean1)
print(mean2)


head(c1)

#######################################################################################
################################### custom analysis  ##################################
#######################################################################################
library("Deducer")
library("Rcmdr")
library("ROCR")

t=as.data.frame(sapply(data[,2:8],as.numeric))
head(t)
t$phen=abs(t$phen-2)
dat<-as.data.frame(impute.knn(data.matrix(t[,1:7]))$data)
head(dat)

dat1<-dat
dat<-dat1

table(dat1$gender)
dat<-subset(dat,gender==1)
head(dat)

model.glm1 <- (glm(phen ~ gender + age,family=binomial(),data=dat,na.action=na.omit))
model.glm2 <- (glm(phen ~ gender + age + tc,family=binomial(),data=dat,na.action=na.omit))
model.glm3 <- (glm(phen ~ gender + age + tc +tg ,family=binomial(),data=dat,na.action=na.omit))
model.glm4 <- (glm(phen ~ gender + age + tc +tg + hdl,family=binomial(),data=dat,na.action=na.omit))
model.glm5 <- (glm(phen ~ gender + age + tc +tg + ldl,family=binomial(),data=dat,na.action=na.omit))
model.glm6 <- (glm(phen ~ gender + age + tc +tg + ldl + hdl,family=binomial(),data=dat,na.action=na.omit))

library("PredictABEL")
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
pred3 <- predRisk(model.glm3)
pred4 <- predRisk(model.glm4)
pred5 <- predRisk(model.glm5)
pred6 <- predRisk(model.glm6)

tiff("ROC.femaile.tiff")
#labels <- c(paste("model", c("A: AUC=0.674 (0.659-0.688)","B: AUC=0.734 (0.721-0.747)","C: AUC=0.779 (0.766-0.791)","D: AUC=0.915 (0.908-0.923)","E: AUC=0.965 (0.961-0.970)", "F: AUC=0.972 (0.968-0.976)")))
#plotROC(data=dat,cOutcome=7,predrisk=cbind(pred1,pred2,pred3,pred4,pred5,pred6),labels=labels)
plotROC(data=dat,cOutcome=7,predrisk=cbind(pred1,pred2,pred3,pred4,pred5,pred6))
dev.off()

length(pred1)
length(pred2)
length(pred3)
length(pred4)
length(pred5)
length(pred6)

par(bg = "white")


# obtain predicted risks
library("PredictABEL")
pred1 <- predRisk(model1)
pred2 <- predRisk(model2)
pred3 <- predRisk(model3)
pred4 <- predRisk(model4)
labels <- c(paste("model", c("A:AUC=0.81 (0.78-0.84)","B:AUC=0.80 (0.77-0.83)","C:AUC=0.75 (0.72-0.79)","D:AUC=0.73 (0.69-0.76)")))
png("ROC.png")
plotROC(data=dataset1[xrow,], cOutcome=1, predrisk=cbind(pred1,pred2,pred3,pred4), labels=labels)
dev.off()
library("pROC")
r1<-roc(dataset1[xrow,]$type~pred1)
r2<-roc(dataset1[xrow,]$type~pred2)
r3<-roc(dataset1[xrow,]$type~pred3)
r4<-roc(dataset1[xrow,]$type~pred4)

f1<-coords(r1,"b")[c(3,2)]
f2<-coords(r2,"b")[c(3,2)]
f3<-coords(r3,"b")[c(3,2)]
f4<-coords(r4,"b")[c(3,2)]

round(rbind(f1,f2,f3,f4),2)


length(dataset1$type)
length(pred1)

getwd()
table(case$ARA,case$gender)
table(case$RF,case$Stomach)

x<-data.frame(ids$anti.RNA.pol3,ids$gut)
y<-which(! apply(x,1,function(x) any(is.na(x))))
x[y,]
t<-oneway.test(x[m,1]~x[m,2],data=x)
t.test(x[m,1]~x[m,2])
boxplot(x[m,1]~x[m,2])
m<-which(x[,1]<20)

boxplot(x[,1])






#######################################################################################
################################### classification  ##################################
#######################################################################################
N=10  # times of cross validation
K=5  # fold cross validation
rlt1<-matrix(NA,N*K,6)
rlt2<-matrix(NA,N*K,6)
j<-1

data=ColNARemove(data)
head(data)
t=as.data.frame(sapply(data[,1:7],as.numeric))
head(t)
t$phen=abs(t$phen-2)
pheno <- t$phen
dat<-as.data.frame(impute.knn(data.matrix(t[,1:6]))$data)
head(dat)
for(n in 1:N){
  index<-CvSampling(dim(dat)[1],K)
  for(i in 1:K){
    train=dat[index[[i]]$train,]
    trainlab=pheno[index[[i]]$train]
    test=dat[index[[i]]$test,]
    testlab=pheno[index[[i]]$test]
    model <- randomForest(as.factor(trainlab)~., data=train,importance=T) # as.factor(y)
    predict1 <- predict(model, train)
    predict2 <- predict(model, test)
    t1<-table(predict1, trainlab)
    t2<-table(predict2, testlab)
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt1[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    
    # for rpart
    model <- rpart(as.factor(trainlab) ~.,data=train,method="class")
    predict1 <- predict(model, train,type="class")
    predict2 <- predict(model, test, type="class")
    sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
    spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
    sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
    spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
    accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
    accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
    rlt2[j,]<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
    print(j)
    j=j+1
    
  }
}

rlt1
colMeans(rlt1)
rlt2
colMeans(rlt2)

colnames(dat)=c("Gender","Age",'T-CHOL',"TG","HDL-C","LDL-C")
fit <- rpart(pheno ~ .,method="class", data=dat[,1:6])
100*fit$variable.importance/sum(fit$variable.importance)


printcp(fit) # display the results
plotcp(fit) # visualize cross-validation results
summary(fit) # detailed summary of splits

# plot tree
pdf("rpart.tree.pdf")
par<-par(mar=c(1,1,1,1))
plot(fit, uniform=TRUE,main="Classification Tree for Rheumatoid Arthritis with lipid profiles",cex.main=0.9,cex=0.9)
text(fit, use.n=TRUE, all=TRUE, cex=.7)
dev.off()
? text


colnames(rlt1)<-c("sen.train","spe.train", "accu.train","sen.test","spe.test","accu.test")
colnames(rlt2)<-c("sen.train","spe.train", "accu.train","sen.test","spe.test","accu.test")  
write.table(rlt1,file="lipidprotein.randomforest.Predition.Result.txt",col.names=NA,row.names=T,quote=F,sep="\t")
write.table(rlt2,file="lipidprotein.rpart.Predition.Result.txt",col.names=NA,row.names=T,quote=F,sep="\t")

##################################




source("class.glasso.r")
source("GSCPackage.R")
library("impute")
library("DAAG")
library("PredictABEL")
library("Deducer")
library("Rcmdr")
library("ROCR")
library("MASS")
library("e1071")
library("class")
library("BayesTree")
library("randomForest")
library("gbm")

# 1, logistic regression
# 2, KNN
# 3, SVM
# 4, bayestree
# 5, random forest
# 6, Boosting
# 7, Fuzzy Rule-based Systems



# missing value must be remove, no other choice
source("GSCPackage.R")
install.packages("gplots")
install.packages("KernSmooth")

data<-RawNARemove(data)
dat<-data[,2:7]
head(dat)
rlt<-impute.knn(data.matrix(dat),k = 3)
phen<-abs(as.numeric(as.factor(data$phen))-2)
phen
input<-data.frame(rlt$data,phen)

dat<-data[,2:(dim(data)[2]-1)]
N=100  # times of cross validation
K=5  # fold cross validation
rlt1<-rlt2<-rlt3<-rlt4<-rlt5<-rlt6<-rlt7<-rlt8<-matrix(NA,N*K,6)
j<-1
input<-data.frame(dat,phen)
for(n in 1:N){
  dat<-CvSampling(dim(input)[1],5)
  for(i in 1:K){
    train=input[dat[[i]]$train,]
    trainlabel<-as.factor(phen[dat[[i]]$train])
    test=input[dat[[i]]$test,]
    testlabel<-phen[dat[[i]]$test]
    ### logistic model
    model<-glm(trainlabel~.,family=binomial,data=train)
    predict(model,train,type="response")
    predict(model,test,type="response")
    predict(model, train)
    x<-predict(model, train)
    rlt1[j,]<-assess(model,train,test,trainlabel,testlabel)
    ## 2. KNN model
    model1<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$train,],cl=pheno[dat[[i]]$train],k=6)
    model2<-knn(train=input[dat[[i]]$train,],test=input[dat[[i]]$test,],cl=pheno[dat[[i]]$train],k=6)
    rlt2[j,]<-assess2(model1,model2,trainlabel,testlabel)    
    ## 3. N*K-fold cross-validation for SVM classification    
    model<-svm(x=train,y=trainlabel)
    rlt3[j,]<-assess3(model,train,test,trainlabel,testlabel)    
    ## 4. N*K-fold cross-validation for BayesTree
    model1 = (apply(pnorm(bart(train,trainlabel,train)$yhat.test),2,mean))
    model2 = (apply(pnorm(bart(train,trainlabel,test)$yhat.test),2,mean))
    rlt4[j,]<-assess4(model1,model2,trainlabel,testlabel)    
    ## 5. N*K-fold cross-validation for randomForest
    model <- randomForest(trainlabel~., data=train)
    rlt5[j,]<-assess3(model,train,test,trainlabel,testlabel)    
    ## 6. N*K-fold cross-validation for Boosting
    model <- gbm(trainlabel~.,data=train,distribution = "adaboost")
    model1 <- predict.gbm(model,train,100,type="response")
    model2 <- predict.gbm(model,test,100)
    rlt6[j,]<-assess2(model1,model2,trainlabel,testlabel)    
    
    j=j+1
    print(j)
  }
}
colMeans(rlt1)
colMeans(rlt2)
colMeans(rlt3)
colMeans(rlt4)
colMeans(rlt5)
colMeans(rlt6)


RINfun=function(yorig){
  yranks=rank(yorig)
  tempp=(yranks-.5)/(length(yranks))
  return(qnorm(tempp))
}

new<-data.frame(sapply(case[,1:17],function(x) RINfun(x)))
head(new)
p<-matrix(0,4,14)
v<-matrix(0,4,14)

for(i in 1:4){
  for(j in 1:14){
    fit<-summary(lm(new[,7+i]~new[,3+j]+age+Gender,data=data.frame(new)))
    p[i,j]<-fit$coefficients[2,4]  
    v[i,j]<-fit$coefficients[2,1]  
    
  }
}
rownames(p)<-c("ch","tg","hdl","ldl")
colnames(p)<-colnames(new)[4:17]
rownames(v)<-c("ch","tg","hdl","ldl")
colnames(v)<-colnames(new)[4:17]
write.table(p,file="correlation.between.lipid.rf.pvalue.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(v,file="correlation.between.lipid.rf.beta.txt",sep="\t",col.names=NA,row.names=T,quote=F)
