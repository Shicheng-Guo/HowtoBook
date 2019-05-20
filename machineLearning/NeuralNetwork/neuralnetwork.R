setwd("")
# The Boston dataset is a collection of data about housing values in the suburbs of Boston. 
# Our goal is to predict the median value of owner-occupied homes (medv) using all the other
# continuous variables available.
set.seed(500)
library(MASS)
data <- Boston
head(data)
apply(data,2,function(x) sum(is.na(x)))
index <- sample(1:nrow(data),round(0.75*nrow(data)))
train <- data[index,]
test <- data[-index,]
lm.fit <- glm(medv~., data=train)
summary(lm.fit)
pr.lm <- predict(lm.fit,test)
MSE.lm <- sum((pr.lm - test$medv)^2)/nrow(test)

# normalization to xlim=c(0,1) or xlim=c(-1,1)
maxs <- apply(data, 2, max) 
mins <- apply(data, 2, min)
scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
train_ <- scaled[index,]
test_ <- scaled[-index,]

library(neuralnet)
n <- names(train_)
f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)
plot(nn)

pr.nn <- compute(nn,test_[,1:13])
pr.nn_ <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)
test.r <- (test_$medv)*(max(data$medv)-min(data$medv))+min(data$medv)
MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(test_)
print(paste(MSE.lm,MSE.nn))

par(mfrow=c(1,2))

plot(test$medv,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(test$medv,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)


library("MASS")
library("readxl")
library("boot")
library("neuralnet")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC")
source("https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/assess.R")
beta<-read.table("methyLevel_54samples.txt",head=T,sep="\t")
phen<-data.frame(read_xlsx("ROC.xlsx",sheet=1))
rownames(beta)<-paste(beta[,1],":",beta[,2],"-",beta[,3],sep="")
beta<-beta[,4:ncol(beta)]
input<-na.omit(data.frame(phen=phen$phen,t(beta)))
set.seed(200)
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
glm.fit <- glm(phen~.,data=input,family=binomial)
xx<-cv.glm(input,glm.fit,cost,K=10)$delta

input<-input[,1:7]
colnames(input)<-c("Phenotypes",'Genetics',"Epigenetics","Environmental","Lifestyle","Image","Pathological")
head(input)
set.seed(450)
cv.error <- NULL
k <- 10
library("plyr") 
pbar <- create_progress_bar('text')
pbar$init(k)
n <- names(input)
f <- as.formula(paste("Phenotypes ~", paste(n[!n %in% "Phenotypes"], collapse = " + ")))
rlt<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  nn <- neuralnet(f,data=train.cv,hidden=c(5,2),act.fct = "logistic",linear.output = FALSE)
  plot(nn,lwd=0.85,cex=1.2)
  pr.nn <- compute(nn,test.cv)
  rlt<-rbind(rlt,data.frame(test.cv[,1],pr.nn$net.result))  
  pbar$step()
}



install.packages("party")
install.packages("reprtree")
library(randomForest)
library(party)
library(reprtree)
rf <- randomForest(Species ~ ., data=iris)
getTree(rf, 1, labelVar=TRUE)
cforest(Species ~ ., data=iris, controls=cforest_control(mtry=2, mincriterion=0))
plot(rf)
reprtree:::plot.getTree(rf)



