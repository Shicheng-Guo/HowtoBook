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

