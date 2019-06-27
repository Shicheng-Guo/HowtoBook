set.seed(49)
library("randomForest")
library("arm")
library("plyr") 
cv.error <- NULL
k <- 10
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]

  RF <- randomForest(as.factor(phen) ~ ., data=input, importance=TRUE,proximity=T)
  imp<-RF$importance
  head(imp)
  imp<-imp[order(imp[,4],decreasing = T),]
  topvar<-match(rownames(imp)[1:10],colnames(input))
  
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]

  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,hidden=c(3),act.fct = "logistic",linear.output = FALSE)
  plot(nn,lwd=0.85,cex=1)
  pr.nn <- compute(nn,test.cv)
  rlt1<-rbind(rlt1,data.frame(phen=train.cv[,1],pred=nn$net.result))  
  rlt2<-rbind(rlt2,data.frame(phen=test.cv[,1],pred=pr.nn$net.result))  
}

