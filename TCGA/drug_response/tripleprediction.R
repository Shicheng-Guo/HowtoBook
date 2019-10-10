load("~/hpc/project/TCGA/pancancer/FPKM/mRNA.t3000.mxnet.RData")
mRNA<-rbind(mRNA$train.cv,mRNA$test.cv)
load("~/hpc/project/TCGA/pancancer/meth450/meth.t6000.mxnet.RData")
meth<-rbind(meth$train.cv,meth$test.cv)
load("~/hpc/project/TCGA/pancancer/miRNA/pancancer.miRNA.drugResponse.RData")

triple<-names(which(sort(table(c(rownames(miRNA),rownames(meth),rownames(mRNA))))==3))

miRNA<-miRNA[match(triple,rownames(miRNA)),]
mRNA<-mRNA[match(triple,rownames(mRNA)),]
meth<-meth[match(triple,rownames(meth)),]

data<-cbind(miRNA,mRNA,meth)
input<-cbind(miRNA[,1],mRNA[,2:ncol(mRNA)],meth[,2:ncol(meth)])
input<-input[,c(1,sample(1:ncol(input),250))]
input[1:5,1:5]
colnames(input)
index <- sample(1:nrow(input),round(0.9*nrow(input)))
train.cv <- data.matrix(input[index,])
test.cv <-  data.matrix(input[-index,])

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

