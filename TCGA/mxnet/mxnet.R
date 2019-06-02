# case-control prediction

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
