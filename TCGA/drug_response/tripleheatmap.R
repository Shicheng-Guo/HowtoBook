# CHG1
library("randomForest")
library("arm")
library("plyr") 
library("PredictABEL")
library("neuralnet")
library("caret")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
source("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/bin/id2phen4.R")

setwd("~/hpc/project/TCGA")

# methylation
load("methdata.pancancer.nomissing.RData")
colnames(input)<-id2phen4(colnames(input))
input<-input[,grep("-01",colnames(input))]
phen<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/pancancer.chemotherapy.response.txt",head=T,sep="\t")
phen$ID4<-paste(phen$bcr_patient_barcode,"-01",sep="")
input<-input[,colnames(input) %in% phen$ID4]
input<-input[,match(unique(colnames(input)),colnames(input))]
rx<-findCorrelation(t(input), cutoff = 0.8, names = F,exact = ncol(t(input)) < 100)
input<-input[-rx,]
phen<-phen[na.omit(unlist(lapply(colnames(input),function(x) match(x,phen$ID)[1]))),]
dim(input)
dim(phen)
sort(table(phen$bcr_patient_barcode))
table(phen$measure_of_response)
levels(phen$measure_of_response)<-c(0,1,1,0)
input<-data.frame(phen=phen$measure_of_response,t(input))

P=apply(input[,2:ncol(input)],2,function(x) summary(glm(as.factor(input[,1])~x,family=binomial))$coefficients[2,4])
input<-input[,c(1,match(names(P[head(order(P),n=2000)]),colnames(input)))]

RF <- randomForest(as.factor(phen) ~ ., data=input, importance=TRUE,proximity=T)
imp<-RF$importance
head(imp)
imp<-imp[order(imp[,4],decreasing = T),]
topvar<-match(rownames(imp)[1:2000],colnames(input))
newinput <- t(input[,topvar])
colnames(newinput)<-input[,1]
newinput[1:5,1:5]
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
pdf("meth.heatmap.randomForest.n2.pdf")
HeatMap(newinput)
dev.off()

# mRNA
load("Pancancer.DrugResponse.V5292.N1462.RData")
input<-newinput


# miRNA

        
        
