````
# for ips methylatin 450K analysis
library("GEOquery")
GSE39279 <- getGEO("GSE39279")
data1 <- as.data.frame(exprs(GSE39279[[1]]))
phen1 <- pData(phenoData(GSE39279[[1]]))

GSE16559 <- getGEO("GSE39279")
data2 <- as.data.frame(exprs(GSE39279[[1]]))
phen2 <- pData(phenoData(GSE39279[[1]]))

GSE85566 <- getGEO("GSE85566")
data3 <- as.data.frame(exprs(GSE85566[[1]]))
phen3 <- pData(phenoData(GSE85566[[1]]))

phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"schizophrenia"
data1=na.omit(data)
PCAPlot(t(data1),phen1,output="GSE41169.scz.normal.pdf",multifigure=T)  # status
PCAPlot(t(data1),phen2,output="GSE41169.gender.pdf",multifigure=T)  # gender
```
