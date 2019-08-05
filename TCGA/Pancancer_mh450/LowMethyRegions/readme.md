```
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/cpg2bed.R")
setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
load("methdata.pancancer.RData")
data<-na.omit(methdata)
Median<-unlist(apply(data,1,function(x) median(x)))
cpgs<-rownames(data)[which(Median<0.2)]
bed<-cpg2bed(cpgs)
```
