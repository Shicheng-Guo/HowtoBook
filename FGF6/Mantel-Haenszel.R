The Haldane OR for the male stratum is 25.83.  The Haldane OR for the female stratum is 35.98. 

The joint OR across both strata is 28.52.  The Mantel-Haenszel test of homogeneity across strata is p=0.728. 

Therefore, there’s no evidence of significant effect differences between males and females.

```
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/hemochromatosis/Manuscript/round1/bak")
list.files()
data<-read.table("FGF6-PheTyp7_Iron_C1.phen.txt",head=T,sep="\t")
head(data)
F1<-subset(data,sex==2)
M1<-subset(data,sex==1)
dim(F1)
dim(M1)
table(M1$phen,M1$FGF6)
table(F1$phen,F1$FGF6)
```

