```
library("ChAMP")
library("outliers")
setwd("/home/sguo/dyh/idat")
myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
save(myLoad,file="myLoad.RData")
# turn off the plot if you run it on the server since of the problem of X11
myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = "sampleSheet.txt", resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
save(myNorm,file="myNorm.RData")
TMR<-champ.TrueMethyl(beta.norm = myNorm$beta, pd = myLoad$pd, adjPVal = 0.5, adjust.method = "BH",compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = TRUE)
```
