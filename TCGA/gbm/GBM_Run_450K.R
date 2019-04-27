library("GEOquery")
GSE41826 <- getGEO("GSE41826")
data <- as.data.frame(exprs(GSE41826[[1]]))
phen <- pData(phenoData(GSE41826[[1]]))
phe<-data.frame(rownames(phen),phen$title,phen$`tissue:ch1`)
GSE41826<-list()
GSE41826$beta=data
GSE41826$phen=phen
save(GSE41826,file="GSE41826.RData")

GSE89707 <- getGEO("GSE89707")
data <- as.data.frame(exprs(GSE89707[[1]]))
phen <- pData(phenoData(GSE89707[[1]]))
GSE89707<-list()
GSE89707$beta=data
GSE89707$phen=phen
save(GSE89707,file="GSE89707.RData")

GSE74486 <- getGEO("GSE74486")
data <- as.data.frame(exprs(GSE74486[[1]]))
phen <- pData(phenoData(GSE74486[[1]]))
GSE74486<-list()
GSE74486$beta=data
GSE74486$phen=phen
save(GSE74486,file="GSE74486.RData")

GSE103659 <- getGEO("GSE103659")
data <- as.data.frame(exprs(GSE103659[[1]]))
phen <- pData(phenoData(GSE103659[[1]]))
GSE103659<-list()
GSE103659$beta=data
GSE103659$phen=phen
save(GSE103659,file="GSE103659.RData")

GSE114534 <- getGEO("GSE114534")
data <- as.data.frame(exprs(GSE114534[[1]]))
phen <- pData(phenoData(GSE114534[[1]]))
GSE114534<-list()
GSE114534$beta=data
GSE114534$phen=phen
save(GSE114534,file="GSE114534.RData")

GSE66351 <- getGEO("GSE66351")
data <- as.data.frame(exprs(GSE66351[[1]]))
phen <- pData(phenoData(GSE66351[[1]]))
GSE66351<-list()
GSE66351$beta=data
GSE66351$phen=phen
save(GSE66351,file="GSE66351.RData")


dim(GSE41826$beta)
dim(GSE89707$beta)
dim(GSE74486$beta)
dim(GSE103659$beta)
dim(GSE114534$beta)
dim(GSE66351$beta)

dim(GSE41826$phen)
dim(GSE89707$phen)
dim(GSE74486$phen)
dim(GSE103659$phen)
dim(GSE114534$phen)
dim(GSE66351$phen)

cas<-grep("Tumor",phe[,2])
con<-grep("_Normal",phe[,2])
pvalue<-apply(data,1,function(x) t.test(x[cas],x[con],paired=T)$p.value)
tvalue<-apply(data,1,function(x) t.test(x[cas],x[con],paired=T)$statistic)

output<-data.frame(data,tvalue,pvalue)
output<-subset(output,pvalue<0.0005)
dim(output)
write.table(output,file="CCA.lncRNA.diff.txt",sep="\t",quote=F,col.names = NA,row.names = T)
