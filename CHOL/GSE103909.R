wget -r -l 1 -nd -e robots=off --reject jpg,html ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103909/suppl/

library("GEOquery")
GSE103909 <- getGEO("GSE103909")
data <- as.data.frame(exprs(GSE103909[[1]]))
phen <- pData(phenoData(GSE103909[[1]]))

phen1<-sapply(strsplit(as.character(phen$characteristics_ch1.7),"[:]"),function(x) as.numeric(unlist(x)[2]))  # status 1:control, 2:scz
phen1[phen1==1]<-"Normal"
phen1[phen1==2]<-"schizophrenia"
phen2<-sapply(strsplit(as.character(phen$characteristics_ch1),"[:]"),function(x) (unlist(x)[2]))  # gender

data1=na.omit(data)
PCAPlot(t(data1),phen1,output="GSE41169.scz.normal.pdf",multifigure=T)  # status
PCAPlot(t(data1),phen2,output="GSE41169.gender.pdf",multifigure=T)  # gender
