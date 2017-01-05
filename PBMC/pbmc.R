library("GEOquery")
# PBMC-111 samples
load("/home/shg047/oasis/monod/GEO/GSE53045_matrix.Rdata")
GSE53045Beta <- as.data.frame(exprs(GSE530451))
phen <- pData(phenoData(GSE530451))
GSE53045NormalPBMC<-GSE53045Beta
dim(GSE53045NormalPBMC)
# PBMC-60 samples
load("/home/shg047/oasis/monod/GEO/GSE35069_matrix.Rdata")
GSE35069Beta <- as.data.frame(exprs(GSE350691))
phen <- pData(phenoData(GSE350691))
GSE35069NormalPBMC<-GSE35069Beta
dim(GSE35069NormalPBMC)
# PBMC-20 samples
load("/home/shg047/oasis/monod/GEO/GSE32148_matrix.Rdata")
GSE32148Beta <- as.data.frame(exprs(GSE321481))
phen <- pData(phenoData(GSE321481))
GSE32148NormalPBMC<-GSE32148Beta[,which(phen$description=="Normal peripheral blood sample")]
dim(GSE32148NormalPBMC)
# PBMC-192
load("/home/shg047/oasis/monod/GEO/GSE36054_matrix.Rdata")
GSE36054Beta <- (exprs(GSE360541))
GSE36054Beta <- as.data.frame(GSE36054Beta)
phen <- pData(phenoData(GSE360541))
GSE36054NormalPBMC<-GSE36054Beta
dim(GSE36054NormalPBMC)
# PBMC-78
load("/home/shg047/oasis/monod/GEO/GSE36064_matrix.Rdata")
GSE36064Beta <- as.data.frame(exprs(GSE360641))
phen <- pData(phenoData(GSE360641))
GSE36064NormalPBMC<-GSE36064Beta
dim(GSE36064NormalPBMC)
# PBMC-689
load("/home/shg047/oasis/monod/GEO/GSE42861_matrix.Rdata")
GSE42861Beta <- as.data.frame(exprs(GSE428611))
phen <- pData(phenoData(GSE428611))
GSE42861NormalPBMC<-GSE42861Beta[,which(phen$characteristics_ch1.1=="disease state: Normal")]
dim(GSE42861NormalPBMC)
data<-cbind(GSE53045NormalPBMC,GSE35069NormalPBMC,GSE32148NormalPBMC,GSE36054NormalPBMC,GSE36064NormalPBMC,GSE42861NormalPBMC)
tmp<-t(apply(data,1,function(x) c(cpg=rownames(x),mean=mean(x,na.rm=T),
median=median(x,na.rm=T),
SD=sd(x,na.rm=T),
Quantile=quantile(x,na.rm=T),
SampleSize=length(na.omit(x))
)
)
)
write.table(tmp,file="Normal.PBMC.GEO.HM450K.Beta.txt",col.names=NA,row.names=T,sep="\t",quote=F)
