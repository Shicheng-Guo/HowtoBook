# Pan-cancer methylation biomarker project

##########################################################################################################################
###################################### Introduction and Usage  #########################################################
##########################################################################################################################

# this code is for differential methylation region detection
# Oct/17/2016

##########################################################################################################################
###################################### Function and Package Load #########################################################
##########################################################################################################################

PairTtestPValue<-function(data,x1,x2,pair=FALSE){
  data<-data.matrix(data)
  output<-matrix(NA,dim(data)[1],6)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(pair==TRUE){
      Valid<-nrow(na.omit(data.frame(data[i,x1],data[i,x2])))
    }else{
      Valid<-100
    }
    if( sum(!is.na(data[i,x1]))>=3 & sum(!is.na(data[i,x2]))>=3 & Valid>3){ 
      tmp1<-try(t.test((data[i,x1]),(data[i,x2]),paired=pair, na.action=na.omit))
      output[i,1]<-format(tmp1$p.value, scientific=TRUE)
      output[i,2]<-round(mean((data[i,x1]))-mean((data[i,x2])),3)
      output[i,3]<-round(mean((data[i,x1])),3)
      output[i,4]<-round(mean((data[i,x2])),3)
      output[i,5]<-round(sd(data[i,x1]),3)
      output[i,6]<-round(sd(data[i,x2]),3)
      print(i)
    }
  }
  rownames(output)<-rownames(data)
  output
}

##########################################################################################################################
###################################### Working Pipeline ##################################################################
##########################################################################################################################

setwd(getwd())
# ESCA
system("gdc-client download -m gdc_manifest.2016-10-17T22-20-50.162580-ESCA.tsv")
system("cp ./*/*.txt ./")

# STAD
system("gdc-client download -m Stomach.gdc_manifest.2017-08-25T02-05-56.749099.txt")
system("cp ./*/*.txt ./")

# CHOL
system("gdc-client download -m CHOL.gdc_manifest.2017-08-25T02-09-04.795450.txt")
system("cp ./*/*.txt ./")

# COAD
system("gdc-client download -m COAD.READ.gdc_manifest.2017-08-25T02-03-23.440521.txt")
system("cp ./*/*.txt ./")


library("stringr")
file<-list.files(pattern="jhu*")
data<-c()
for(i in file){
  tmp<-read.table(i,head=T,skip=1,row.names=1,sep="\t",check.names = FALSE,as.is=T)
  data<-cbind(data,tmp[,1])
  print(i)
}

#load("PancancerMethMatrix_March2016.RData")
#load("PancancerMethMatrix_March2016.Test.RData")
# colnames(data)<-unlist(lapply(colnames(data),function(x) gsub("[.]","-",x)))
rownames(data)<-rownames(tmp)
idv<-as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
colnames(data)<-idv
cancertype<-unique(unlist(lapply(file,function(x) unlist(strsplit(x,"_|.Human"))[2])))
save(data,file=paste(cancertype,"meth.RData",sep="."))

# Identify Paired Tumor-Adjacent Samples(01/11) 
IDV<-unique(as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*")))
pairidv<-c()
for (i in 1:length(IDV)){
  t1<-paste(IDV[i],"-01",sep="") 
  t2<-paste(IDV[i],"-11",sep="")
  if(all(any(grepl(t1,file)),any(grepl(t2,file)))){
    pairidv<-c(pairidv,t1,t2)
  }
}



# Note: CHOL only have two normal tissues. 
#===========================================================
# Pair-wise DMS test
if(length(pairidv)>3){
  newdata<-data[,match(pairidv,colnames(data))]
  newdata<-newdata+matrix(rnorm(nrow(newdata)*ncol(newdata),0.0001,0.0001),dim(newdata)[1],dim(newdata)[2])   # row=gene, col=inv
  type<-substr(colnames(newdata),14,15)
  x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
  x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
  Rlt1<-PairTtestPValue(newdata,x1,x2,pair=TRUE)
}
# non-pair-wise DMS test
newdata<-data+matrix(rnorm(nrow(data)*ncol(data),0.0001,0.0001),dim(data)[1],dim(data)[2])   # row=gene, col=inv
type<-substr(colnames(data),14,15)
#  x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive (not stable since 06,09 will be existed in filename)
#  x2<-which(type==names(table(type))[2])   # type 2, normal or resistant (not stable since 06,09 will be existed in filename)
x1<-which(type=="01")   # type 1, cancer or sensitive
x2<-which(type=="11")   # type 2, normal or resistant
Rlt2<-PairTtestPValue(newdata,x1,x2,pair=F)
#===========================================================
cancertype<-unique(unlist(lapply(file,function(x) unlist(strsplit(x,"_|.Human"))[2])))
type<-substr(colnames(data),14,15)
x1<-which(type=="01")   # type 1, cancer or sensitive
x2<-which(type=="11")   # type 2, normal or resistant

xtable<-table(c(which(rowMeans(data[,x1])>0.3),which(rowMeans(data[,x2])<0.2)))
newdata<-data.frame(cancer=rowMeans(data[,x1]),normal=rowMeans(data[,x2]))
newdata<-newdata[as.numeric(names(xtable[which(xtable==2)])),]
head(newdata)
dim(newdata)

# remove PBMC hypermethylated regions
pbmc<-read.table("/media/NAS3_volume2/shg047/HM450/TCGA/Normal.PBMC.GEO.HM450K.Beta.txt",sep="\t",row.names=1,head=T,check.names=F)
PBMC<-pbmc[match(rownames(newdata),rownames(pbmc)),]
PBMCSubset<-subset(PBMC,mean<0.1 & median<0.1)
dim(PBMCSubset)
target<-newdata[match(rownames(PBMCSubset),rownames(newdata)),]
dim(target)
# remove Normal tissue hypermethylated regions
# GSI
# Number of significant in all the cancer
##########################################################################################################################
###################################### Annotation and Plot################################################################
##########################################################################################################################
# Annotation 
hm450anno<-read.table("/media/Home_Raid1/shg047/db/hg19/GPL13534.map")
result<-data.frame(target,hm450anno[match(rownames(target),hm450anno[,4]),])
dim(result)
outputfile=paste("TCGA-",cancertype,".DMS.txt",sep="")
write.table(result,file=outputfile,col.names=T,row.names=F,sep="\t",quote=F)
