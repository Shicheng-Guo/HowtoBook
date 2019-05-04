#########################  Simulation  ######################################


#########################  Simulation  ######################################
library("gplots")

color.map <-function(cancertype){
    cancertype<-as.numeric(cancertype)
    cancertype[is.na(cancertype)]<-0
    cancertype
  }




#########################  Subfunction  ######################################
#########################  Subfunction  ######################################


# missing remove and imputation


# data operation
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# data-mining

ClustValidtion<-function(ClusteResult,ClinicNumericFact){
  library("fpc")
  rlt<-list()
  ClusteResult<-as.matrix(ClusteResult)
  ClinicNumericFact<-as.matrix((ClinicNumericFact))
  data<-cbind(ClinicNumericFact,ClusteResult)
  data2<-RawNARemove(data)
  ClinicNumericFact<-data2[,1:dim(ClinicNumericFact)[2]]
  ClusteResult<-data2[,(dim(as.matrix(ClinicNumericFact))[2]+1):(dim(as.matrix(ClinicNumericFact))[2]+dim(ClusteResult)[2])]
  CorrectRand<-Vi<-matrix(NA,dim(ClusteResult)[2],dim(as.matrix(ClinicNumericFact))[2])
  for (i in 1:dim(ClusteResult)[2]){
    for (j in 1:dim(as.matrix(ClinicNumericFact))[2]){
      result<-cluster.stats(NULL,ClinicNumericFact[,j],ClusteResult[,i],compareonly=T)
      CorrectRand[i,j]<-result$corrected.rand
      Vi[i,j]<-result$vi
    }
  }
  rlt$correctRand<-CorrectRand
  rlt$vi<-Vi
  rlt
}
skmean<-function(data, classnum){
  # suitable to big data especially for col > row matrix, defort to cluste matrix to maximun 10 clusters
  cluster<-matrix(NA, dim(data)[1], classnum-1)
  for (i in 2:classnum){
    ask<-skmeans(data,i)
    file<-paste(data,".sk",i,".RData",sep="")
    cluster[,i-1]<-ask$cluster
    save(ask, file=file)
    print(i)
  }
  cluster
}
sparHierarchiCl<-function(mydata,method="complete",cut=2){
  library("sparcl")
  rlt<-list()
  set.seed(1)
  # Do tuning parameter selection for sparse hierarchical clustering
  perm.out <- HierarchicalSparseCluster.permute(data.matrix(mydata), wbounds=c(1.1,seq(2,15,by=1)),nperms=5)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw,method=method)
  rlt$sparsedata <-mydata[,which(sparsehc$ws !=0)] 
  #drawHeatmap2(m)
  d <- dist(data.matrix(rlt$sparsedata),method = "euclidean") #distance matrix
  fit <- hclust(d, method="ward")
  #plot(fit) # display dendogram
  clusters <- cutree(fit, cut) #cut tree into 5 clusters
  #draw dendogram with red borders around the 5 cluster
  #rect.hclust(fit, cut, border="red") 
  #retrun output
  rlt$featurenumber<-sum(sparsehc$ws !=0)
  rlt$weight<-sparsehc$ws
  rlt$group<-clusters;
  rlt$sparsehc=sparsehc
  rlt$bestw<-perm.out$bestw
  rlt$ws<-sparsehc$ws
  rlt
}
sparKmeanCl<-function(mydata,cut=2){
  library("sparcl")
  set.seed(1)
  km.perm <- KMeansSparseCluster.permute(mydata,cut,wbounds=seq(1.1,10,len=15),nperms=5)
  # run sparse k-means
  km.out <- KMeansSparseCluster(mydata,cut,wbounds=km.perm$bestw)
  rlt<-list()
  rlt$sparsedata<-mydata[,which(km.out[[1]]$ws !=0)]
  rlt$featurenumber<-sum(km.out[[1]]$ws !=0)
  rlt$sparsekm<-km.out
  rlt$group<-km.out[[1]]$Cs
  rlt$ws<-km.out[[1]]$ws
  rlt
}

sdsparse<-function(dat,quantile){
  sd<-apply(dat,2,function(x) sqrt(var(x,na.rm=T))/mean(x,na.rm = T))
  dat<-dat[,which(sd>quantile(sd,quantile,na.rm=T))]
  dat
}
tauSparse<-function(dat,pheno,tau=0.8){  
  library(gplots)
  data<-as.list(data.frame(dat))
  aggre<-aggregate(data,by=list(pheno),FUN=median)
  tau2<-as.numeric();
  for (i in 2:(ncol(aggre)-1)){
    tau1<-(nrow(aggre)-sum(aggre[,i])/max(aggre[,i]))/(nrow(aggre)-1)  
    tau2<-c(tau2,tau1)
  }
  
  color.map <-function(cancertype){ 
    cancertype<-as.numeric(cancertype)
    cancertype[is.na(cancertype)]<-0
    cancertype
  }
  
  patientcolors <- as.character(unlist(lapply(pheno, color.map)))
  data<-as.data.frame(data)
  pdf("heatmap.pdf",width=17,height=17)  #increase it when fig is big
  heat<-heatmap.2(as.matrix(data[,which(tau2>tau)]),col=topo.colors(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors)
  dev.off()
  hightaucolumn<-data[,which(tau2>tau)]
  hightaucolumn
}
CharToNumberFactor<- function(Variable){ 
  Variable<-as.numeric(Variable)
  Variable[is.na(Variable)]<-0
  Variable
}
pvalueSparse<-function(dat,pheno,phenocode1,phenocode2,pvalue=0.005){  
  library(gplots)
    p<-apply(dat,2, function(x) wilcox.test(x[which(pheno==phenocode1)], x[which(pheno==phenocode2)])$p.value)  
    color.map <- function(cancertype){ 
    cancertype<-as.numeric(cancertype)
    cancertype[is.na(cancertype)]<-0
    cancertype
  }
  patientcolors <- as.character(unlist(lapply(phen, color.map)))
  data<-as.data.frame(dat)
  pdf("heatmap.v=0.00005.pdf",width=17,height=17)  #increase it when fig is big
  heat<-heatmap.2(as.matrix(data[,which(p<0.005)]),col="heat.colors", scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors,main="p-value=0.005")
  dev.off()
  hightaucolumn<-data[,which(p<pvalue)]
  hightaucolumn
}
frr<-function(rr1,rr2,p,q,lamdi,frr){
#familial reltaive risk
r1<-10   #relative risk(estimated by odds ratios) for hetrozygotes relative to common homozygotes
r2<-10   #relative risk(estimated by odds ratios) for rare homozygotes to common homozygotes
p<-0.3   #risk allelefrequency
q<-1-p
lamdi0<-8.48  #overall familial relatve risk 
frr<-(p*(p*r2+q*r1)^2+q*(p*r1+q*r2)^2)/(p^2*r2+2*p*q*r1+q^2)^2
proportion<-log(lamdi)/log(lamdi0)
frr
proportion
}
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
}   

ColNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[1]
  NaCol<-which(apply(data,2,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,2,function(x) all(x==0))==T)
  NaCOL<-c(NaCol,zero)
  if(length(NaCOL)>0){
    data1<-data[,-NaCOL]
  }else{
    data1<-data;
  }
  data1
}   

Spmarker<-function(mfile,pheno,tau=0.8,type='txt'){  
  library(gplots)
  library(impute)
  #usage:
  #source("Spmarker.R")
  #Spmarker("mtcars.txt","saminfo.txt",tau=0.3)
  #ע?⣺ ????tau????̫?󣺻???ʾ?? `x' must be a numeric matrix
  #saminfo.txt ?б?????Category??Ϊ????
  
  cat('Reading Molecular Data File\n')
  dat <- read.table(mfile,header=T,comment.char='#',fill=T,sep='\t', as.is=T)
  cat('Reading Sample Information File\n')
  saminfo <- read.table(pheno, header=T, sep='\t',comment.char='#')
  if(sum(colnames(saminfo)=="Category")!=1){return('ERROR: Sample Information File does not have a Category column!')}
  
  dat<-impute.knn(as.matrix(dat[,2:ncol(dat)]))   
  Category<-saminfo$Category
  data<-as.list(data.frame(dat$data,Category))
  aggre<-aggregate(data,by=list(Category),FUN=median)
  tau2<-as.numeric();
  for (i in 2:(ncol(aggre)-1)){
    tau1<-(nrow(aggre)-sum(aggre[,i])/max(aggre[,i]))/(nrow(aggre)-1)	
    tau2<-c(tau2,tau1)
  }
  
  color.map <- function(cancertype){ 
    cancertype
  }
  patientcolors <- as.character(unlist(lapply(saminfo$Category, color.map)))
  
  data<-as.data.frame(data)
  pdf("heatmap.pdf",width=17,height=17)  #increase it when fig is big
  heat<-heatmap.2(as.matrix(data[,which(tau2>tau)]),col=topo.colors(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors)
  dev.off()
}



# cluster analysis and validation
fhclust<-function(mydata,method="euclidean",cut=4,bootstrap=FALSE){
  # 1, Ward Hierarchical Clustering(uni-cluster for row)  Jan 31st 2013
  d <- dist(mydata, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  # plot(fit) # display dendogram
  clusters <- cutree(fit, k=cut) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 cluster
  # rect.hclust(fit, k=cut, border="red") 
  if(bootstrap==TRUE){
    library(pvclust)
    # Ward Hierarchical Clustering with Bootstrapped p values
    fit <- pvclust(mydata, method.hclust="ward",method.dist="euclidean")
    plot(fit) # dendogram with p values
    # add rectangles around groups highly supported by the data
    pvrect(fit, alpha=.95)   
  }
  list<-names(table(clusters))
  rlt<-list()
  rlt<-lapply(list,function(x) names(which(clusters==x)))
  rlt
}
skmean<-function(data, classnum){
  # 2, K-Means Clustering with 5 clusters  Jan 31st 2013
  # suitable to big data especially for col > row matrix, defort to cluste matrix to maximun 10 clusters
  cluster<-matrix(NA, dim(data)[1], classnum-1)
  for (i in 2:classnum){
    ask<-skmeans(data,i)
    file<-paste(data,".sk",i,".RData",sep="")
    cluster[,i-1]<-ask$cluster
    save(ask, file=file)
    print(i)
  }
  cluster
}
fmclust<-function(mydata){
  # 3, Model Based Clustering   Jan 31st 2013
  library(mclust)
  fit <- Mclust(mydata)
  plot(fit, "BIC") # plot results
  print(fit) # display the best model 
  rlt<-list()
  rlt$group<-fit$classification
}
sparcl<-function(mydata,method="complete",cut=4){
  # 4, select feature with laso and then do hierarchical cluster to treat spare problem  #Feb 1 2013
  library("sparcl")
  set.seed(1)
  # Do tuning parameter selection for sparse hierarchical clustering
  perm.out <- HierarchicalSparseCluster.permute(x, wbounds=c(1.5,2:9),nperms=5)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw,method="complete")
  m<-x[,which(sparsehc$ws !=0)]
  #drawHeatmap2(m)
  d <- dist(m, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  #plot(fit) # display dendogram
  clusters <- cutree(fit, k=cut) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 cluster
  rect.hclust(fit, k=cut, border="red") 
  #retrun output
  rlt<-list()
  rlt$featurenumber<-sum(sparsehc$ws !=0)
  rlt$group<-clusters;
  
} 
fbiclust<-function(mydata){
  # 1, biclust Clustering
  library("biclust")
  set.seed(1)
  bics <- biclust(mydata,BCPlaid(), back.fit = 2, shuffle = 3, fit.model = ~m + a + b,iter.startup = 5, iter.layer = 30, verbose = TRUE)
  rlt<-list()
  r1<-fbicorder(bics, cols=FALSE, rev=FALSE)
  r2<-fbicorder(bics, cols=TRUE, rev=FALSE)
  rlt$row=r1
  rlt$col=r2
  rlt
}
fisa2<-function(mydata){
  # 2, isa2 cluster analysis 
  library("isa2")
  isa.result <- isa(mydata)
  
  bc <- isa.biclust(isa.result)
  # drawHeatmap(mydata, bc, 1)
  #bubbleplot(mydata, bc)
  parallelCoordinates(mydata, bc, number=1)
  
  names(isa.result)
  ## Find the best bicluster for each block in the input
  best <- apply(cor(isa.result$rows, data[[2]]), 2, which.max)
  ## Check correlation
  sapply(seq_along(best),function(x) cor(isa.result$rows[,best[x]], data[[2]][,x]))
  ## The same for the columns
  sapply(seq_along(best),function(x) cor(isa.result$columns[,best[x]], data[[3]][,x]))
  ## Plot the data and the modules found
  if (interactive()) {
    layout(rbind(1:2,3:4))
    image(data[[1]], main="In-silico data")
    sapply(best, function(b) image(outer(isa.result$rows[,b],
                                         isa.result$columns[,b]),
                                   main=paste("Module", b)))
  }
}
clusterplot<-function(mydata,group1){
  #1, Centroid Plot against 1st 2 discriminant functions
  library(fpc)
  plotcluster(mydata, fit$cluster) 
}
clustercompare<-function(mydata,group1,group2){
  #2, compare different cluster reslut
  library(fpc)
  d<-dist(mydata)
  cluster.stats(d,group1 ,group2)  
}
fbicorder<-function(bicResult, cols=TRUE, rev=FALSE){
#### Function to order variables or objects that appear in a bicluster[function bicorder]
  i<-numeric()
  res<-c()
  order<-vector();
  if(!cols){
    le<-dim(bicResult@RowxNumber)[2]
    for(i in 1:le){
      order[which(bicResult@RowxNumber[,i])[!(which(bicResult@RowxNumber[,i]) %in% res)]]<-i;
      res<-c(res,which(bicResult@RowxNumber[,i])[!(which(bicResult@RowxNumber[,i]) %in% res)])
    }
    count<-1:dim(bicResult@RowxNumber)[1]
    order[count[!(count %in% res)]]<-0;
  }else{
    le<-dim(bicResult@NumberxCol)[1]
    for(i in 1:le){
      order[which(bicResult@NumberxCol[i,])[!(which(bicResult@NumberxCol[i,]) %in% res)]]<-i;
      res<-c(res,which(bicResult@NumberxCol[i,])[!(which(bicResult@NumberxCol[i,]) %in% res)])
    }
    count<-1:dim(bicResult@NumberxCol)[2]
    order[count[!(count %in% res)]]<-0;
  }
  
  if(rev) order<-rev(order)
  order
}



# mean test 
MatrixTtest<-function(matrix,pheno,pvalue){ 
}
MatrixWilcoxtest<-function(matrix,pheno,pvalue){ 
}



#plot (Heatmap, Manhattan)
# scalefactorforcolor is used to scale the color so that the heatmap is more beautiful

heatmap.draw<-function(data,saturationvalue,pheno,heatmapfilename,scalefactorforcolor){
library("gplots")
patientcolors <- as.character(unlist(lapply(pheno, color.map)))
data<-as.data.frame(data)
data[abs(data)>saturationvalue]<-saturationvalue
pdf(heatmapfilename,width=17,height=17)  #increase it when fig is big
heat<-heatmap.2(as.matrix(scalefactorforcolor*data),col=redgreen(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors)
dev.off()
}

#Genetic Epidemilogy 




## Code for testing FPCA scores; run2
library(rrcov)
## Get the indices of genes functional PCA scores.  The header contains gene names and same gene's
## entries are consecutive.
group.by.gene <- function(header){
  last.key <- ""
  idx.set <- vector('integer')
  gene.ids <- list()
  idx <- 1
  i<-0
  for (name in strsplit(header, '_', fixed=TRUE)){
    i=i+1
    if(length(name)>2){
      next
    }
    key <- name[1] #paste(name[1], name[2], sep='_')
    val <- as.integer(name[length(name)])
    if (key == last.key){
      stopifnot(last.val+1 == val)
      idx.set[length(idx.set)+1] <- idx
    } else {
      stopifnot(1 == val)
      gene.ids[[last.key]] <- idx.set
      idx.set <- c(idx)
      last.key <- key
    }
    last.val <- val
    idx <- idx + 1
    #print(i)
  }
  gene.ids[[last.key]] <- idx.set
  gene.ids <- gene.ids[-1]
  return(gene.ids)
}
## Get the samples whose names are in a sample subset.
select.samples <- function(sample.names, subset.samples){
  without.suffix <- sapply(strsplit(sample.names, '-', fixed=TRUE), function(x){paste(x[1:3], collapse='-')})
  idx <- without.suffix %in% subset.samples
  return( sample.names[idx] )
}
## sample.file is assumed to be drug_response.csv 
#rlt<-compute.pval("","")
#write.table(rlt,file="",col.names=F)
#rlt<-compute.pval("FPCAScore.N.11.txt","scalefactor.txt")


compute.pval <- function(score.file, sample.file){
  ## separate FPCA scores into 2 groups
  fpca.scores <- read.table(score.file, header=TRUE, stringsAsFactor=FALSE)
  drug.response <- read.table(sample.file,head=T,sep="\t")
  case.samples <- subset(drug.response, Status == '1')$Sampleid
  con.samples <- subset(drug.response, Status == "11")$Sampleid
  case.idx <-(rownames(fpca.scores) %in% case.samples)
  con.idx <- (rownames(fpca.scores) %in%  con.samples)
  ## Group columns by genes and compute test statistics one-by-one
  genes <- group.by.gene(colnames(fpca.scores))
  pvals <- vector('numeric', length=length(genes))
  for (i in 1:length(genes)){
    cmat <- fpca.scores[, genes[[i]] ,drop=FALSE]
    tryCatch(pvals[i] <- T2.test(cmat[case.idx,], cmat[con.idx,] )$p.value,
             error = function(errMsg){
               print(paste("Error with gene", names(genes)[i], errMsg))
               pvals[i] <- NaN
             })
    print(i)
  }
  names(pvals) <- names(genes)
  return(pvals)
  print(length(genes))
}


#======================
#collect information when sparse cluster has been completed.
#ReCollection("FPCAScore.N.41.txt",sd.quantile=0.3,cut=2,heatmapdataquantile=0.97) #  good
#ReCollection("FPCAScore.N.31.txt",sd.quantile=0.3,cut=2,heatmapdataquantile=0.97) #  good
#ReCollection("FPCAScore.N.21.txt",sd.quantile=0.3,cut=2,heatmapdataquantile=0.97) #  good
#=====================
ReCollection<-function(input,sd.quantile=0.8,cut=2,heatmapdataquantile=0.99){
library("gplots")
library("sparcl")
rdata<-paste(input,".RData",sep="")
load(rdata)
out1<-paste("cluster.result",input,"hirar.sdq=",sd.quantile,"K=",cut,"RData",sep=".")
out2<-paste("cluster.result",input,"kmean.sdq=",sd.quantile,"K=",cut,"RData",sep=".")
load(out1)   #For Second Time
load(out2)   #For Second Time
result<-list()
result$cl<-as.numeric(substr(rownames(fpcascore),14,15))
result$cl1<-rlt1$group
result$cl2<-rlt2$group

# redraw heatmap with optional genes
out3<-paste(input,"hirar.sdq=",sd.quantile,"K=",cut,"Q>",heatmapdataquantile,"pdf",sep=".")
heatmap.draw(rlt1$sparsedata[,which(rlt1$ws>quantile(rlt1$ws,heatmapdataquantile))],result$cl,out3)
out4<-paste(input,"kmean.sdq=",sd.quantile,"K=",cut,"Q>",heatmapdataquantile,"pdf",sep=".")
heatmap.draw(rlt2$sparsedata,result$cl,out4)
}





#==========================
# use to plot heatmap, red is low expression, green is high expression. Be careful to set saturationvalue, if this value is very large, the resolution of the color is small.
#load("FPCAScore.N.11.txt.RData")
#phen<-as.numeric(substr(rownames(fpcascore),14,15))
#phenocode1="1"
#phenocode2="11"
#pvaluequatile=0.008
#heatmappdfname="Q0.01.pdf"
#pvalueSparse(fpcascore,saturationvalue=100,phen,"1","11",0.001,"Q0.01.pdf")
pvalueSparse<-function(dat,saturationvalue,phen,phenocode1,phenocode2,pvaluequatile,heatmappdfname){
  library(gplots)
  p<-apply(dat,2, function(x) try(wilcox.test(x[which(phen==phenocode1)], x[which(phen==phenocode2)], paired = TRUE)$p.value))
  color.map <- function(cancertype){
    cancertype<-as.numeric(cancertype)
    cancertype[is.na(cancertype)]<-0
    cancertype
  }
  patientcolors <- as.character(unlist(lapply(phen, color.map)))
  data<-as.data.frame(dat)
  pdf(heatmappdfname,width=17,height=17)  #increase it when fig is big
  input<-data[,which(p<quantile(p,pvaluequatile,na.rm=T))]
  input[(abs(input)>saturationvalue)]<-saturationvalue
  heat<-heatmap.2(data.matrix(input),col=redgreen(75), scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,keysize=0.5,RowSideColors=patientcolors)
  dev.off()
}
# this is funtion is special for TCGA Cancer-Normal Subset (load FunctionPCA RData to plot pvalue based heatmap)
#===============
#heat.map.draw("FPCAScore.N.11.txt",pvaluequatile=0.008,saturationvalue=30)
#heat.map.draw("FPCAScore.N.21.txt",pvaluequatile=0.008,saturationvalue=30)
#heat.map.draw("FPCAScore.N.31.txt",pvaluequatile=0.008,saturationvalue=30)
#heat.map.draw("FPCAScore.N.41.txt",pvaluequatile=0.008,saturationvalue=30)
#===============

heat.map.draw<-function(file,pvaluequatile,saturationvalue=100){
rdata<-paste(file,".RData",sep="")
load(rdata)
phen<-as.numeric(substr(rownames(fpcascore),14,15))
phenocode1="1"
phenocode2="11"
pvaluequatile=0.008
heatmappdfname=paste(file,"heatmap","pq=",pvaluequatile,"pdf",sep=".")
pvalueSparse(fpcascore,saturationvalue,phen,"1","11",pvaluequatile,heatmappdfname)
}
#=========================================================






###############################################################
# heatmap.3

# mat <- matrix(1:100, byrow=T, nrow=10)
# column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
# column_annotation<-as.matrix(column_annotation)
# colnames(column_annotation) <- "True"
# heatmap.3(mat, ColSideColors=column_annotation)

# mat <- matrix(1:100, byrow=T, nrow=10)
# column_annotation1 <- sample(c("red", "blue", "green"), 10, replace=T)
# column_annotation2 <- sample(c("red", "blue", "green"), 10, replace=T)
# column_annotation <- cbind(column_annotation1,column_annotation2)    
# column_annotation<-as.matrix(column_annotation)
# colnames(column_annotation) <- c("True","False")
# heatmap.3(mat, ColSideColors=column_annotation)
# dev.off()

heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
                             else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
                             else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(colnames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}
