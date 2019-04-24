# install.packages("openxlsx")
library("openxlsx")
library("devtools")
library("GSEABase")
library("ComplexHeatmap")

ci95<-function(x){
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  m<-round(mean(x),2)
  d<-round(mean(x)-error,2)
  u<-round(mean(x)+error,2)
  paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
}


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data[1:3, 1:3]
data[is.na(data)] = ""

####### Most frequent mutation top50 genes
data<-data[head(order(unlist(apply(data,1,function(x) sum(x !=""))),decreasing = T),50),]
input<-as.matrix(data)
head(input)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  Insertion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#008000", col = NA))
  },
  Complex = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "green", col = NA))
  },
  Others = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "purple", col = NA))
  }
)

col = c("SNV" = "red", "Deletion" = "blue", "Insertion" = "#008000","Complex"="green","Others"="purple")
AT=c("SNV", "Deletion", "Insertion","Complex","Others")
oncoPrint(input, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = AT, labels = AT),
          row_names_gp = gpar(fontsize = 10))

oncoPrint(t(input), get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = AT, labels = AT),
          row_names_gp = gpar(fontsize = 10))

####### most freqent mutation genes
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data[1:3, 1:3]
data[is.na(data)] = ""
head(data)
target<-unique(c("TTN","TP53","KRAS","SPOP","MUC16","SYNE1","KMT2D","SPTA1","KMT2C","RYR2","LRP1B","HMCN1","OBSCN","CSMD1","CSMD3","FAT3","ATM","RYR3","CACNA1E","MUC17"))
input<-data[c(na.omit(match(target,rownames(data)))),]
input<-as.matrix(input)
head(input)
oncoPrint(input, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = AT, labels = AT))


####### cancer Driver mutation genes
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data[1:3, 1:3]
data[is.na(data)] = ""
head(data)
target<-unique(c("TP53","SPOP","FOXA1","KMT2D","KMT2C","ATM","PTEN","KDM6A","ZFHX3","FAT4","PIK3CA","CDK12","CTNNB1","ARID2","GRIN2A","APC","RNF213","AKAP9","PTPRC","PTPRC","BRCA2","TMPRSS2","NCOR1",
                 "TPR","SPEN","MED12","MECOM","SETD2","BRAF","KMT2A","TBX3","USP6","AFF3","CNOT3","CHD4","MTOR","ARID1A","ZNF521","NCOR2","NRG1","JAK1","CREBBP","TRRAP","CDH11","FHIT",
                 "KAT6B","IF4","RNF43","RUNX11","HRAS","EP300","BCL11B","MYH9","HNF1A","SF3B1","NUMA1","ERBB4","MAP2K4","TET1","IL6ST","MYH11","NTRK1","MET","NSD1","TRIM33","ERCC2","SMAD4"))
input<-data[c(na.omit(match(target,rownames(data)))),]
input<-as.matrix(input)
head(input)
oncoPrint(input, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = AT, labels = AT),
          row_names_gp = gpar(fontsize = 10))
oncoPrint(t(input), get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = AT, labels = AT),
          row_names_gp = gpar(fontsize = 10))

#################################################################################################################
#################################################################################################################
# install.packages("phangorn")
install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
install.packages("geiger")
install.packages("fs")
library("phangorn")
library("devtools")
library("ape")
text.string<-"(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"
vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE,edge.width=2)
library(phytools)
roundPhylogram(vert.tree)
plot(unroot(vert.tree),type="unrooted",no.margin=TRUE,lab4ut="axial",edge.width=2)

install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")
library(ape)
library(phangorn)
library(seqinr)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data<-t(data.matrix(data))
iid<-unlist(lapply(strsplit(rownames(data),"_"),function(x) x[1]))
par(mfrow=c(4,2))
for(i in unique(iid)){
iiid<-which(iid %in% i)
if(length(iiid)>2){
newdata<-data[iiid,]
newdata[newdata ==1]=0
newdata[newdata >1] = 1
newdata<-rbind(newdata,Germline=0)
input<-as.phyDat(newdata, type="USER", levels = c(0, 1))
pratchet <- pratchet(input)
plot(root(pratchet,outgroup="Germline"),cex=2)
}
}

data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data<-t(data.matrix(data))
newdata<-data
newdata[newdata ==1]=0
newdata[newdata >1] = 1
newdata<-rbind(newdata,Germline=0)
input<-as.phyDat(newdata, type="USER", levels = c(0, 1))
pratchet <- pratchet(input)
plot(root(pratchet,outgroup="Germline"),cex=0.85)

roundPhylogram(root(pratchet,outgroup="Germline"),cex=0.25)
# roundPhylogram(pratchet)
# plot(pratchet)
# plot(unroot(pratchet),type="unrooted",no.margin=TRUE,lab4ut="axial",edge.width=2,cex=1.2)

#################################################################################################################
#################################################################################################################
#install.packages("VennDiagram")
library(gplots)
library(VennDiagram)

count<-read.table("mutation.counts.txt")
count<-count[1:(nrow(count)-1),]
label<-unlist(lapply(strsplit(as.character(count[,2]),"[.]"),function(x) x[1]))
vector<-as.numeric(count[,1])
names(vector)<-label
barplot(vector,col="blue")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data<-data.matrix(data)
sid<-unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1]))
for(id in iid){
inr<-which(sid %in% id)
if(length(inr)==3){
input<-data[,inr]
input[input ==1]=0
input[input>1] = 1
head(input)
input<-input[-which(rowSums(input)==0),]
colnames(input)
T1<-rownames(input)[which(input[,1]>0)]
T2<-rownames(input)[which(input[,2]>0)]
T3<-rownames(input)[which(input[,3]>0)]
vd <- venn.diagram(list(T1=T1, T2=T2, T3=T3),fill = 2:4,filename=NULL)
pdf(paste("Figure",id,"venn.pdf",sep="."))
grid.draw(vd)
dev.off()
}
}

for(id in iid){
  inr<-which(sid %in% id)
  if(length(inr)==2){
    input<-data[,inr]
    input[input ==1]=0
    input[input>1] = 1
    head(input)
    input<-input[-which(rowSums(input)==0),]
    colnames(input)
    dim(input)
    T1<-rownames(input)[which(input[,1]>0)]
    T2<-rownames(input)[which(input[,2]>0)]
    vd <- venn.diagram(list(T1=T1, T2=T2),fill = 2:3,filename=NULL)
    pdf(paste("Figure",id,"venn.pdf",sep="."))
    grid.draw(vd)
    dev.off()
  }
}

for(id in iid){
  inr<-which(sid %in% id)
  if(length(inr)==4){
    input<-data[,inr]
    input[input ==1]=0
    input[input>1] = 1
    head(input)
    input<-input[-which(rowSums(input)==0),]
    colnames(input)
    T1<-rownames(input)[which(input[,1]>0)]
    T2<-rownames(input)[which(input[,2]>0)]
    T3<-rownames(input)[which(input[,3]>0)]
    T4<-rownames(input)[which(input[,4]>0)]
    vd <- venn.diagram(list(T1=T1, T2=T2),fill = 2:5,filename=NULL)
    pdf(paste("Figure",id,"venn.pdf",sep="."))
    grid.draw(vd)
    dev.off()
  }
}

#################################################################################################################
#################################################################################################################
tps<-c()
Rowname<-c()
for(id in unique(iid)){
  print(id)
  inr<-which(sid %in% id)
  if(length(inr)>2){
  input<-data[,inr]
  input[input ==1]=0
  input[input>1] = 1
  head(input)
  input<-input[-which(rowSums(input)==0),]
  colnames(input)
  head(input)
  trunk<-sum(rowSums(input)==3)
  private<-sum(rowSums(input)==1)
  share<-nrow(input)-trunk-private
  tps<-rbind(tps,c(trunk,share,private))
  Rowname<-c(Rowname,id)
  }
}

rownames(tps)<-Rowname
colnames(tps)<-c("trunk","share","private")
TPS<-tps/(rowSums(tps))
write.table(tps,file="trunk.share.private.tps.num.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.table(TPS,file="trunk.share.private.tps.prop.txt",sep="\t",col.names = NA,row.names = T,quote=F)

#################################################################################################################
#################################################################################################################

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data<-data.matrix(data)
sid<-unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1]))
sid
output<-c()
for(id in unique(iid)){
    print(id)
    inr<-which(sid %in% id)
    input<-data[,inr]
    input[input ==1]=0
    input[input>1] = 1
    head(input)
    idata<-unlist(apply(data.frame(input),1,function(x) sum(x)>0))
    output<-cbind(output,idata)
}
colnames(output)<-unique(iid)
output[output=="TRUE"]<-1
output[output=="FALSE"]<-0
head(output)
write.table(output,file="prostate.mutationProfile.txt",sep="\t",col.names = NA,row.names = T,quote=F)



#################################################################################################################
#################################################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")
library("fgsea")
library("ReactomePA")
library("clusterProfiler")
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)

data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
de <-c("4312","8318","10874","55143","55388","991")
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
emapplot(x)
cnetplot(x, categorySize="pvalue", foldChange=geneList)


y <- gsePathway(geneList, nPerm=10000,pvalueCutoff=0.2,pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y)
head(res)
emapplot(y, color="pvalue")
gseaplot(y, geneSetID = "R-HSA-69242")
viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)



#################################################################################################################
#################################################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CancerMutationAnalysis", version = "3.8")

