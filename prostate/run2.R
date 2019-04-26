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


################################################################################
library(gplots)
library(VennDiagram)
count<-read.table("prostat.type.txt",head=T,row.names = 1)
iid<-unlist(lapply(strsplit(as.character(colnames(count)),"[.]"),function(x) x[1]))
data<-data.matrix(count)
head(data)
data[data==1]<-0
data[data>0]<-1
colnames(data)<-label
head(data)
sid<-unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1]))
sid

for(id in unique(sid)){
  inr<-which(sid %in% id)
  if(length(inr)==3){
    input<-data[,inr]
    head(input)
    rmv<-which(rowSums(input)==0)
    length(rmv)/nrow(input)
    input<-input[-which(rowSums(input)==0),]
    colnames(input)
    T1<-rownames(input)[which(input[,1]>0)]
    T2<-rownames(input)[which(input[,2]>0)]
    T3<-rownames(input)[which(input[,3]>0)]
    vd <- venn.diagram(list(T1=T1, T2=T2, T3=T3),fill = 2:4,filename=NULL)
    pdf(paste("Figure",id,"point.venn.pdf",sep="."))
    grid.draw(vd)
    dev.off()
  }
}


for(id in unique(sid)){
  inr<-which(sid %in% id)
  if(length(inr)==2){
    input<-data[,inr]
    head(input)
    rmv<-which(rowSums(input)==0)
    length(rmv)/nrow(input)
    input<-input[-which(rowSums(input)==0),]
    colnames(input)
    T1<-rownames(input)[which(input[,1]>0)]
    T2<-rownames(input)[which(input[,2]>0)]
    vd <- venn.diagram(list(T1=T1, T2=T2),fill = 2:3,filename=NULL)
    pdf(paste("Figure",id,"point.venn.pdf",sep="."))
    grid.draw(vd)
    dev.off()
  }
}


count<-read.table("prostat.type.txt",head=T,row.names = 1)
iid<-unlist(lapply(strsplit(as.character(colnames(count)),"[.]"),function(x) x[1]))
sid<-unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1]))
data<-data.matrix(count)
head(data)
data[data==1]<-0
data[data>0]<-1
colnames(data)<-iid
data<-t(data)
newdata<-rbind(data,Germline=0)
newdata[,1]
input<-as.phyDat(newdata, type="USER", levels = c(0, 1))
pratchet <- pratchet(input)
pdf("Figure_Prostate_Root_Phylogenetic_treemode.SNV.pdf")
plot(root(pratchet,outgroup="Germline"),cex=0.85)
roundPhylogram(root(pratchet,outgroup="Germline"),cex=0.25)
dev.off()


iid<-unlist(lapply(strsplit(rownames(data),"_"),function(x) x[1]))
for(i in unique(iid)){
  iiid<-which(iid %in% i)
  if(length(iiid)>2){
    newdata<-data[iiid,]
    newdata[1:2,1:2]
    newdata<-rbind(newdata,Germline=0)
    input<-as.phyDat(newdata, type="USER", levels = c(0, 1))
    pratchet <- pratchet(input)
    treeRatchetBL<- acctran(pratchet,input)
    write.tree(treeRatchetBL,file=paste(i,"pratchet.tree",sep="."))
    pdf(paste(i,"pratchet.raw.tree.pdf",sep="."))
    plot(root(treeRatchetBL,outgroup="Germline"),cex=0.95)
    dev.off()
    pdf(paste(i,"pratchet.round.tree.pdf",sep="."))
    roundPhylogram(root(treeRatchetBL,outgroup="Germline"),cex=0.25)
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
###########################     SNV based Sum(sub-Sample)       #########################################
#################################################################################################################
sheet7= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 7)
sheet7= as.data.frame(sheet7)
data<-xl2matrix(sheet7)
head(data)
data[data=="_"]<-0
data[data!=0]<-1
head(data)
iid<-unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1]))
iid
head(data)
output1<-c()
output2<-c()
for(id in unique(iid)){
    print(id)
    inr<-which(iid %in% id)
    input<-data[,inr]
    head(data.frame(input))
    idata1<-unlist(apply(data.frame(input),1,function(x) sum(as.numeric(x))>0))
    idata2<-unlist(apply(data.frame(input),1,function(x) sum(as.numeric(x))))
    output1<-cbind(output1,idata1)
    output2<-cbind(output2,idata2)
    
}
colnames(output1)<-unique(iid)
colnames(output2)<-unique(iid)
output1[output1=="TRUE"]<-1
output1[output1=="FALSE"]<-0
head(output1)
head(output2)
write.table(output1,file="prostate.mutationProfile.SNV.01.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.table(output2,file="prostate.mutationProfile.SNV.num.txt",sep="\t",col.names = NA,row.names = T,quote=F)

##########################################################################################################
##########################################################################################################
##########################################################################################################
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/vcf")
data<-read.table("symbol.snv.txt",head=T,sep="\t",row.names = 1)
data<-data.matrix(data)
iid<-unlist(lapply(strsplit(colnames(data),"_"),function(x) x[1]))
iid
output1<-c()
output2<-c()
for(id in unique(iid)){
  print(id)
  inr<-which(iid %in% id)
  input<-data[,inr]
  input[input ==1]=0
  input[input>1] = 1
  head(input)
  idata1<-unlist(apply(data.frame(input),1,function(x) sum(x)>0))
  idata2<-unlist(apply(data.frame(input),1,function(x) sum(x)))
  output1<-cbind(output1,idata1)
  output2<-cbind(output2,idata2)
  
}
colnames(output1)<-unique(iid)
colnames(output2)<-unique(iid)
output1[output1=="TRUE"]<-1
output1[output1=="FALSE"]<-0
head(output1)
head(output2)
write.table(output1,file="prostate.mutationProfile.01.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.table(output2,file="prostate.mutationProfile.num.txt",sep="\t",col.names = NA,row.names = T,quote=F)

data<-output1
data<-t(data)
data<-data[-match(c("P7","P10","P11","P15"),rownames(data)),]
newdata<-rbind(data,Germline=0)
newdata[,1]
input<-as.phyDat(newdata, type="USER", levels = c(0, 1))
pratchet <- pratchet(input)
treeRatchetBL<- acctran(pratchet,input)
write.tree(treeRatchetBL,file=paste(i,"symbol.pratchet.remove.singlesample.tree",sep="."))
pdf(paste(i,"pratchet.raw.rm.single.tree.pdf",sep="."))
plot(root(treeRatchetBL,outgroup="Germline"),cex=0.95)
dev.off()
pdf(paste(i,"pratchet.round.rm.single.sysmboltree.pdf",sep="."))
roundPhylogram(root(treeRatchetBL,outgroup="Germline"),cex=0.25)
dev.off()
#################################################################################################################
#################################################################################################################
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("ReactomePA", version = "3.8")
BiocManager::install("urltools", version = "3.8")
BiocManager::install("backports", version = "3.8")
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("igraph", version = "3.8")
BiocManager::install("fgsea", version = "3.8")
BiocManager::install("DOSE", version = "3.8")
BiocManager::install("ggraph", version = "3.8")
BiocManager::install("enrichplot", version = "3.8")
BiocManager::install("doing", version = "3.8")
BiocManager::install("data.table", version = "3.8")
library("devtools")
library("backports")
library("fgsea")
library("ReactomePA")
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("ReactomePA")

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
library("CancerMutationAnalysis")
data(WoodBreast07)
data(WoodColon07)
data(JonesPancreas08)
data(ParsonsGBM08)
data(ParsonsMB11)

head(GeneAlterGBM)
head(GeneCovGBM)
head(GeneCovGBM)
ScoresGBM <- cma.scores(cma.alter = GeneAlterGBM,
                        cma.cov = GeneCovGBM,
                        cma.samp = GeneSampGBM,
                        passenger.rates = BackRatesGBM["MedianRates",])
head(ScoresGBM)
set.seed(188310)
FdrGBM <-  cma.fdr(cma.alter = GeneAlterGBM,
                   cma.cov = GeneCovGBM,
                   cma.samp = GeneSampGBM,
                   scores = "logLRT",
                   passenger.rates = BackRatesGBM["MedianRates",],
                   showFigure=TRUE,
                   cutoffFdr=0.1,
                   M = 5)

head(FdrGBM[["logLRT"]])

#################################################################################################################
##################### Mapp to COSMIC database Profile ############################################################
#################################################################################################################
devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)





#################################################################################################################
#################################################################################################################
library("readxl")
xl2matrix<-function(x){
  rownames(x) <-x[,1]
  x<-x[,2:ncol(x)]
  x
}
sheet1= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 1)
sheet7= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 7)
sheet9= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 9)
sheet8= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 8)
sheet10= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 10)

sheet1= as.data.frame(sheet1)
sheet7= as.data.frame(sheet7)
sheet8= as.data.frame(sheet8)
sheet9= as.data.frame(sheet9)
sheet10= as.data.frame(sheet10)

sheet1<-xl2matrix(sheet1)
sheet7<-xl2matrix(sheet7)
sheet9<-xl2matrix(sheet9)
sheet8<-xl2matrix(sheet8)
sheet10<-xl2matrix(sheet10)

head(sheet1)
head(sheet7)
head(sheet9)
head(sheet8)
head(sheet10)
sheet7[sheet7=="_"]<-0
sheet7[sheet7!=0]<-1
############################################
Rlt<-c()
head(sheet8)
for(i in 1:nrow(sheet8)){
  rlt<-summary(glm(as.numeric(sheet8[i,])~sheet1$Age))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet8)
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(lm(colSums(sheet8)~sheet1$Age))
write.table(Rlt,file="sheet7.SNV.vs.age.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
## Mutation with fPSA/tPSA
Rlt<-c()
for(i in 1:nrow(sheet8)){
  rlt<-summary(glm(as.numeric(sheet8[i,])~sheet1$fPSA_tPSA))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet8)
hist(log(Rlt[,4],10))
head(Rlt[order(Rlt[,4]),])
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(lm(colSums(sheet8)~sheet1$fPSA_tPSA))
write.table(Rlt,file="sheet7.SNV.vs.fpsa.tpsa.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
## Mutation with fPSA/tPSA
Rlt<-c()
for(i in 1:nrow(sheet8)){
  rlt<-summary(glm(as.numeric(sheet8[i,])~sheet1$Gleason_Score_Value))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet8)
hist(log(Rlt[,4],10))
head(Rlt[order(Rlt[,4]),])
sum(Rlt[,1]>0)
sum(Rlt[,1]<0)
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(lm(colSums(sheet8)~sheet1$Gleason_Score_Value))
plot(colSums(sheet8)~sheet1$Gleason_Score_Value)
t.test(colSums(sheet8)[which(sheet1$Gleason_Score_Value>7.5)],colSums(sheet8)[which(sheet1$Gleason_Score_Value<7.5)])
write.table(Rlt,file="sheet7.SNV.vs.Gleason_Score.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
## Mutation with fPSA/tPSA
Rlt<-c()
for(i in 1:nrow(sheet8)){
  rlt<-summary(glm(as.numeric(sheet8[i,])~sheet1$Metastasis_binary))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet8)
hist(log(Rlt[,4],10))
head(Rlt[order(Rlt[,4]),])
sum(Rlt[,1]>0)
sum(Rlt[,1]<0)
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(glm(colSums(sheet8)~sheet1$Metastasis_binary))
plot(colSums(sheet8)~sheet1$Gleason_Score_Value)
t.test(colSums(sheet8)[which(sheet1$Gleason_Score_Value>7.5)],colSums(sheet8)[which(sheet1$Gleason_Score_Value<7.5)])
write.table(Rlt,file="sheet7.SNV.vs.Metastasis.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)


Rlt<-c()
head(sheet10)
for(i in 1:nrow(sheet10)){
  rlt<-summary(glm(as.numeric(sheet10[i,])~sheet1$Age))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet10)
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(lm(colSums(sheet10)~sheet1$Age))
write.table(Rlt,file="sheet10.Symbol.vs.age.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
## Mutation with fPSA/tPSA
Rlt<-c()
for(i in 1:nrow(sheet10)){
  rlt<-summary(glm(as.numeric(sheet10[i,])~sheet1$fPSA_tPSA))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet10)
hist(log(Rlt[,4],10))
head(Rlt[order(Rlt[,4]),])
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(lm(colSums(sheet10)~sheet1$fPSA_tPSA))
write.table(Rlt,file="sheet10.Symbol.vs.fpsa.tpsa.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
## Mutation with fPSA/tPSA
Rlt<-c()
for(i in 1:nrow(sheet10)){
  rlt<-summary(glm(as.numeric(sheet10[i,])~sheet1$Gleason_Score_Value))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet10)
hist(log(Rlt[,4],10))
head(Rlt[order(Rlt[,4]),])
sum(Rlt[,1]>0)
sum(Rlt[,1]<0)
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(lm(colSums(sheet10)~sheet1$Gleason_Score_Value))
plot(colSums(sheet10)~sheet1$Gleason_Score_Value)
t.test(colSums(sheet10)[which(sheet1$Gleason_Score_Value>7.5)],colSums(sheet10)[which(sheet1$Gleason_Score_Value<7.5)])
write.table(Rlt,file="sheet10.Symbol.vs.Gleason_Score.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)
## Mutation with fPSA/tPSA
Rlt<-c()
for(i in 1:nrow(sheet10)){
  rlt<-summary(glm(as.numeric(sheet10[i,])~sheet1$Metastasis_binary))$coefficients[2,]
  Rlt<-rbind(Rlt,rlt)
}
rownames(Rlt)<-rownames(sheet10)
hist(log(Rlt[,4],10))
head(Rlt[order(Rlt[,4]),])
sum(Rlt[,1]>0)
sum(Rlt[,1]<0)
sum(Rlt[,1]>0)/sum(Rlt[,1]<0)
summary(glm(colSums(sheet10)~sheet1$Metastasis_binary))
plot(colSums(sheet10)~sheet1$Gleason_Score_Value)
t.test(colSums(sheet10)[which(sheet1$Gleason_Score_Value>7.5)],colSums(sheet10)[which(sheet1$Gleason_Score_Value<7.5)])
write.table(Rlt,file="sheet10.Symbol.vs.Metastasis.pvalue.txt",sep="\t",quote=F,col.names = NA,row.names = T)

data<-sheet8
data<-t(data)
newdata<-rbind(data,Germline=0)
newdata[,1]
input<-as.phyDat(newdata, type="USER", levels = c(0, 1))
pratchet <- pratchet(input)
treeRatchetBL<- acctran(pratchet,input)
write.tree(treeRatchetBL,file=paste(i,"pratchet.tree",sep="."))
pdf(paste(i,"pratchet.raw.tree.pdf",sep="."))
plot(root(treeRatchetBL,outgroup="Germline"),cex=0.95)
dev.off()
pdf(paste(i,"pratchet.round.tree.pdf",sep="."))
roundPhylogram(root(treeRatchetBL,outgroup="Germline"),cex=0.25)
dev.off()

pratchet <- pratchet(input)
pdf("Figure_Prostate_Root_Phylogenetic_treemode.SNV.pdf")
plot(root(pratchet,outgroup="Germline"),cex=0.85)
roundPhylogram(root(pratchet,outgroup="Germline"),cex=0.25)
dev.off()

