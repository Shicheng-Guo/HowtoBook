

*(2019-08-01) I re-selected the most interesting proteins to be detected in LIHC samples with 1) over expression in cancer and worse outcome with over expression 2) low expression in cancer (D>0) and better-outcome when over-expression (HR<1)
```
data1<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lihc/protein_array/TCGA-HR-OS-Meta-Pvalue-2019.txt",head=T,sep="\t",row.names=1)
data2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lihc/protein_array/TCGA-HR-OS-LIHC_CHOL-Pvalue-2019.txt",head=T,sep="\t",row.names=1)
data3<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lihc/protein_array/TCGA-SMD-DGE-Meta-Pvalue-2019.txt",head=T,sep="\t",row.names=1)
head(data1)
head(data2)
head(data3)
X1<-data1[match(rownames(data3),rownames(data1)),]
X2<-data2[match(rownames(data3),rownames(data2)),]
xdata<-data.frame(data3,X1,X2)
head(xdata)
RLT<-subset(xdata,xdata[,8]<10^-3 & xdata[,17]<10^-3 &  xdata[,26]<10^-3 & xdata[,1]*(xdata[,10]-1)>0)
dim(RLT)
write.table(RLT,file="Minghua_ProteinArray_P2_D20190801.txt",sep="\t",col.names = NA,row.names = T,quote=F)
getwd()
target<-c("CENPA","TPX2","ZWINT","HJURP","CDCA5","BUB1","BUB1B")
new<-xdata[match(target,xdata[,27]),]
```
