data=read.table("Nash.mf.matrix.txt",head=T,row.names=1,check.names=F)
colnames(data)=unlist(lapply(strsplit(colnames(data),split="_"),function(x) paste(x[1],x[2],sep="_")))
colnames(data)
threshold=30
pdata1=data[,c(2:4,18)]
pdata1[,4]=pdata1[,4]*100
head(pdata1)
F1=which(rowMeans(pdata1)>threshold)
length(F1)
rownames(F1)
sum((data[F1,c(31)])>threshold,na.rm=T)/length(na.omit((data[F1,c(31)])>threshold))
sum((data[F1,c(33)])>threshold,na.rm=T)/length(na.omit((data[F1,c(33)])>threshold))
sum((data[F1,c(35)])>threshold,na.rm=T)/length(na.omit((data[F1,c(35)])>threshold))
sum((data[F1,c(1)])>threshold,na.rm=T)/length(na.omit((data[F1,c(1)])>threshold))
sum(rowMeans(data[F1,c(31,33,35)])>threshold,na.rm=T)/length(na.omit(rowMeans(data[F1,c(31,33,35)])>threshold))

pdata2=data[F1,30:35]
newdata=pdata2[which(pdata2[,6]<pdata2[,4] & pdata2[,4]<pdata2[,2] & pdata2[,2]>pdata2[,1] & pdata2[,4]>pdata2[,3] & pdata2[,6]>pdata2[,5]),]
colnames(newdata)=unlist(lapply(strsplit(colnames(newdata),split="_"),function(x) x[2]))
dim(newdata)
target=subset(newdata,newdata[,2]>60)
input=data[match(rownames(target),rownames(data)),]
colnames(input)=unlist(lapply(strsplit(colnames(input),split="_"),function(x) paste(x[1],x[2],sep="_")))


# LINE-1
data=read.table("LINE-1.mf.txt",head=T,row.names=1,check.names=F)
LINE=colMeans(data,na.rm=T)
names(LINE)=unlist(lapply(strsplit(colnames(data),split="_"),function(x) paste(x[1],x[2],sep="_")))
t.test(LINE[2:4],LINE[c(31,33,35)])

colnames(data)
threshold=30
pdata1=data[,c(2:4,18)]
pdata1[,4]=pdata1[,4]*100
head(pdata1)
F1=which(rowMeans(pdata1)>threshold)
length(F1)
rownames(F1)
sum((data[F1,c(31)])>threshold,na.rm=T)/length(na.omit((data[F1,c(31)])>threshold))
sum((data[F1,c(33)])>threshold,na.rm=T)/length(na.omit((data[F1,c(33)])>threshold))
sum((data[F1,c(35)])>threshold,na.rm=T)/length(na.omit((data[F1,c(35)])>threshold))
sum((data[F1,c(1)])>threshold,na.rm=T)/length(na.omit((data[F1,c(1)])>threshold))
