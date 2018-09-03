TtestPValue<-function(data,x1,x2,pair=FALSE){
  data<-data.matrix(data)
  noise<-matrix(rnorm(nrow(data)*ncol(data),0,0.00001),nrow(data),ncol(data))
  data=data+noise

  output<-matrix(NA,dim(data)[1],6)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    print(i)
    out<-data.frame()
    if(pair==TRUE){
      Valid<-nrow(na.omit(data.frame(data[i,x1],data[i,x2])))
    }else{
      Valid<-100
    }
    if( sum(!is.na(data[i,x1]))>=3 & sum(!is.na(data[i,x2]))>=3 & Valid>3){
      tmp1<-try(t.test((data[i,x1]),(data[i,x2]),paired=F, na.action=na.omit))
      output[i,1]<-format(tmp1$p.value, scientific=TRUE)
      output[i,2]<-round(mean(data[i,x1],na.rm=T)-mean(data[i,x2],na.rm=T),3)
      output[i,3]<-round(mean(data[i,x1],na.rm=T),3)
      output[i,4]<-round(mean(data[i,x2],na.rm=T),3)
      output[i,5]<-round(sd(data[i,x1],na.rm=T),3)
      output[i,6]<-round(sd(data[i,x2],na.rm=T),3)
    }
  }
  rownames(output)<-rownames(data)
  output
}

args = commandArgs(trailingOnly=TRUE)
sam<-read.table("~/hpc/epimarker/bedgraph/SampleMatrix.txt",head=T,sep="\t")
data=read.table(args[1],head=T,sep="\t",row.names=1,check.names=F)
chr=unlist(strsplit(args[1],split=".tab.matrix.rlt"))
filename=unlist(lapply(strsplit(colnames(data),split="[.]"),function(x) x[1]))
newsam=sam[match(filename,sam[,1]),]
blood=which(newsam[,3]=="blood")
solid=which(newsam[,3]=="tissue")
rlt<-TtestPValue(data,blood,solid)
output1=paste(args[1],".pvalue.rlt",sep="")
output2=paste(args[1],"Sig.pvalue.rlt",sep="")
colnames(rlt)=c("Pvalue","delta","M1","M2","SD1","SD2")
Rlt<-data.frame(rlt)
write.table(Rlt,file=output1,sep="\t",quote=F,col.names=NA,row.names=T)

bres<-subset(Rlt,as.numeric(as.character(M1))<0.3 & as.numeric(as.character(M2))>0.6 & as.numeric(as.character(Pvalue))<0.0000001)
sres<-subset(Rlt,as.numeric(as.character(M1))>0.6 & as.numeric(as.character(M2))<0.3 & as.numeric(as.character(Pvalue))<0.0000001)
Sig<-rbind(bres,sres)
write.table(Sig,file=output2,sep="\t",quote=F,col.names=NA,row.names=T)
