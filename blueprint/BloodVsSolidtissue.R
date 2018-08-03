
 awk -F"\t" '{print $1"\t"$3}' S.table.txt > blueprint.id.txt
 sam<-read.table("blueprint.id.txt",sep="\t")
 mitometh<-read.table("matrix.txt",head=T,row.names=1,as.is=T,check.names=F)
 blood<-which(!is.na(match(unlist(strsplit(colnames(mitometh),split=".IRF5.tab")),sam[,1])))
 bloodata<-mitometh[1,blood]
 quantile(bloodata,seq(0,1,0.01),na.rm=T)
 exclude=colnames(bloodata)[which(bloodata[1,]>0.1)]
 sam[match(unlist(strsplit(exclude,split=".IRF5.tab")),sam[,1]),]
  
 colnames(mitometh)<-unlist(lapply(colnames(mitometh),function(x) substr(x,1,6)))
 rownames(mitometh)<-lapply(rownames(mitometh),function(x) unlist(strsplit(x,"[:-]"))[2])
 boxplot(t(mitometh),outline = F,horizontal = T)
 par(las=2,cex=0.6,mar=c(6,6,1,4))
 boxplot(t(mitometh),outline = T,horizontal = T)
 par(las=2,cex=0.6,mar=c(6,6,1,4))
 boxplot((mitometh),outline = T,horizontal = T)
