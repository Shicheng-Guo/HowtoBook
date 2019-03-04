HeatMap<-function(data){
  
  #  note: this function include correlation based heatmap (pearson or spearman)
  #  data: row is gene and column is sample
  #  colname and rowname cannot be NULL  
  #  Usage example:
  #  test<- matrix(runif(100),nrow=20)
  #  colnames(test)=c("A","A","A","B","B")
  #  rownames(test)=paste("Gene",1:20,sep="")
  #  HeatMap(test)
  
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  
  Hclust=function(x){hclust(x,method="complete")}
  Distfun=function(x){as.dist((1 - cor(t(x),method = "spearman")))}
  
  ColSideColors=sidecol(colnames(data))
  
  heatmap.2(data,trace="none", 
            hclust=Hclust,
            distfun=Distfun, 
            cexRow = 0.5, cexCol = 0.5,
            ColSideColors=ColSideColors,
            density.info="none",col=colors,
            Colv=T,Rowv = TRUE,
            keysize=0.9, margins = c(5, 10)
  )
}

test<- matrix(runif(100),nrow=20)
colnames(test)=c("A","A","A","B","B")
rownames(test)=paste("Gene",1:20,sep="")
HeatMap(test)
