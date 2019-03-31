setwd("/mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/GEO")

library("metafor")
library("GEOquery")

GPL96<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/GEO/GPL96.anno.txt",sep="\t",head=T)

GSE55457 <- getGEO("GSE55457")
data1 <- as.data.frame(exprs(GSE55457[[1]]))
phen1 <- pData(phenoData(GSE55457[[1]]))

data <- getGEO("GSE55584")
data2 <- as.data.frame(exprs(data[[1]]))
phen2 <- pData(phenoData(data[[1]]))

data <- getGEO("GSE55235")
data3 <- as.data.frame(exprs(data[[1]]))
phen3 <- pData(phenoData(data[[1]]))

newphen1<-phen1$'clinical status:ch1'
newphen2<-phen2$'clinical status:ch1'
newphen3<-phen3$'disease state:ch1'

newphen1<-gsub("normal control","GSE55457_Normal",newphen1)
newphen1<-gsub("rheumatoid arthritis","GSE55457_RA",newphen1)
newphen1<-gsub("osteoarthritis","GSE55457_Normal",newphen1)

newphen2<-gsub("osteoarthritis","GSE55584_Normal",newphen2)
newphen2<-gsub("rheumatoid arthritis","GSE55584_RA",newphen2)

newphen3<-gsub("healthy control","GSE55235_Normal",newphen3)
newphen3<-gsub("osteoarthritis","GSE55235_Normal",newphen3)
newphen3<-gsub("synovial tissue isolated from osteoarthritic joint","GSE55235_Normal",newphen3)
newphen3<-gsub("rheumatoid arthritis","GSE55584_RA",newphen3)

input<-data.frame(data1,data2,data3)
Seq<-c(newphen1,newphen2,newphen3)

Symbol<-GPL96[match(rownames(input),GPL96$ID),2]

P<-c()
beta<-c()
for(i in 1:nrow(input)){
  print(rownames(input)[i])
  mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
  sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
  num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
  m1i=mean[seq(1,length(mean),by=2)]
  m2i=mean[seq(2,length(mean),by=2)]
  sd1i=sd[seq(1,length(mean),by=2)]
  sd2i=sd[seq(2,length(mean),by=2)]
  n1i=num[seq(1,length(mean),by=2)]
  n2i=num[seq(2,length(mean),by=2)]
  Source<-unlist(lapply(strsplit(names(m1i),"_"),function(x) x[1]))
  output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
  output$source=Source
  output<-na.omit(output)
  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
  res <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
  P<-c(P,res$pval)
  beta<-c(beta,res$beta)
  filename=gsub("/","_",Symbol[i])
  if(res$pval<0.00000000005){
    pdf(paste(filename,".pdf",sep=""))
    plot(res)
    text(0, -0.1, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (Q = ",
                                                .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                                ", p = ", .(formatC(res$QEp, digits=2, format="E")), "; ", I^2, " = ",
                                                .(formatC(res$I2, digits=1, format="f")), "%)")))
    
    text(0, -0.25, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (beta = ",
                                                 .(formatC(res$beta, digits=2, format="f")), ", se = ", .(formatC(res$se, digits=2, format="f")),
                                                 ", zval = ", .(formatC(res$zval, digits=2, format="f")), "; ", P, " = ",
                                                 .(formatC(res$pval, digits=2, format="E")), ")")))
    text(0, -0.4, pos=4, cex=0.75, rownames(input)[i])
    dev.off()
    system("mv *.pdf ./meta")
  }
}







