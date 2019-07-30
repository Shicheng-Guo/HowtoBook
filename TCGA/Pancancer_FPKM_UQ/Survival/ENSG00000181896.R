library("meta")
library("metafor")
library("survival")
library("survminer")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")

load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")
TCGAProjects=c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
panc<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/panc.txt",head=T)
phen1=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)

i=grep("ENSG00000181896",rownames(input))

idx<-which(phen$phen2==1 | phen$phen2==11)
phen<-phen[idx,]
input<-rnaseqdata[,idx]
idx<-which(phen$pid %in% paste("TCGA-",TCGAProjects,sep=""))
phen<-phen[idx,]
input<-input[,idx]
input[1:5,1:5]
input<-log(input+1,2)
input<-RawNARemove(input)
input<-RawZeroRemove(input)
Seq<-paste(phen$project_id,phen$phen2,sep="-")
rlt<-c()
coll<-c()
out<-c()
print(i)
mean<-tapply(as.numeric(input[i,]),Seq,function(x) mean(x,na.rm=T))
sd<-tapply(as.numeric(input[i,]),Seq,function(x) sd(x,na.rm=T))
num<-tapply(as.numeric(input[i,]),Seq,function(x) length(x))
m1i=mean[seq(1,length(mean),by=2)]
m2i=mean[seq(2,length(mean),by=2)]
sd1i=sd[seq(1,length(mean),by=2)]
sd2i=sd[seq(2,length(mean),by=2)]
n1i=num[seq(1,length(mean),by=2)]
n2i=num[seq(2,length(mean),by=2)]
Source<-unlist(lapply(strsplit(names(m1i),"-"),function(x) x[2]))
output<-data.frame(cbind(n1i,m1i,sd1i,n2i,m2i,sd2i))
output$source=Source
output<-na.omit(output)
es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output)
rlt<-rbind(rlt,c(i,md$beta,md$pval,md$ci.lb,md$ci.ub,md$I2,md$tau2))
studlab=unlist(lapply(rownames(output),function(x) unlist(strsplit(x,"-"))[2]))
coll<-c(coll,i)
m<-metagen(yi,seTE=vi,data = es,comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="SMD")
print(rownames(input)[i])
pdf(paste(rownames(input)[i],".SMD.PANC.pdf",sep=""))
forest(m,leftlabs =studlab,
       lab.e = "Intervention",
       pooled.totals = FALSE,
       smlab = "",studlab=studlab,
       text.random = "Overall effect",
       print.tau2 = FALSE,
       col.diamond = "blue",
       col.diamond.lines = "black",
       col.predict = "red",
       print.I2.ci = TRUE,
       digits.sd = 2,fontsize=8,xlim=c(-6,1))
dev.off()

fixedEffect<-c(m$TE.fixed,m$lower.fixed,m$upper.fixed,m$pval.fixed)
randomEffect<-c(m$TE.random,m$lower.random,m$upper.random,m$pval.random)
out<-rbind(out,c(fixedEffect,randomEffect))


OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
idx<-which(c(phen$phen2==1))
phen<-phen[idx,]
input<-rnaseqdata[,idx]
input[1:5,1:5]
idx<-na.omit(match(OS$submitter_id,phen$phen3))
input<-input[,idx]
phen<-phen[idx,]
phen<-data.frame(phen,OS[match(phen$phen3,OS$submitter_id),])
phen$censored<-as.numeric(!phen$censored)
phen$week=phen$time/7
Z<-c()
for(z in quantile(input[i,],seq(0, 1, 0.1))[2:9]){
  HR<-c()
  for(TCGAProject in TCGAProjects){
    newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
    xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
    dat<-data.frame(Rna=newdata[i,],xphen)
    dat$Rna[dat$Rna<=z]<-0
    dat$Rna[dat$Rna>z]<-1
    hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
    hr1=hr.fit$coefficients[1,]
    hr2=hr.fit$conf.int[1,]
    HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
  }
  HR<-na.omit(HR)
  HR <- HR[!is.infinite(rowSums(HR[,6:7])),]
  P<-t.test((1-HR[,6])*(1-HR[,7]))$p.value
  print(c(z,P))
  m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
  Z<-rbind(Z,c(z,m$pval.fixed))
  print(c(z,m$pval.fixed))
}
  
thres<-Z[which.min(Z[,2]),1]
HR<-c()
for(TCGAProject in TCGAProjects){
  newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
  xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
  dat<-data.frame(Rna=newdata[i,],xphen)
  dat$Rna[dat$Rna<=mean(newdata[i,],na.rm=T)]<-0
  dat$Rna[dat$Rna>mean(newdata[i,],na.rm=T)]<-1
  hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
  hr1=hr.fit$coefficients[1,]
  hr2=hr.fit$conf.int[1,]
  HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
  if(subpanel==TRUE){
    fit <- survfit(Surv(week,censored)~Rna, data = dat)
    survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
                      fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                      palette = c("#E7B800","#2E9FDF"),
                      legend = "bottom",legend.title = rownames(input)[i],
                      legend.labs = c("Low-expression","High-expression"))
    ggsave(file = paste(rownames(input)[i],"_",TCGAProject,"_KM.pdf",sep=""), survp$plot)
  }
}
print(i)
rownames(HR)<-TCGAProjects
m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
pdf(paste(rownames(input)[i],".OS.HR.PANC.pdf",sep=""))
  print(rownames(input)[i])
  forest(m,leftlabs = rownames(HR),
       lab.e = "Intervention",
       pooled.totals = FALSE,
       smlab = "",studlab=rownames(HR),
       text.random = "Overall effect",
       print.tau2 = FALSE,
       col.diamond = "blue",
       col.diamond.lines = "black",
       col.predict = "red",
       print.I2.ci = TRUE,
       digits.sd = 2,fontsize=9)
dev.off()

fixedEffect<-c(exp(m$TE.fixed),exp(m$lower.fixed),exp(m$upper.fixed),m$pval.fixed)
randomEffect<-c(exp(m$TE.random),exp(m$lower.random),exp(m$upper.random),m$pval.random)
out<-rbind(out,c(fixedEffect,randomEffect))
colnames(out)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
rownames(out)<-c("SMD","HR")

