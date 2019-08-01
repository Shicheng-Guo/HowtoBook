source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/GscTools.R")
library("meta")
library("metafor")
library("survival")
library("survminer")

# i=grep("ENSG00000266601",rownames(input))
load("~/hpc/methylation/Pancancer/RNA-seq/rnaseqdata.pancancer.RData")
TCGAProjects=c("LIHC","CHOL")
panc<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/panc.txt",head=T)
phen1=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/TCGA-clinical-11093.tsv",header = T,sep="\t")
phen2=read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/extdata/File_metadata2.txt",header = T,sep="\t")
head(phen1)
head(phen2)
colnames(rnaseqdata)<-unlist(lapply(strsplit(colnames(rnaseqdata),"/"),function(x) x[2]))
phen<-data.frame(phen2,phen1[match(phen2$cases.0.case_id,phen1$case_id),])
phen$file_name=gsub(".gz","",phen$file_name)
# prepare phenotype information
phen<-phen[match(colnames(rnaseqdata),phen$file_name),]
phen$phen4<-id2phen4(phen$cases.0.samples.0.submitter_id)
phen$phen3<-id2phen3(phen$cases.0.samples.0.submitter_id)
phen$phen2<-id2bin(phen$cases.0.samples.0.submitter_id)
phen$pid<-phen$project_id
head(phen)

OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
# match survival information
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

out2<-c()
for(i in 1:nrow(input)){
  HR<-c()
  for(TCGAProject in TCGAProjects){
    newdata<-input[,phen$project_id==paste("TCGA-",TCGAProject,sep="")]
    xphen<-phen[phen$project_id==paste("TCGA-",TCGAProject,sep=""),]
    dat<-data.frame(Rna=newdata[i,],xphen)
    thres<-mean(dat[,1],na.rm=T)
    dat$Rna[dat$Rna<=thres]<-0
    dat$Rna[dat$Rna>thres]<-1
    hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
    hr1=hr.fit$coefficients[1,]
    hr2=hr.fit$conf.int[1,]
    HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
    
    # fit <- survfit(Surv(week,censored)~Rna, data = dat)
    # survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
    #                   fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
    #                   palette = c("#2E9FDF","#E7B800"),
    #                   legend = "bottom",legend.title = rownames(input)[i],
    #                   legend.labs = c("Low-expression","High-expression"))
    # ggsave(file = paste(rownames(input)[i],"_",TCGAProject,"_KM.pdf",sep=""), survp$plot)
  }
  print(i)
  rownames(HR)<-TCGAProjects
  m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
  fixedEffect<-c(exp(m$TE.fixed),exp(m$lower.fixed),exp(m$upper.fixed),m$pval.fixed)
  randomEffect<-c(exp(m$TE.random),exp(m$lower.random),exp(m$upper.random),m$pval.random)
  out2<-rbind(out2,c(fixedEffect,randomEffect))
}

# pdf(paste(rownames(input)[i],".OS.HR.PANC.pdf",sep=""))
# print(rownames(input)[i])
# forest(m,leftlabs = rownames(HR),
#        lab.e = "Intervention",
#        pooled.totals = FALSE,
#        smlab = "",studlab=rownames(HR),
#        text.random = "Overall effect",
#        print.tau2 = FALSE,
#        col.diamond = "blue",
#        col.diamond.lines = "black",
#        col.predict = "red",
#        print.I2.ci = TRUE,
#        digits.sd = 2,fontsize=9,xlim=c(0.2,3))
# dev.off()

colnames(out2)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
rownames(out2)<-row.names(input)
out2<-data.frame(out2,ENSG2Symbol(rownames(out2)))
write.table(out2,file="TCGA-HR-OS-LIHC_CHOL-Pvalue-2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)
out2<-out2[order(out2$pval.random,decreasing = F),]
write.table(out2,file="TCGA-HR-OS-LIHC_CHOL-Pvalue.Sort.2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)

