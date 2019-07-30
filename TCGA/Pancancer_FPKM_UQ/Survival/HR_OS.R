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
    thres<-mean(newdata[i,],na.rm=T)
    dat$Rna[dat$Rna<=thres]<-0
    dat$Rna[dat$Rna>thres]<-1
    hr.fit<-summary(coxph(Surv(week,censored)~Rna,dat))
    hr1=hr.fit$coefficients[1,]
    hr2=hr.fit$conf.int[1,]
    HR<-rbind(HR,c(hr1,hr2[3],hr2[4]))
  }
  print(i)
  rownames(HR)<-TCGAProjects
  m<-metagen(HR[,1],seTE=HR[,3],comb.fixed = TRUE,comb.random = TRUE,prediction=F,sm="HR")
  fixedEffect<-c(exp(m$TE.fixed),exp(m$lower.fixed),exp(m$upper.fixed),m$pval.fixed)
  randomEffect<-c(exp(m$TE.random),exp(m$lower.random),exp(m$upper.random),m$pval.random)
  out2<-rbind(out2,c(fixedEffect,randomEffect))
}
colnames(out2)<-c("TE.fixed","lower.fixed","upper.fixed","pval.fixed","TE.random","lower.random","upper.random","pval.random")
rownames(out2)<-row.names(input)[1:4]
out2<-data.frame(out2,ENSG2Symbol(rownames(out2)))
write.table(out2,file="TCGA-HR-OS-Meta-Pvalue-2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)
out2<-out2[order(out2$pval.random,decreasing = F),]
write.table(out2,file="TCGA-HR-OS-Meta-Pvalue.Sort.2019.txt",sep="\t",col.names = NA,row.names = T,quote=F)

