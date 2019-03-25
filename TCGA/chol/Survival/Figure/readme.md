```
library("survival")
library("survminer")

OS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/OverallSurvivalTime.txt",head=T,sep="\t")
data<-input[,which(id2bin(colnames(input))==1)]
newdata<-data[,na.omit(match(OS$submitter_id,id2phen3(colnames(data))))]
colnames(newdata)<-id2phen3(colnames(newdata))

phen<-OS[match(colnames(newdata),OS$submitter_id),]
head(phen)
phen$censored<-as.numeric(! phen$censored)
phen$week=phen$time/7
head(phen)

ENST2Symbol<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/ENSG.ENST.ENSP.Symbol.hg19.bed")
ENSG<-unlist(lapply(strsplit(as.character(rownames(newdata)),split="[.]"),function(x) x[1]))
Symbol<-ENST2Symbol[match(ENSG,ENST2Symbol$V7),5]

HR<-c()
for(i in 1:nrow(newdata)){
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
  HR<-rbind(HR,hr)
  print(i)
}
HR1<-HR

HR<-HR1
NewRowName<-paste(Symbol,rownames(newdata),sep="_")
HRV=cbind(HR,Symbol,rownames(newdata))
write.table(HRV,file="../../RNAseq_Suvival.HR.txt",sep="\t",quote=F,row.names = T,col.names = NA)


HRVV<-read.table(file="../../RNAseq_Suvival.HR.txt",sep="\t",head=T)

DATA<-HRVV[which(as.numeric(HRVV[,6])<10^-50),]

for(j in 1:nrow(DATA)){
gene=DATA[j,7]
i=match(gene,Symbol)
dat<-data.frame(Rna=newdata[i,],phen)
dat$Rna[dat$Rna<=median(dat$Rna)]<-0
dat$Rna[dat$Rna>median(dat$Rna)]<-1
coxph(Surv(week,censored)~Rna,dat)
hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
fit <- survfit(Surv(week,censored)~Rna, data = dat)
survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
                  fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                  palette = c("#E7B800","#2E9FDF"),
                  legend = "bottom",legend.title = paste(gene,"HR=",round(hr[2],2),sep=" "),
                  legend.labs = c("Low-Expression","High-Expression"))
ggsave(file = paste("../../Survival_Figure/",gene,".Pancancer.pdf",sep=""), survp$plot)
}

CHOL<-c("ADCY2","AHRR","ARRB2","BCAN","DLX5","DMRTA2","FGF8","GRIN1","GSC","HIST1H3G","HOXA1","HOXA9","HOXD12","HOXD4","IFFO1","LHX1","LHX8","LINC00461","MNX1","NR2E1","OTX1","PEG3","PITX2","POU4F3","RBMY1F","SATB2-AS1","SHOX2","SIM1","SLC2A14","SPEG","TFAP2E","TRH","TRIM58","UMOD","ZNF518B")
ESCA2<-c("ACSS3","AKAP12","BCAT1","CD200","CDO1","CYP2A13","DMRTA21","DMRTA22","FEZF2","FNDC1","FOXB1","FOXE1","GHSR","GRIK4","GSC","GSX1","HOXD9","ITGA8","LINC00466","LINE-1","LncRNA","MIR129-2","NKAIN4","NKX1-2","NOTO","NR2E1","NTM","PCDH10","PDE4B","PEX5L","RBMY1F","RUNDC3B","RXFP3","RYR2","SHOX2","SLC6A2","SOX11","SPAG6","SRCIN1","THY1","TLX1","TRH","VSTM2A","WNT7B","ZNF273","ZNF385D","ZNF415","ZNF583","ZNF781")
ESCA1<-c("TRIM71","USP44","SCOC","HMX3","ZNF878","ZNF345","TACC2","ADHFE1","TMEM132C","ZNF385B","GFRA1","DPY19L2P4","cg15830431","KIAA0226L","cg20655070","ZNF71","KCNA6","SLITRK5","KCNA3","LRAT","ELMO1","ZNF793","ZNF542","ZNF844","TLX2","cg05249644","CH25H","ZNF418","ZNF570","TBX4","SALL1","SPATA32","EOMES","ZNF790","ZIK1","LINE-1","cg19396867","APC","chrM","CNR1","ZNF132","EYA4","ZNF461","CHST2","ZNF470","ZNF569","TFPI2","RNLS")
FULL<-c(CHOL,ESCA1,ESCA2)

for(j in FULL){
  gene=j
  i=match(gene,Symbol)
  dat<-data.frame(Rna=newdata[i,],phen)
  dat$Rna[dat$Rna<=median(dat$Rna)]<-0
  dat$Rna[dat$Rna>median(dat$Rna)]<-1
  coxph(Surv(week,censored)~Rna,dat)
  hr<-summary(coxph(Surv(week,censored)~Rna,dat))$coefficients[1,]
  fit <- survfit(Surv(week,censored)~Rna, data = dat)
  survp<-ggsurvplot(fit, data = dat,conf.int = F,pval = TRUE,
                    fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                    palette = c("#E7B800","#2E9FDF"),
                    legend = "bottom",legend.title = paste(gene,"HR=",round(hr[2],2),sep=" "),
                    legend.labs = c("Low-Expression","High-Expression"))
  ggsave(file = paste("../../Survival_Figure/FULL/",gene,".Pancancer.pdf",sep=""), survp$plot)
}


```
