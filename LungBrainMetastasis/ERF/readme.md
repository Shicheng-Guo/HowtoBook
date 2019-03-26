```
library("survival")
library("survminer")
library("meta")
setwd("~/hpc/project/LungBrainMetastasis/submission")
gene="ERF"
erf<-read.table("ERF.tsv",head=T,sep="\t")
erf$censored<-abs(as.numeric(erf$censored)-2)
erf$label<-abs(as.numeric(erf$label)-1)
head(erf)
coxph(Surv(time,censored)~label,erf)
hr<-summary(coxph(Surv(time,censored)~label,erf))$coefficients[1,]
fit <- survfit(Surv(time,censored)~label, data = erf)
survp<-ggsurvplot(fit, data = erf,conf.int = F,pval = TRUE,
                  fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                  palette = c("#E7B800","#2E9FDF"),
                  legend = "bottom",legend.title = paste(gene,"HR=",round(hr[2],2),sep=" "),
                  legend.labs = c("Wile","Mutation"))
ggsave(file = paste("./",gene,".Mutation.lungBrainMetastasis.pdf",sep=""), survp$plot)
```
