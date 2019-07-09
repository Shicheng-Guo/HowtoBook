install.packages("survminer")
library("survival")
library("survminer")
library("openxlsx")
library(survminer)
require("survival")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis")
data<-read.xlsx("Result.xlsx",sheet=1,rowNames=T)
P<-c()
for(i in 1:13){
xdata<-data[,c(2,3,4,6,7,8,9,12,13,15)]
xdata<-xdata[-i,]
model <- coxph( Surv(OS,status) ~ .,data =xdata)
ggforest(model,data=xdata)
fit<-summary(model)
P<-c(P,fit$coefficients[8,5])
}
i=11
xdata<-data[,c(2,3,4,6,7,8,9,12,13,15)]
xdata<-xdata[-11,]
model <- coxph( Surv(OS,status) ~ .,data =xdata)
fit<-summary(model)
ggforest(model,data=xdata)
library("meta")
library("ggplot2")

input<-data.frame(fit$conf.int)
write.table(input[,c(1,3,4)],file="ShareRatio_CoxRegression_95CI.txt",col.names = NA,row.names = T,sep="\t",quote=F)
input[3,4]<-10^25
input$Group<-rownames(input)
input$Condition<-1:nrow(input)
input$RiskRatio<-input[,1]
input$LowerLimit<-input[,3]
input$UpperLimit<-input[,4]

p = ggplot(data=input,aes(x = Group,y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+ 
  facet_wrap(~Condition,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p + scale_y_log10()+ guides(col = guide_legend(reverse = TRUE))
