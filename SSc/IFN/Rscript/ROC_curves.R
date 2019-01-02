#ROC curve - Fig. 3B
load("CD4_RGSS_data_aftercombatM.RData")
cg.comparison = read.csv("cg.comparison_aftercombatM.csv", header = T,stringsAsFactors = F)
DMS = cg.comparison[which(cg.comparison[,3] < 0.01),]
seq1 = which(DMS$UCSC_RefGene_Name == "MX1")
MX1 = DMS[seq1,]
MX1_beta = as.data.frame(CD4_RGSS_data_aftercombatM[MX1$X,])
seq2 = which(DMS$UCSC_RefGene_Name == "OAS1")
OAS1 = DMS[seq2,]
OAS1_beta = as.data.frame(CD4_RGSS_data_aftercombatM[OAS1$X,])
seq3 = which(DMS$UCSC_RefGene_Name == "USP18")
USP18 = DMS[seq3,]
USP18_beta = as.data.frame(CD4_RGSS_data_aftercombatM[USP18$X,])
seq4 = which(DMS$UCSC_RefGene_Name == "IFIT1")
IFIT1 = DMS[seq4,]
IFIT1_beta = as.data.frame(CD4_RGSS_data_aftercombatM[IFIT1$X,])
seq5 = which(DMS$UCSC_RefGene_Name == "IRF7")
IRF7 = DMS[seq5,]
IRF7_beta = as.data.frame(CD4_RGSS_data_aftercombatM[IRF7$X,])
seq6 = which(DMS$UCSC_RefGene_Name == "RSAD2")
RSAD2 = DMS[seq6,]
RSAD2_beta = as.data.frame(CD4_RGSS_data_aftercombatM[RSAD2$X,])
IFN_data = rbind(MX1_beta,OAS1_beta,USP18_beta,IFIT1_beta,IRF7_beta,RSAD2_beta)
IFN_data = as.data.frame(t(IFN_data))
clinical = read.csv("CD4_RGSS_clinical.csv", header = T)
IFN_ROC = cbind(IFN_data,clinical[,2:3])
write.csv(IFN_ROC,file="IFN_CD4_aftercombatM_ROC.csv",quote=F)

library(readr)
library(pROC)
library(ggsci)
library(ggplot2)
library(ggthemes)
dat_cd4 = read_csv("IFN_CD4_aftercombatM_ROC.csv")
dat_cd4 = dat_cd4[,-1]
dat_cd4_1 = dat_cd4
colnames(dat_cd4_1)[22] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = dat_cd4_1[,-23]
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_IFNsites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_IFNsites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curves - Supplementary Fig. 5
dat_cd4_2 = dat_cd4
colnames(dat_cd4_2)[22] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)
Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-23]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_IFNsites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_IFNsites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 6A
dat_cd4_1 = dat_cd4[,c(14:15,22)]
colnames(dat_cd4_1)[3] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_IFIT1sites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_IFIT1sites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 6B
dat_cd4_1 = dat_cd4[,c(16:18,22)]
colnames(dat_cd4_1)[4] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_IRF7sites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_IRF7sites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 6C
dat_cd4_1 = dat_cd4[,c(1:6,22)]
colnames(dat_cd4_1)[7] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_MX1sites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_MX1sites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 6D
dat_cd4_1 = dat_cd4[,c(7:9,22)]
colnames(dat_cd4_1)[4] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_OAS1sites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_OAS1sites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 6E
dat_cd4_1 = dat_cd4[,c(10:13,22)]
colnames(dat_cd4_1)[5] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_USP18sites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_USP18sites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 6F
dat_cd4_1 = dat_cd4[,c(19:22)]
colnames(dat_cd4_1)[4] = "Type"
dat_cd4_1$Type = ifelse(dat_cd4_1$Type == "Case", 1, 0)
dat_cd4_1 = na.exclude(dat_cd4_1)
glm.data = glm(Type~., data=dat_cd4_1, family="binomial")
Methy.rocobj  = roc(dat_cd4_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_RSAD2sites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_RSAD2sites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 7
dat_cd4_2 = dat_cd4[,c(14:15,22:23)]
colnames(dat_cd4_2)[3] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-4]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_IFIT1sites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_IFIT1sites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 8
dat_cd4_2 = dat_cd4[,c(16:18,22:23)]
colnames(dat_cd4_2)[4] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-5]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_IRF7sites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_IRF7sites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 9
dat_cd4_2 = dat_cd4[,c(1:6,22:23)]
colnames(dat_cd4_2)[7] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-8]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_MX1sites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_MX1sites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 10
dat_cd4_2 = dat_cd4[,c(7:9,22:23)]
colnames(dat_cd4_2)[4] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-5]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_OAS1sites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_OAS1sites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 11
dat_cd4_2 = dat_cd4[,c(10:13,22:23)]
colnames(dat_cd4_2)[5] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-6]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_USP18sites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_USP18sites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 12
dat_cd4_2 = dat_cd4[,c(19:23)]
colnames(dat_cd4_2)[4] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
dat_cd4_2 = na.exclude(dat_cd4_2)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,-5]
  glm.data = glm(Type~., data=temp, family="binomial", control = list(maxit = 50))
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_RSAD2sites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_RSAD2sites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Fig. 3C-3F
load("CD4_RGSS_data_aftercombatM.RData")
cg.comparison = read.csv("cg.comparison_aftercombatM.csv", header = T,stringsAsFactors = F)
DMS = cg.comparison[which(cg.comparison[,3] < 0.01),]
seq = which(DMS$UCSC_RefGene_Name == "IFI44L")
IFI44L = DMS[seq,]
IFI44L_beta = as.data.frame(CD4_RGSS_data_aftercombatM[IFI44L$X,])
IFI44L_data = as.data.frame(t(IFI44L_beta))
clinical = read.csv("CD4_RGSS_clinical.csv", header = T)
IFI44L_ROC = cbind(IFI44L_data,clinical[,2:3])
write.csv(IFI44L_ROC,file="IFI44L_CD4_aftercombatM_ROC.csv",quote=F)

library(readr)
library(pROC)
library(ggsci)
library(ggplot2)
library(ggthemes)
dat_cd4 = read_csv("IFI44L_CD4_aftercombatM_ROC.csv")
dat_cd4 = dat_cd4[,-1]
dat_cd4_3 = dat_cd4
colnames(dat_cd4_3)[7] = "Type"
dat_cd4_3$Type = ifelse(dat_cd4_3$Type == "Case", 1, 0)

Diseases = unique(dat_cd4_3$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_3, dat_cd4_3$Disease_Type == Diseases[i])
  temp = temp[,-8]
  glm.data = glm(Type~., data=temp, family="binomial")
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_6CpGsites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_6CpGsites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 15A
dat_cd4_3 = dat_cd4
colnames(dat_cd4_3)[7] = "Type"
dat_cd4_3$Type = ifelse(dat_cd4_3$Type == "Case", 1, 0)
dat_cd4_3 = dat_cd4_3[,-8]
glm.data = glm(Type~., data=dat_cd4_3, family="binomial")
Methy.rocobj  = roc(dat_cd4_3$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_6CpGsites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_6CpGsites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 13A
dat_cd4_2 = dat_cd4[,c(4,7)]
colnames(dat_cd4_2)[2] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)
glm.data = glm(Type~., data=dat_cd4_2, family="binomial")
Methy.rocobj  = roc(dat_cd4_2$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD4_1CpGsite_Total", ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD4_1CpGsite_Total", ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 14A-14D
dat_cd4_2 = dat_cd4[,c(4,7:8)]
colnames(dat_cd4_2)[2] = "Type"
dat_cd4_2$Type = ifelse(dat_cd4_2$Type == "Case", 1, 0)

Diseases = unique(dat_cd4_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd4_2, dat_cd4_2$Disease_Type == Diseases[i])
  temp = temp[,c(1,2)]
  glm.data = glm(Type~., data=temp, family="binomial")
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD4_1CpGsite_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD4_1CpGsite_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Fig. 3G-3H
load("CD8_GS_data_aftercombatM.RData")
cg.comparison = read.csv("cg.comparison_aftercombatM.csv", header = T,stringsAsFactors = F)
DMS = cg.comparison[which(cg.comparison[,3] < 0.01),]
seq = which(DMS$UCSC_RefGene_Name == "IFI44L")
IFI44L = DMS[seq,]
IFI44L_beta = as.data.frame(CD8_GS_data_aftercombatM[IFI44L$X,])
IFI44L_data = as.data.frame(t(IFI44L_beta))
site = as.data.frame(CD8_GS_data_aftercombatM[match("cg06872964",row.names(CD8_GS_data_aftercombatM)),])
colnames(site) = "cg06872964"
IFI44L_all = cbind(IFI44L_data, site)
clinical = read.csv("CD8_GS_clinical.csv", header = T)
IFI44L_ROC = cbind(IFI44L_all,clinical[,2:3])
write.csv(IFI44L_ROC,file="IFI44L_CD8_aftercombatM_ROC.csv",quote=F)

dat_cd8 = read_csv("IFI44L_CD8_aftercombatM_ROC.csv")
dat_cd8 = dat_cd8[,-1]
dat_cd8_1 = dat_cd8[,c(1:2,4:5)]
colnames(dat_cd8_1)[3] = "Type"
dat_cd8_1$Type = ifelse(dat_cd8_1$Type == "Case", 1, 0)

Diseases = unique(dat_cd8_1$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd8_1, dat_cd8_1$Disease_Type == Diseases[i])
  temp = temp[,1:3]
  glm.data = glm(Type~., data=temp, family="binomial")
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p=ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD8_2CpGsites_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD8_2CpGsites_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}

#ROC curve - Supplementary Fig. 15B
dat_cd8_1 = dat_cd8[,c(1:2,4)]
colnames(dat_cd8_1)[3] = "Type"
dat_cd8_1$Type = ifelse(dat_cd8_1$Type == "Case", 1, 0)
glm.data = glm(Type~., data=dat_cd8_1, family="binomial")
Methy.rocobj  = roc(dat_cd8_1$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD8_2CpGsites_Total",  ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD8_2CpGsites_Total",  ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 13B
dat_cd8_2 = dat_cd8[,c(3,4)]
colnames(dat_cd8_2)[2] = "Type"
dat_cd8_2$Type = ifelse(dat_cd8_2$Type == "Case", 1, 0)
glm.data = glm(Type~., data=dat_cd8_2, family="binomial")
Methy.rocobj  = roc(dat_cd8_2$Type, glm.data$fitted.values,smooth = F)
Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
idx.bestcutoff = which.max(Methy.rocdata$combine)
best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")

p = ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
  geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
  xlab("1- Specificity")+ylab("Sensitivity")+ 
  geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
  scale_color_npg()+
  annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
  annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
  annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
  scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(plot.title=element_text(colour="black",size=18,face="bold"))+
  theme(axis.line = element_line(color = 'black',size=1.05))+ 
  theme(axis.text.x=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
ggsave(filename = paste("CD8_1CpGsite_Total", ".pdf",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="pdf")
ggsave(filename = paste("CD8_1CpGsite_Total", ".tiff",sep=""), plot = p, 
       width=7, height = 6.5, units="in", device ="tiff")

#ROC curve - Supplementary Fig. 14E-14F
dat_cd8_2 = dat_cd8[,3:5]
colnames(dat_cd8_2)[2] = "Type"
dat_cd8_2$Type = ifelse(dat_cd8_2$Type == "Case", 1, 0)

Diseases = unique(dat_cd8_2$Disease_Type)
for(i in 1:length(Diseases)){
  temp = subset(dat_cd8_2, dat_cd8_2$Disease_Type == Diseases[i])
  temp = temp[,c(1,2)]
  glm.data = glm(Type~., data=temp, family="binomial")
  Methy.rocobj  = roc(temp$Type, glm.data$fitted.values,smooth = F)
  Methy.rocdata = data.frame(Sens = Methy.rocobj$sensitivities, Spec = Methy.rocobj$specificities)
  Methy.rocdata$combine = Methy.rocdata[,1] + Methy.rocdata[,2]
  idx.bestcutoff = which.max(Methy.rocdata$combine)
  best.sens = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,1]), digits = 2, format = "f")
  best.spec = formatC(as.numeric(Methy.rocdata[idx.bestcutoff,2]), digits = 2, format = "f")
  AUC = formatC(as.numeric(Methy.rocobj$auc[[1]]), digits = 2, format = "f")
  
  p=ggplot(Methy.rocdata, aes(x = 1-Sens, y = Spec)) + 
    geom_line(size=2.5,colour ="#DC0000B2")+theme_pander()+ 
    xlab("1- Specificity")+ylab("Sensitivity")+ 
    geom_abline(intercept=0,slope=1 ,colour="black",linetype=4,size=1.2)+
    scale_color_npg()+
    annotate("text", x = 0.70, y = 0.30, label = paste("AUC  = ",AUC,""),size = 10)+
    annotate("text", x = 0.70, y = 0.22, label = paste("Sens  = ",best.spec, ""),size = 10)+
    annotate("text", x = 0.70, y = 0.14, label = paste("Spec  = ",best.sens, ""),size = 10)+
    scale_x_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    scale_y_continuous(limits=c(0,1.0), breaks=seq(0,1,0.25))+
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
    theme(plot.title=element_text(colour="black",size=18,face="bold"))+
    theme(axis.line = element_line(color = 'black',size=1.05))+ 
    theme(axis.text.x=element_text(size=15,face="bold"))+
    theme(axis.title.x=element_text(size=20,vjust=2,face="bold"))+
    theme(axis.text.y=element_text(size=15,face="bold"))+
    theme(axis.title.y=element_text(size=20,vjust=0,face ="bold"))
  ggsave(filename = paste("CD8_1CpGsite_", Diseases[i], ".pdf",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="pdf")
  ggsave(filename = paste("CD8_1CpGsite_", Diseases[i], ".tiff",sep=""), plot = p, 
         width=7, height = 6.5, units="in", device ="tiff")
  
}
