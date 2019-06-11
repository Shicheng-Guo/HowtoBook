sheet1= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 1)
sheet2= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 2)
sheet1= as.data.frame(sheet1)
sheet2= as.data.frame(sheet2)
head(sheet1)
head(sheet2)

index<-unlist(lapply(sheet2$Filename,function(x) unlist(strsplit(x,"_"))[1]))
Mean<-tapply(sheet2$`Strelka+Vardict`,index,function(x) mean(x))
write.table(Mean,file="Mean.MutationLoad.txt",sep="\t",quote=F,col.names = NA,row.names = T)


sheet1= read_excel("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/prostate/Result.xlsx",sheet = 1)
sheet1= as.data.frame(sheet1)
head(sheet1)
summary(lm(Age~MLxMillion,sheet1))$coefficients[2,4]
summary(lm(tPSA~MLxMillion,sheet1))$coefficients[2,4]
summary(lm(fPSA~MLxMillion,sheet1))$coefficients[2,4]
summary(lm(fPSA_tPSA~MLxMillion,sheet1))$coefficients[2,4]
summary(lm(Gleason_Score_Value~MLxMillion,sheet1))$coefficients[2,4]
summary(glm(Metastasis_binary~MLxMillion,sheet1,family=binomial))$coefficients[2,4]

