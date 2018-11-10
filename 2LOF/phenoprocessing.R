setwd("/gpfs/home/guosa/hpc/project/pmrp/phen")

file=list.files(pattern="*_SampleIDs.txt")
file

# Obesity
phen1<-read.table("IBDCH_Michigen_ID.txt",head=T)
phen2<-read.table("IBDCH_Phetyp10_Obesity_SampleIDs.txt",head=T)
which(is.na(match(phen1[,1],phen2[,1])))
phen=cbind(phen1,phen2[match(phen1[,1],phen2[,1]),])

# RA
file
phen1<-read.table("IBDCH_Michigen_ID.txt",head=T)
phen2<-read.table(file[2],head=T,sep="\t")
head(phen1)
head(phen2)
match(phen1[,1],phen2[,1])
phen=cbind(phen1,phen2[match(phen1[,1],phen2[,1]),])
write.table(phen,"IBDCH_Phetyp1_RA_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)

# PA
file
phen1<-read.table("IBDCH_Michigen_ID.txt",head=T)
phen2<-read.table(file[3],head=T,sep="\t")
head(phen1)
head(phen2)
match(phen1[,1],phen2[,1])
phen=cbind(phen1,phen2[match(phen1[,1],phen2[,1]),])
write.table(phen,"IBDCH_Phetyp3_PA_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)

# Tyroid
phen1<-read.table("IBDCH_Michigen_ID.txt",head=T)
phen2<-read.table(file[4],head=T,sep="\t")
head(phen1)
head(phen2)
match(phen1[,1],phen2[,1])
phen=cbind(phen1,phen2[match(phen1[,1],phen2[,1]),])
write.table(phen[,c(1,2,3,4,5,7,8)],"IBDCH_Phetyp5_Thyroid_C1_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)
write.table(phen[,c(1,2,3,4,6,7,8)],"IBDCH_Phetyp5_Thyroid_C2_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)


# SSc
file
phen1<-read.table("IBDCH_Michigen_ID.txt",head=T)
phen2<-read.table(file[5],head=T,sep="\t")
head(phen1)
head(phen2)
match(phen1[,1],phen2[,1])
phen=cbind(phen1,phen2[match(phen1[,1],phen2[,1]),])
head(phen)
write.table(phen[,c(1,2,3,4,5,8,9)],"IBDCH_Phetyp6_SSc_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)
write.table(phen[,c(1,2,3,4,6,8,9)],"IBDCH_Phetyp6_ANA_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)
write.table(phen[,c(1,2,3,4,7,8,9)],"IBDCH_Phetyp6_ENA_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)

# IRON
file
phen1<-read.table("IBDCH_Michigen_ID.txt",head=T)
phen2<-read.table(file[6],head=T,sep="\t")
head(phen1)
head(phen2)
match(phen1[,1],phen2[,1])
phen=cbind(phen1,phen2[match(phen1[,1],phen2[,1]),])
head(phen)
write.table(phen[,c(1,2,3,4,5,7,8)],"IBDCH_Phetyp7_Iron_C1_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)
write.table(phen[,c(1,2,3,4,6,7,8)],"IBDCH_Phetyp7_Iron_C2_rev2_SampleIDs.Michigen.txt",sep="\t",col.names = T,row.names=F,quote=F)





