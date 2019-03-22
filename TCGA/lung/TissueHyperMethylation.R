##### 
rltquantile<-apply(input,1,function(x) quantile(x,na.rm = T))
newrltquantile<-t(rltquantile[,which(rltquantile[2,]>0.6)])

map<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/hg19/GPL13534_450K_hg19.bed",sep="\t")
newrlt<-data.frame(map[match(rownames(newrltquantile),map[,4]),],tcga_quantile=newrltquantile)
head(newrlt)
write.table(newrlt,file="/home/guosa/hpc/methylation/methbase/tcga.hyper.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

# system('wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')
# system("gunzip GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz")
data<-read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct",head=T,sep="\t",skip=2)
rnaseq<-data.matrix(data[,3:ncol(data)])
ratio<-rnaseq[,ncol(rnaseq)]/(rowMeans(rnaseq[,1:(ncol(rnaseq)-1)]))
newdata<-data.frame(data,ratio)
BOG<-subset(newdata,ratio>20 & Whole.Blood>1)
BOGHyerMeth<-merge(newrlt,BOG,by.x="V5",by.y="Description")
newBOGHyperMeth<-BOGHyerMeth[,c(2:4,1,5:ncol(BOGHyerMeth))]

write.table(newBOGHyperMeth,file="newBOGHyperMeth.hg19.bed",sep="\t",quote=F,row.names = F,col.names = F)
