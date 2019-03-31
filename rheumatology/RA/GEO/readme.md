
## Differential Gene expression (DEG) for FLS from RA 

Postive control:  CSF2RB, PSM9 over expressed in RA while JUND, NFIL3 low expressed in RA

339 DEGs were derived from 11 RA FLS vs 11 OA by RNA-seq published in Nature Communication in 2018

Suppose LLDT-8 affect the therapy of RA, then LLDT-8 should recover parts of these 339 genes. 

Xinpeng mentioned a paper published in ART in 2016 about lncRNA expression in FLS from RA patients. 

Interesting thing is they mentioned 10 samples were collected in the study. However, when I check the data in GEO: GSE83147, only 3 samples were used for the lncRNA-array. 

https://arthritis-research.biomedcentral.com/articles/10.1186/s13075-016-1129-4#Sec2

Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83147

wget -e robots=off -nH -r -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83147/suppl/


```
system("wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/GWAS-reported-gene.txt")
system("wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/RA_FLS_DEG_11vs11_Nature_Commnication.txt")
data1<-read.table("GWAS-reported-gene.txt",head=T,sep="\t")
data2<-read.table("RA_FLS_DEG_11vs11_Nature_Commnication.txt",head=T,sep="\t")
GwasGene<-unlist(strsplit(as.character(data1[,2]),"[,| ]"))
DMERGene<-as.character(data2[,1])
unique(GwasGene[which(GwasGene %in% DMERGene)])
```
combine GWAS-SNP-Gene and DEG in NC, we found 4 genes are overlapped: "CDH6"  "CFTR"  "ADCY8" "ACOXL"


```
a=/gpfs/home/guosa/hpc/db/hg19/refGene.hg19.bed
b=/gpfs/home/guosa/nc/GWAS-RA-378.GRCH37.bed
bedtools window -a $b -b $a -w 500000 | awk '{print $9}' | sort -u  >  GWAS-378SNP-50K-Gene.txt
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/RA/RA_FLS_DEG_11vs11_Nature_Commnication.txt
awk '{print $1}' RA_FLS_DEG_11vs11_Nature_Commnication.txt | sort -u > DMER.GENE.List.txt
comm GWAS-378SNP-50K-Gene.txt  DMER.GENE.List.txt -1 -2
```
take 1M as the extension distance, we can find 35 overlap genes between GWAS-SNP and DEGs. This LINC00565 is quite interesting, it is significantly over-expressed in RA FLS and is a significant lincRNA identified by GWAS study.  We can check whether It can be downregulated by LLDT-8 and then we can report this story. 

ABCA8
ACOXL
ACYP1
AKR1C3
ARHGEF3-AS1
BCL2L11
BIRC3
CDH6
CFTR
CREB5
DSCR8
ESM1
FAM216B
FRAS1
GALNT15
GCSAM
GPR157
GRTP1-AS1
IL12A
IL12A-AS1
KCNJ15
LIMK2
LINC00565
LPIN2
MCF2L
MINK1
NUDT6
PGR
RTTN
SH2D5
SLC2A5
SLC9C1
SMOC2
TYR
VIT



totally GEO dataset which can be combined together:

mRNA-array:
1. GSE36700
2. GSE83147
3. GSE55457
4. GSE55584
5. GSE55235	
6. GSE7669
7. GSE1919
8. GSE2053

methylation array:
1. GSE46650
2. GPL13534

Rheumatoid arthritic synovium response to methotrexate

Rheumatoid arthritic synovium response to tocilizumab

Rheumatoid arthritic synovium response to rituximab

RNA-seq dataset: GSE29127



Some other dataset:
```
FULLNAME	DATASET
Rheumatoid Arthritis	GSE12051
Rheumatoid Arthritis	GSE1919
Rheumatoid Arthritis	GSE1402
Rheumatoid Arthritis	GSE48006
Rheumatoid Arthritis	GSE55468
Rheumatoid Arthritis	GSE48780
Rheumatoid Arthritis	GSE55457
Rheumatoid Arthritis	GSE55584
Rheumatoid Arthritis	GSE55235
Rheumatoid Arthritis	GSE50646
Rheumatoid Arthritis	GSE49604
Rheumatoid Arthritis	GSE42296
Rheumatoid Arthritis	GSE45536
Rheumatoid Arthritis	GSE39340
Rheumatoid Arthritis	GSE39428
Rheumatoid Arthritis	GSE36700
Rheumatoid Arthritis	GSE29746
Rheumatoid Arthritis	GSE27390
Rheumatoid Arthritis	GSE22956
Rheumatoid Arthritis	GSE17755
Rheumatoid Arthritis	GSE22108
Rheumatoid Arthritis	GSE15258
Rheumatoid Arthritis	GSE12653
Rheumatoid Arthritis	GSE15573
Rheumatoid Arthritis	GSE13026
Rheumatoid Arthritis	GSE12021
Rheumatoid Arthritis	GSE9027
Rheumatoid Arthritis	GSE3698
Rheumatoid Arthritis	GSE3824
Rheumatoid Arthritis	GSE1919
Rheumatoid Arthritis	GSE1402
```
When we combined different GEO dataset and we found: 
```

input<-read.table("candidate.txt",sep="\t",head=T)
data<-read.table("GSE55457.txt",sep="\t",head=T)
out<-data.frame(input,data[match(input$rightID,data$Gene.symbol),])
x1<-subset(out,regulation.m=="down"  & as.numeric(as.character(logFC))>0)
x2<-subset(out,regulation.m=="up"  & as.numeric(as.character(logFC))<0)
xx<-rbind(x1,x2)
write.table(xx,file="GSE55457.rlt.txt",sep="\t",quote=F,col.names = T,row.names = T)


input<-read.table("candidate.txt",sep="\t",head=T)
data<-read.table("GSE55235.txt",sep="\t",head=T)
out<-data.frame(input,data[match(input$rightID,data$Gene.symbol),])
x1<-subset(out,regulation.m=="down"  & as.numeric(as.character(logFC))>0)
x2<-subset(out,regulation.m=="up"  & as.numeric(as.character(logFC))<0)
xx<-rbind(x1,x2)
write.table(xx,file="GSE55235.rlt.txt",sep="\t",quote=F,col.names = T,row.names = T)

input<-read.table("candidate.txt",sep="\t",head=T)
data<-read.table("GSE7669.txt",sep="\t",head=T)
out<-data.frame(input,data[match(input$rightID,data$Gene.symbol),])
x1<-subset(out,regulation.m=="down"  & as.numeric(as.character(logFC))>0)
x2<-subset(out,regulation.m=="up"  & as.numeric(as.character(logFC))<0)
xx<-rbind(x1,x2)
write.table(xx,file="GSE7669.rlt.txt",sep="\t",quote=F,col.names = T,row.names = T)


input<-read.table("candidate.txt",sep="\t",head=T)
data<-read.table("GSE55584.txt",sep="\t",head=T)
out<-data.frame(input,data[match(input$rightID,data$Gene.symbol),])
x1<-subset(out,regulation.m=="down"  & as.numeric(as.character(logFC))>0)
x2<-subset(out,regulation.m=="up"  & as.numeric(as.character(logFC))<0)
xx<-rbind(x1,x2)
write.table(xx,file="GSE55584.rlt.txt",sep="\t",quote=F,col.names = T,row.names = T)

input<-read.table("candidate.txt",sep="\t",head=T)
data<-read.table("GSE49604.txt",sep="\t",head=T)
out<-data.frame(input,data[match(input$rightID,data$Gene.symbol),])
x1<-subset(out,regulation.m=="down"  & as.numeric(as.character(logFC))>0)
x2<-subset(out,regulation.m=="up"  & as.numeric(as.character(logFC))<0)
xx<-rbind(x1,x2)
write.table(xx,file="GSE49604.rlt.txt",sep="\t",quote=F,col.names = T,row.names = T)
```




