1. merge CpG-SNP list
```
for i in `ls *cpgsnp.bed`
do
awk '{print $4}' $i >> cpgSNP.hg19.list.txt
done
```
2. extract CpG-SNPs from 1000 Genome Project
```
cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr$i --extract /home/guosa/hpc/db/hg19/plan2/commonCpGSNP/cpgSNP.hg19.list.txt --make-bed --out chr$i.cpgSNP
done
```
```
perl swith2csv.pl fastphase_hapguess_switch.out > AHRR.diplo.txt
```
```
data<-read.table("AHRR.hap.txt",row.names=1)
write.table(table(data),file="AHRR.dip.txt",sep="\t",quote=F,row.names=T,col.names=NA)
```

```
setwd("/home/guosa/hpc/methylation/methtype")
data<-read.table("AHRR.hap.txt",row.names=1)
write.table(table(data),file="AHRR.dip.txt",sep="\t",quote=F,row.names=T,col.names=NA)

rw=unlist(lapply(strsplit(rownames(table(data)),split=""),function(x) sum(x %in% "C")))
cw=unlist(lapply(strsplit(colnames(table(data)),split=""),function(x) sum(x %in% "C")))
rrw=tapply(rowSums(table(data)),rw,sum)
ccw=tapply(colSums(table(data)),cw,sum)
vv=rbind(data.frame(v=rrw,w=names(rrw)),data.frame(v=ccw,w=names(ccw)))
tapply(vv$v,vv$w,sum)
```
