While the initial number of variants was stated as being ‘over 500,000’, the subsequent steps of quality control (i.e. genotype call rates, deviation from Hardy-Weinberg Equilibrium (HWE)) were just mentioned in the passing, and how these affected the final distribution of the variants based on frequency (i.e. how many is <= 1% (rare variants), 1-10% (moderately common variants) and >= 10% (common variants)) wasn’t mentioned. Can the authors include this information within text? This would be useful to get a sense on how many rare variants there were left before going into phasing as well as a the final number of variants used in the analysis.
```
plink --bfile ../FinalRelease_QC_20140311_Team1_Marshfield --mind 0.05 --geno 0.95 --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO --freq --out plink1
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO --extract fSNP.bed --range --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO.fSNP
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield.MIND.GENO.fSNP --freq --out plink2

data<-read.table("plink1.frq",sep="",head=T)
L1<-nrow(subset(data,MAF==0))
L2<-nrow(subset(data,MAF>0 & MAF<0.001))
L3<-nrow(subset(data,MAF>0.001 & MAF<0.005))
L4<-nrow(subset(data,MAF>0.005 & MAF<0.01))
L5<-nrow(subset(data,MAF>0.01 & MAF<0.05))
L6<-nrow(subset(data,MAF>0.05))
x1<-c(L1,L2,L3,L4,L5,L6)

data<-read.table("plink2.frq",sep="",head=T)
L1<-nrow(subset(data,MAF==0))
L2<-nrow(subset(data,MAF>0 & MAF<0.001))
L3<-nrow(subset(data,MAF>0.001 & MAF<0.005))
L4<-nrow(subset(data,MAF>0.005 & MAF<0.01))
L5<-nrow(subset(data,MAF>0.01 & MAF<0.05))
L6<-nrow(subset(data,MAF>0.05))
x2<-c(L1,L2,L3,L4,L5,L6)
```
