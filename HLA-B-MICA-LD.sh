# hg19 genomic region for HLA-B and MICA
# chr6:31,232,075-31,391,038

cd /home/local/MFLDCLIN/guosa/hpc/db/1000Genome

awk '$7=="CEU"' integrated_call_samples_v2.20130502.ALL.ped > CEU.txt
plink --bfile chr6 --make-bed --keep CEU.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CEU.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()


awk '$7=="CHB"' integrated_call_samples_v2.20130502.ALL.ped > CHINA.txt
awk '$7=="CHS"' integrated_call_samples_v2.20130502.ALL.ped >> CHINA.txt
plink --bfile chr6 --make-bed --keep CHINA.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CHINA.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()

