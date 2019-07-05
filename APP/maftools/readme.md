```
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")


install.packages("BiocManager")
BiocManager::install("maftools")

setwd("D:\\LungBrain")
laml<-annovarToMaf(annovar="T1.hg19_multianno.txt", Center = NULL, refBuild = "hg19",
             tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
             MAFobj = F, sampleAnno = NULL)
laml = read.maf(maf = laml)
pdf("Lung.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
getGeneSummary(laml)
oncoplot(maf = laml, top = 10)
oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))
```
