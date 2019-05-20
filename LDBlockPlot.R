# LD block with R
# LD block with R
install.packages("gaston")
library("gaston")
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
# Compute LD
ld.x <- LD(x, c(1,ncol(x)))
# Plot a tiny part of the LD matrix
LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20] )
# Customize the plot
LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20], 
         graphical.par = list(cex =1.5, bg = "gray"), 
         polygon.par = list(border = NA), write.ld = NULL )

