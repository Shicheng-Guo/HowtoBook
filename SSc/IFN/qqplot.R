# http://review.frontiersin.org/review/433093/18/232809/#tab/History
# Shicheng Guo
install.packages("Haplin")
library("Haplin")
cd4<-read.table("CD4.txt",head=T)
pQQ(cd4[,2], nlabs =nrow(cd4), conf = 0.95) 
cd8<-read.table("CD8.txt",head=T)
pQQ(cd8[,2], nlabs =nrow(cd8), conf = 0.95) 
