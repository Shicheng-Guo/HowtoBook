# http://GettingGeneticsDone.blogspot.com/
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html
# 

# Define the function
install.packages("Haplin")
library("Haplin")
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

#Generate some fake data that deviates from the null

set.seed(42)
par(mfrow=c(2,2))
pvalues=runif(10000)

# NULL UNIF distribution
pQQ(pvalues, nlabs =200, conf = 0.95, mark = F) 

# UNIF distribution with outliner
pvalues[sample(10000,10)]=pvalues[sample(10000,10)]/5000
pQQ(pvalues, nlabs =2000, conf = 0.95) 
