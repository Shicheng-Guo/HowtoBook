# Allele based analysis
# Allele frequencies of the CTG repeat in controls and schizophrenic patients, including subtypes, course specifiers and positive family history, and results of Fisher's exact test and Monte Carlo method


calcOddsRatio <- function(mymatrix,alpha=0.05,referencerow=3){
  rlt<-list()
  numrow <- nrow(mymatrix)
  myrownames <- rownames(mymatrix)
  
  for (ii in 1:numrow){
    rowname <- myrownames[ii]
    DiseaseUnexposed <- mymatrix[referencerow,1]
    ControlUnexposed <- mymatrix[referencerow,2]
    if (ii != referencerow){ 
      DiseaseExposed <- mymatrix[ii,1]
      ControlExposed <- mymatrix[ii,2] 
      totExposed <- DiseaseExposed + ControlExposed
      totUnexposed <- DiseaseUnexposed + ControlUnexposed
      probDiseaseGivenExposed <- DiseaseExposed/totExposed
      probDiseaseGivenUnexposed <- DiseaseUnexposed/totUnexposed
      probControlGivenExposed <- ControlExposed/totExposed
      probControlGivenUnexposed <- ControlUnexposed/totUnexposed
      
      pvalue=chisq.test(matrix(c(DiseaseExposed,ControlExposed,DiseaseUnexposed,ControlUnexposed),2,2,byrow=T))$p.value
      # calculate the odds ratio
      oddsRatio <- (probDiseaseGivenExposed*probControlGivenUnexposed)/(probControlGivenExposed*probDiseaseGivenUnexposed)
      # calculate a confidence interval
      confidenceLevel <- (1 - alpha)*100
      sigma <- sqrt((1/DiseaseExposed)+(1/ControlExposed)+(1/DiseaseUnexposed)+(1/ControlUnexposed))
      # sigma is the standard error of our estimate of the log of the odds ratio
      z <- qnorm(1-(alpha/2))
      lowervalue <- oddsRatio * exp(-z * sigma)
      uppervalue <- oddsRatio * exp( z * sigma)
      OR=paste(rownames(mymatrix)[ii],"      ",round(oddsRatio,2)," (",round(lowervalue,2),",",round(uppervalue,2),")","          P=",pvalue,sep="")
      print(OR)
    }
  }
  
}

input<-data.frame(control=c(38,72,78,12,0),case=c(42,67,80,14,1))
input
rownames(input)<-c(6,9,10,11,13)
input
calcOddsRatio(input,referencerow=3)
