for(i in 1:100){
x<-paste("cp ",out1[i,1],"-* ./pick",sep="")
system(x)
}
