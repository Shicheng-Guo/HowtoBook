data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/US-CN/pubmed.txt")
data<-data[order(data[,2],decreasing = T),]
data
input<-as.vector(data[,3])
names(input)<-as.character(data[,1])
input
barplot(input,col=3:nrow(data),ylim=c(0,2000))

