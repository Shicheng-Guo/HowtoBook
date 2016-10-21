
data<-read.table("output.mf",head=T,sep="\t",row.names=1)
sam<-read.table("SamConfig.txt",sep="\t",head=T)

idx<-as.character(sam[match(colnames(data),sam[,4]),8])
idx[1:10]<-"iPS"

library("ggplot2")
input<-data.frame(me5=as.numeric(R2),idx)
input<-input[-grep("Day",input$idx),]
input$idx<-as.character(input$idx)

p <- ggplot(input, aes(factor(idx), me5))
p + geom_boxplot(aes(fill = factor(idx)))
p
dev.off()

R2=data[match("chr10:79857968-79858049",rownames(data)),]
boxplot(R2~idx,input)
