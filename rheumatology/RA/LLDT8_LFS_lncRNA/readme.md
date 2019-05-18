```
library("ggplot2")
library("extrafont")
input<-data.frame(item=c("DMSO","TNF","LLDT8"),mean=c(14.7,20.51,31.93),sd=c(4.02,4.25,4.25))
input$item <- factor(input$item,levels = c("DMSO", "TNF", "LLDT8"))
p<- ggplot(input, aes(x=item, y=mean, fill=item)) 
p<- p+ geom_bar(stat="identity", color="black",position=position_dodge(),size=1.1) 
p<- p+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),size=1.1) 
p<- p+labs(title="", x="", y = "Apoptosis Ratio (TUNEL)")
p<- p+ theme_classic()
p<-p+theme(text = element_text(size=25,family="TT Arial"))
p<- p+scale_fill_manual(values=c('red','green','blue'))
print(p)
```
