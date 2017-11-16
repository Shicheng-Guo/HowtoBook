
#install.packages("quantmod")
#install.packages("TTR")
library(quantmod)
library(TTR)
myATR        <- function(x) ATR(HLC(x))[,'atr']
mySMI        <- function(x) SMI(HLC(x))[, "SMI"]
myADX        <- function(x) ADX(HLC(x))[,'ADX']
myAroon      <- function(x) aroon(cbind(Hi(x),Lo(x)))$oscillator
myBB         <- function(x) BBands(HLC(x))
myChaikinVol <- function(x) Delt(chaikinVolatility(cbind(Hi(x),Lo(x))))[, 1]
myEMA        <- function(x,n) EMA(HLC(x)[,3],n)
myEMV        <- function(x) EMV(cbind(Hi(x),Lo(x)),Vo(x))[,2]
myMACD       <- function(x) MACD(Cl(x),percent=F,nFast=12, nSlow=26, nSig=9)
myMFI        <- function(x) MFI(HLC(x), Vo(x))
mySAR        <- function(x) SAR(cbind(Hi(x),Cl(x))) [,1]
myVolat      <- function(x) volatility(OHLC(x),calc="garman")[,1]
myCCI        <- function(x) CCI(x,n = 20, maType, c = 0.015)
myRSI        <- function(x) RSI(x[,4], n=12)

beginning <- as.Date(Sys.Date()-540)
date<-Sys.Date()

setwd("/media/Home_Raid1/shg047/bak/bak2")
input<-read.table("bak.db",sep="\t")

#setwd("/home/sguo/Dropbox/stock/mystock")
#input<-read.table("/home/sguo/Dropbox/stock/bak.db",sep="\t")

#setwd("C:\\Users\\shg047\\Dropbox\\stock\\mystock")
#input<-read.table("C:\\Users\\shg047\\Dropbox\\stock\\bak.db",sep="\t")

for(i in 1:nrow(input)){
  symbol=as.character(input[i,1])
  try(Object<-getSymbols(Symbols=as.character(symbol),src="google",env=NULL,from=beginning,to=date,auto.assign=FALSE),silent = T)
  if(nrow(na.omit(Object))<350) next
  if(any(is.na(tail(Object,1)[1,1:3]))) next
  if(nrow(na.omit(Object))<150) next
  if((tail(Object,1)[1,5])<100000) next
  if((tail(Object,1)[1,1])>500) next
  if((tail(Object,1)[1,1])<10) next
  
  #myEMAD5<-myEMA(na.omit(Object),n=5)
  #myEMAD10<-myEMA(na.omit(Object),n=10)
  #myEMAD20<-myEMA(na.omit(Object),n=20)
  #myEMAD60<-myEMA(na.omit(Object),n=60)
  myEMAD120<-myEMA(na.omit(Object),n=120)
  #myEMAD250<-myEMA(na.omit(Object),n=250)
  #if(tail(Object,1)[1,4]<tail(myEMAD250,1)[1,1]) next
  #if(tail(Object,1)[1,4]<tail(myEMAD120,1)[1,1]) next
  #if(tail(Object,1)[1,4]<tail(myEMAD60,1)[1,1]) next
  
  #if(! all(tail(myEMAD5,10)[,1]>tail(myEMAD10,10)[,1])) next
  #if(! all(tail(myEMAD10,10)[,1]>tail(myEMAD30,10)[,1])) next
  #if(! all(tail(myEMAD30,10)[,1]>tail(myEMAD60,10)[,1])) next
  BBobject<-myBB(na.omit(Object))
  head(BBobject)
  if(tail(Object,1)[1,4]<tail(myEMAD120,1)) next
  if(mean(tail(BBobject,3)[,3]-tail(BBobject,3)[,1])/mean(tail(BBobject,30)[,3]-tail(BBobject,30)[,1])<1.1) next
  #if(! tail(myEMAD5,1)[1,1]>=tail(myEMAD10,1)[1,1] && tail(myEMAD5,5)[1,1]<=tail(myEMAD10,5)[1,1]) next
  if(var(tail(BBobject,30)[1:27,3]-tail(BBobject,30)[1:27,1])>0.3) next
  if(mean(tail(Object[,4],3),na.rm=T)/mean(tail(Object[,4],30),na.rm=T)<1.05) next
  #if(mean(tail(BBobject,2)[,3]-tail(BBobject,2)[,1])/mean(tail(BBobject,30)[,3]-tail(BBobject,30)[,1])>1.20) next
  #if(! tail(myEMAD10,1)[1,1]>=tail(myEMAD20,1)[1,1] && tail(myEMAD10,5)[1,1]<=tail(myEMAD20,5)[1,1]) next
  # daily to weekly or monthly
  # Object1<-to.weekly(Object)
  # Object1<-to.monthly(Object)
  #MACDobject<-myMACD(Cl(Object1))    
  #cond1<- as.numeric(MACDobject[nrow(MACDobject),1]) >= -1 && as.numeric(MACDobject[nrow(MACDobject)-2,1]) < -2
  # MSIDobject<-mySMI(na.omit(Object))
  # BBobject<-myBB(na.omit(Object1))         
  #CCIobject<-CCI(na.omit(HLC(Object1)))
  #myRSIValue<-myRSI(na.omit(Object1))
  #tail(CCIobject)
  #tail(myRSIValue)
  #if(Object[nrow(Object),5]>100000 && cond7>=1.1 && cond7<=1.5){
  #if(Object[nrow(Object),5]>100000 && tail(CCIobject,n=1)[1]< -95 && tail(myRSIValue,n=1)[1]<45){
  print(c(i,symbol))
  file=paste(symbol,"png",sep=".")
  png(file,width = 8, height = 8, units = 'in', res = 1200)
  chartSeries(Object,subset='last 24 months',main=symbol,TA=c(addVo(),addMACD(),addBBands(),addCCI(),addRSI(n = 12, wilder = TRUE),addEMA(n=60, col="blue"),addEMA(n=120, col="yellow")))
  write.table(symbol,file="mycadidate.bb.txt",sep="\t",col.names=F,row.names=F,quote=F,append = T)
  dev.off()
} 


