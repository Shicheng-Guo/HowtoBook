# Factor or not in lm and glm
# for gender (binary), no difference beween factor or not
# However, factor will take a little more time

test1<-function(x){
x1<-rnorm(150)
x2<-rep(c(1,2),75)
y<-rnorm(150)
lm(y~x1+x2)
}

test2<-function(x){
  x1<-rnorm(150)
  x2<-rep(c(1,2),75)
  y<-rnorm(150)
  lm(y~x1+as.factor(x2))
}

start_time <- Sys.time()
for(i in 1:1000){
  time(test1(1))
}
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
for(i in 1:1000){
  time(test2(1))
}
end_time <- Sys.time()
end_time - start_time
