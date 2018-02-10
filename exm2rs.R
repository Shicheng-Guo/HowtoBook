setwd("/mnt/bigdata/Genetic/Projects/shg047/hemochromatosis/haplotype")
file=list.files(pattern="*sel2.map")
db<-read.table("/mnt/bigdata/Genetic/Projects/shg047/db/InfiniumCoreExome-24v1-2_A1_b144_rsids.txt",head=T)
for(i in 1:length(file)){
  map<-read.table(file[i],head=F,sep="\t")
  exm_id=unlist(lapply(map[,2],function(x) unlist(strsplit(as.character(x),split="_"))[1]))
  tmp<-cbind(map,db[match(exm_id,db[,1]),])
  write.table(tmp,paste(file[i],"map",sep="."),sep=" ",quote=F)
}
