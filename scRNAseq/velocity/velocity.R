library(velocyto.R)

ldat <- read.loom.matrices("data/velocyto/onefilepercell_A1_unique_and_others_J2CH1.loom")
ldat <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/ldat.rds"))
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("_unique.bam","",gsub(".*:","",colnames(x)))
  x
})
cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))
emb <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/embedding.rds"))
hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')

# exonic read (spliced) expression matrix
emat <- ldat$spliced;
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning;
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))
length(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))
fit.quantile <- 0.05;
rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = fit.quantile)
pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1))
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Chga',cell.emb=emb,cell.colors=cell.colors)
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile=fit.quantile, diagonal.quantiles = TRUE)
pca.velocity.plot(rvel,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,1,1,1))
rvel1 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,deltaT2 = 1,kCells = 1, fit.quantile=fit.quantile)
pca.velocity.plot(rvel1,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,1,1,1))
show.velocity.on.embedding.cor(emb,vel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)
show.velocity.on.embedding.cor(emb,vel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=20,arrow.lwd=2)



require(BSgenome.Mmusculus.UCSC.mm10)
ip.mm10 <- find.ip.sites('data/genes.gtf',Mmusculus,'mm10')
ip.mm10 <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/ip.mm10.rds"))
gene.info <- read.gene.mapping.info("data/velocyto/dump/onefilepercell_A1_unique_and_others_J2CH1.hdf5",internal.priming.info=ip.mm10,min.exon.count=10);
gene.info <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/gene.info.rds"))
emat <- ldat$spliced; nmat <- ldat$unspliced; smat <- ldat$spanning
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 7)
gvel <- global.velcoity.estimates(emat, nmat, rvel, base.df=gene.info$gene.df, smat=smat, deltaT=1, kCells=5, kGenes = 15, kGenes.trim = 5, min.gene.cells = 0, min.gene.conuts = 500)
pca.velocity.plot(gvel,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,1,1))par(mfrow=c(1,2), mar = c(2.5,2.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
arrow.scale=3; cell.alpha=0.4; cell.cex=1; fig.height=4; fig.width=4.5;
#pdf(file='tsne.rvel_gvel.plots.pdf',height=6,width=12)
show.velocity.on.embedding.cor(emb,rvel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1,main='gene-relative esitmate',do.par=F)
show.velocity.on.embedding.cor(emb,gvel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1,main='gene-structure estimate',do.par=F)
par(mfrow=c(1,2), mar = c(2.5,2.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
x <- tSNE.velocity.plot(rvel,nPcs=15,cell.colors=cell.colors,cex=0.9,perplexity=200,norm.nPcs=NA,pcount=0.1,scale='log',do.par=F,main='gene-relative estimate')
x <- tSNE.velocity.plot(gvel,nPcs=15,cell.colors=cell.colors,cex=0.9,perplexity=200,norm.nPcs=NA,pcount=0.1,scale='log',do.par=F,main='gene-structure estimate')
x <- show.velocity.on.embedding.eu(emb,rvel,n=40,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,nPcs=30,sigma=2.5,show.trajectories=TRUE,diffusion.steps=400,n.trajectory.clusters=15,ntop.trajectories=1,embedding.knn=T,control.for.neighborhood.density=TRUE,n.cores=40) 


