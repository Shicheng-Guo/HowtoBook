#################################
#     define functions          #
#################################
##!/bin/Rscript

#install.packages("Seurat","dplyr","cowplot")
library(Seurat)
library(dplyr)
library(cowplot)


#The integrate function. Change the dimensions as per need####
integrate<-function(data){
  immune.anchors <- FindIntegrationAnchors(object.list = data, dims = 1:30)
  immune.integrated <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
  DefaultAssay(immune.integrated) <- "integrated"
  immune.integrated <- ScaleData(immune.integrated, verbose = FALSE)
  immune.integrated <- RunPCA(immune.integrated, npcs = 30, verbose = FALSE)
  immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:30)
  immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:20)
  immune.integrated <- FindClusters(immune.integrated, resolution = 0.5)
  return(immune.integrated)
}

#Creating seurat objects from raw files. Change min # of RNA feature ###

seurat_object <-function(file_name,p.dir,s,g) {
  object<-Read10X(data.dir = paste0(p.dir,"/",file_name))
  colnames(x = object) <- paste(s, colnames(x = object), sep = '_')
  object <- CreateSeuratObject(object, project = "My_project", min.cells = 5)
  object@meta.data$object <- file_name
  object <- subset(object, subset = nFeature_RNA > 500)
  object <- NormalizeData(object)
  object <- ScaleData(object, display.progress = F)
  object <- FindVariableFeatures(object,selection.method = "vst", nfeatures = 2000)
  return(object)
}



my_Dimplot <-function(data){
  p1 <- DimPlot(data, reduction = "umap", group.by = "group")
  p2 <- DimPlot(data, reduction = "umap", label = TRUE)
  p123<-plot_grid(p1, p2)
  png(filename="Dim_PLot.png")
  return(plot(p123))
  dev.off()
}

my_dimheat<-function(data,markers){
  a2<-DoHeatmap(object=data,features=markers)
  png(filename="heatmap.png")
  return(plot(a2))
  dev.off()
}

my_vlnplot<-function(data,markers){
  
  a3<-VlnPlot(object = data, features = markers)
  png(filename="VlnplotCombinedall.png")
  return(plot(a3))
  dev.off()
}



my_tsne<-function(data){
  p1 <- TSNEPlot(data, pt.size = 0.5, group.by = "group")
  p2 <- TSNEPlot(data, pt.size = 0.8)
  p3<-plot_grid(p1, p2)
  
  png(filename="TSNE_Combinedall.png")
  return(plot(p3))
  dev.off()
}


##############################
#Script, change the PATH
##############################

files=list.files("Your path to folder with matrix files")
p.dir="Your path to folder with matrix files"

# you can create group label as per the variable names below, this list is manual as of now ###
#new.group <- c(CRERATE LABEL/GROUPS for each sample)
names<-files


## Calling create seurat object function on each file ####
j<-0
for(i in names){
  {
    j<-j+1
    assign(i,seurat_object(i,p.dir,paste0("S",j),paste0("mygene.",j))) 
  }
}


### assigning the group names created above, comment out if not needed ####
gp<-0
for (i in names) {
  gp<-gp+1
  Object = get(paste0(i))
  Object@meta.data$group=new.group[gp]
  assign(i,Object)
  print(new.group[gp])
}


##creating the list for integartion ####
Hall=NULL
for(i in names){
  Object = get(paste0(i))
  Hall <- append(Object, Hall)
}

###define you gene names for the markers ###

mymarkers =c("","") # fill your genes of interest

immune.combinedall <- integrate(Hall)
#my_Dimplot(immune.combinedall)
#my_dimheat(immune.combinedall,mymarkers)
#my_vlnplot(immune.combinedall,mymarkers)
#my_tsne(immune.combinedall)


#save.image(file="MyData.RData")

saveRDS(immune.combinedall,file="Integratedmodel_3_0.rds")













