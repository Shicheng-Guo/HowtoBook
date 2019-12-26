
GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                           stats = gene_list,
                           minSize=15,
                           maxSize=600,
                           nperm=10000) %>% 
                  as.data.frame() %>% 
                  dplyr::filter(padj < !!pval)
  print(dim(fgRes))
    
## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  print(dim(rbind(ups,downs)))
  
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  ### Collapse redundant pathways
  Up = fgsea::collapsePathways(keepups, pathways = myGO, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns, myGO, gene_list,  nperm = 500, pval.threshold = 0.05) 
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
           c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
	
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")

  filtRes = rbind(head(fgRes, n = 10),tail(fgRes, n = 10 ))
  
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
  geom_point(size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                      "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()

  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}


library(dplyr)
library(gage)
library(fgsea)

S4table = read.csv("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/GSEA/journal.pone.0145322.s007.csv", header=TRUE, skip =1) %>%  filter(Gene.Symbol != "")
gene_list = S4table$DESeq2.Log2.Fold.Change
names(gene_list) = S4table$Gene.Symbol
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)

# http://software.broadinstitute.org/gsea/downloads.jsp
system("wget http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/msigdb.v7.0.symbols.gmt -O msigdb.v7.0.symbols.gmt")
system("wget http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/c5.bp.v7.0.symbols.gmt -O c5.bp.v7.0.symbols.gmt")

GO_file = "msigdb.v7.0.symbols.gmt"
GO_file = "c5.bp.v7.0.symbols.gmt"
GO_file = "c2.cp.kegg.v7.0.symbols.gmt"

pval=0.05
res = GSEA(gene_list=gene_list,GO_file=GO_file, pval = 0.05)
ggsave(rlt$Plot,file="gsea.jpg")
png("gsea.png")
res$Plot
dev.off()
