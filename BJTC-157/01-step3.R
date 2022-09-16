library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
rm(list=ls())
sce.all.filt <- readRDS("./00.data/01.single cell/scRNA2.rds")


#tSNE
pc.num=1:20
sce.all.filt = RunTSNE(sce.all.filt, dims = pc.num)
embed_tsne <- Embeddings(sce.all.filt, 'tsne')
DimPlot(sce.all.filt, reduction = "tsne",group.by = 'Sample')
tx = sce.all.filt@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx = sce.all.filt@meta.data$Sample)
if(T){
p=ggplot(tx, aes(x=tSNE_1, y=tSNE_2, color=tx)) + geom_point(size=1.8) + 
  scale_color_manual(values=c("H18" = "darkorange1", "H21" = "goldenrod4", 
                              "H23" = "blue1", "H28" = "dodgerblue", 
                              "H30" = "darkviolet", "H37" = "chartreuse2",
                              "H38" = "firebrick1"))
theme= theme(panel.grid =element_blank()) +   ## 删去网格
  theme(panel.border = element_blank(),panel.background =  element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black"))
p1=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("./01.single cell/tSNE.pdf", plot = p1, width = 8, height = 6) 
#ggsave("./01.single cell/tSNE.png", plot = p, width = 8, height = 6)
}

#UMAP
sce.all.filt <- RunUMAP(sce.all.filt, dims = pc.num)
embed_umap <- Embeddings(sce.all.filt, 'umap')
DimPlot(sce.all.filt, reduction = "umap",group.by = 'Sample') 
tx1 = sce.all.filt@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx1 = sce.all.filt@meta.data$Sample)
if(T){
  p1=ggplot(tx1, aes(x=UMAP_1, y=UMAP_2, color=tx1)) + geom_point(size=1.8) + 
    scale_color_manual(values=c("H18" = "darkorange1", "H21" = "goldenrod4", 
                                "H23" = "blue1", "H28" = "dodgerblue", 
                                "H30" = "darkviolet", "H37" = "chartreuse2",
                                "H38" = "firebrick1"))
  theme= theme(panel.grid =element_blank()) +   ## 删去网格
    theme(panel.border = element_blank(),panel.background =  element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))
  p1=p1+theme+guides(colour = guide_legend(override.aes = list(size=5)))
  ggsave("./01.single cell/UMAP.pdf", plot = p1, width = 8, height = 6) 
}
saveRDS(sce.all.filt, "./00.data/01.single cell/scRNA3.rds")
