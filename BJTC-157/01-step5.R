library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(stringr) 
rm(list=ls())
sce<- readRDS("./00.data/01.single cell/scRNA4.rds")

immune <- c("B cell","T cell","TAM")
sce@meta.data$immune_annotation <-ifelse(sce@meta.data$Type %in% immune ,'immune','non-immune')
phe=sce@meta.data
cells.use <- row.names(sce@meta.data)[which(phe$immune_annotation=='immune')]
sce <-subset(sce, cells=cells.use)  
tx = sce@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% cbind(group = sce@meta.data$Type)
if(T){
  p=ggplot(tx, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(size=1.8) + 
    scale_color_manual(values=c("B cell" = "firebrick1", "T cell" = "green4", 
                                "TAM" = "deepskyblue"))
  theme= theme(panel.grid =element_blank()) +   ## 删去网格
    theme(panel.border = element_blank(),panel.background =  element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))
  p3=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
  ggsave("./01.single cell/immune_cell.pdf", plot = p3, width = 8, height = 6) 
  #ggsave("./01.single cell/tSNE.png", plot = p, width = 8, height = 6)
}
diff.wilcox = FindAllMarkers(sce)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p4 <- DoHeatmap(sce,top10$gene,size=3,group.by = 'Type')
ggsave("./01.single cell/immune_heatmap.pdf", plot = p1, width = 8, height = 6) 
P <- (p1|p2)/(p3|p4)+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/picture.png", plot = P, width = 12, height = 7) 
saveRDS(sce, "./00.data/01.single cell/scRNA5.rds")
save(p1,p2,p3,p4,file = './00.data/01.single cell/plot1.RData')
library("jpeg")
library(ggpubr)
library(ggplot2)
A1 <- readJPEG('./01.single cell/celltype.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
