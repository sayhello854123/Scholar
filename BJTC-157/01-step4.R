library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(stringr) 
rm(list=ls())
sce<- readRDS("./00.data/01.single cell/scRNA2.rds")


#tSNE
pc.num=1:20
sce = RunTSNE(sce, dims = pc.num)
embed_tsne <- Embeddings(sce, 'tsne')

table(sce@meta.data$Type)
immune <- c("B cell","T cell","TAM")
sce@meta.data$immune_annotation <-ifelse(sce@meta.data$Type %in% immune ,'immune','non-immune')
table(sce@meta.data$immune_annotation)
tx = sce@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% cbind(group = sce@meta.data$immune_annotation)
if(T){
  p=ggplot(tx, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(size=1.8) + 
    scale_color_manual(values=c("non-immune" = "#48D1CC", "immune" = "#F08080" 
                                ))
  theme= theme(panel.grid =element_blank()) +   ## 删去网格
    theme(panel.border = element_blank(),panel.background =  element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))
  p1=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
  ggsave("./01.single cell/immune_group.pdf", plot = p1, width = 8, height = 6) 
  #ggsave("./01.single cell/tSNE.png", plot = p, width = 8, height = 6)
}
##差异表达
diff.wilcox <- FindAllMarkers(sce)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv")
write.csv(top10, "top10_diff_genes_wilcox.csv")

tx = sce@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% cbind(group = sce@meta.data$Type)
if(T){
  p=ggplot(tx, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(size=1.8) + 
    scale_color_manual(values=c("H18" = "darkorange1", "H21" = "goldenrod4", 
                                "H23" = "blue1", "H28" = "dodgerblue", 
                                "H30" = "darkviolet", "H37" = "chartreuse2",
                                "H38" = "firebrick1"))
  theme= theme(panel.grid =element_blank()) +   ## 删去网格
    theme(panel.border = element_blank(),panel.background =  element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))
  p2=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
  ggsave("./01.single cell/celltype.pdf", plot = p2, width = 8, height = 6) 
  #ggsave("./01.single cell/tSNE.png", plot = p, width = 8, height = 6)
}

saveRDS(sce, "./00.data/01.single cell/scRNA4.rds")
