rm(list = ls())
options(stringsAsFactors = F)
gc()
library(future)
library(ggpointdensity)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(SingleR)
library(MAST)

load('./00.data/single_data/scRNA_celltype.Rdata')
Idents(scRNA) <- 'celltype'
scRNA <- subset(scRNA, idents = c("Hepatic stellate cell"))
cellinfo <- scRNA@meta.data
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims = pc.num) %>% FindClusters()
system.time(save(scRNA, file = "./00.data/single_data/scRNA_HSC.Rdata"))  

DimPlot(scRNA,label = T)
Idents(scRNA)="Sample"
diff <- FindAllMarkers(scRNA,test.use = "MAST")
all.markers = diff %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top10 <- all.markers%>%group_by(cluster)%>%top_n(n=50,wt=avg_log2FC)
DoHeatmap(scRNA,features = top10$gene)+NoLegend()+scale_fill_viridis()
celltype_marker=c(
"Lpl","Plagl1",#项目需要的基因
  "Eln","Mfap2",'Atg5','Pten'#自噬
)
FeaturePlot(scRNA,features = unique(celltype_marker), coord.fixed = T, order = T,
            reduction = "umap",cols = viridis(10),split.by ='Sample')

data <- cbind(Embeddings(object=scRNA[['umap']]),FetchData(scRNA,'Sample'))
p <- ggplot(data = data, mapping = aes(x = UMAP_1,
                                       y = UMAP_2)) + 
  geom_pointdensity() + #密度散点图（geom_pointdensity）
  scale_color_viridis()+theme_bw()+facet_wrap(~Sample,ncol = 2)