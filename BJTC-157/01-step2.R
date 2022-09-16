library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())
sce.all.filt <- readRDS("./00.data/01.single cell/scRNA1.rds")
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
nExp <- round(ncol(sce.all.filt) * 0.04)
sce.all.filt <- doubletFinder_v3(sce.all.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)

#使用Seurat包中的FindNeighbors函数计算构建SNN图。
sce.all.filt=FindNeighbors(sce.all.filt, dims = 1:20, k.param = 60, prune.SNN = 1/15)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt=FindClusters(sce.all.filt, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
apply(sce.all.filt@meta.data[,grep("RNA_snn_res",colnames(sce.all.filt@meta.data))],2,table)
sce.all.filt <- FindClusters(sce.all.filt, resolution = 0.8)
table(sce.all.filt@meta.data$seurat_clusters)
metadata <- sce.all.filt@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
saveRDS(sce.all.filt, "./00.data/01.single cell/scRNA2.rds")

scRNA <- RunPCA(sce.all.filt, features = VariableFeatures(sce.all.filt),npcs = 100) 
scRNA <- JackStraw(scRNA, num.replicate = 100, dims = 100)
scRNA<- ScoreJackStraw(scRNA, dims = 1:100)
plot2 <- ElbowPlot(scRNA, ndims=30, reduction="pca") 
plot5 <- JackStrawPlot(scRNA, dims = 1:30)
plot <- plot2+plot5
ggsave("./01.single cell/01.Seurat/pca.pdf", plot = plot, width = 14, height = 6)
