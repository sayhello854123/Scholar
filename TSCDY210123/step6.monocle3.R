rm(list = ls())
options(stringsAsFactors = F)
gc()
library(monocle3)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(future)
library(RColorBrewer)
library(viridis)
library(Seurat)

load('./00.data/single_data/scRNA_HSC.Rdata')

data <- GetAssayData(scRNA, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
#tsne降维
cds <- reduce_dimension(cds,preprocess_method = "PCA",reduction_method = 'UMAP')

##从seurat导入整合过的UMAP坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')

### Monocle3聚类分区
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")

### 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
##细胞按拟时排序
cds <- order_cells(cds) 
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
ggsave('./06.monocle3/m1.jpg',p1,width = 8,height = 6,dpi = 1000)
ggsave('./06.monocle3/m1.pdf',p1,width = 8,height = 6,dpi = 1000)

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
#Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "./06.monocle3/Trajectory_genes.csv", row.names = F)
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters", 
                         min_expr=0.5, ncol = 2)
p2 <- plot_genes_in_pseudotime(cds[celltype_marker,], color_cells_by="Sample", 
                         min_expr=0.5, ncol = 2)
ggsave('./06.monocle3/m2.jpg',p2,width = 12,height = 8,dpi = 1000)
ggsave('./06.monocle3/m2.pdf',p2,width = 12,height = 8,dpi = 1000)
