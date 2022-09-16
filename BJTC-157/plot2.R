library(GSVA)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)
library(ggimage)
library(grid)
library(ggplot2)
library(patchwork)
library(EBImage)
library(imager)
library("jpeg")
library(ggpubr)
library(ComplexHeatmap)
cellMarker <- data.table::fread('./00.data/02.gene get/T Cell/all.txt',data.table = F)
colnames(cellMarker)[2] <- "celltype"
cellMarker <- lapply(split(cellMarker,cellMarker$celltype), function(x){
  dd = x$gene
  unique(dd)
})
sce<- readRDS("./00.data/01.single cell/scRNA5.rds")
table(sce@meta.data$Type)
sce@meta.data$T_cell <-ifelse(sce@meta.data$Type  %in% "T cell" ,'t cell','non-tcell')
table(sce@meta.data$T_cell)
phe=sce@meta.data
cells.use <- row.names(sce@meta.data)[which(phe$T_cell=='t cell')]
sce_t <-subset(sce, cells=cells.use)
a <- as.data.frame(sce_t@assays$RNA@scale.data)
expr <- as.matrix(a)
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
write.csv(gsva_data,file = '23.csv')
rownames(gsva_data)[6] <- c('Naive')
type <- as.data.frame(cbind(colnames(sce_t),sce_t@meta.data$Sample,sce_t@meta.data$Type))
rownames(type) <- type[,1]
type <- type [,-1]
colnames(type) <- c('Patient','Type')
write.csv(type ,file = '24.csv')
sameSample=intersect(colnames(gsva_data),row.names(type))
gsva_data=gsva_data[,sameSample]
type=type[sameSample,]
#Type=type[order(type[,var]),] 
rt=gsva_data[,row.names(type)]
col_fun = circlize::colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(
  Patient = anno_simple(type$Patient, col = col_fun),
  annotation_name_side = "left")
pheatmap(gsva_data,cluster_cols =FALSE,top_annotation = ha,  
show_colnames=F)
