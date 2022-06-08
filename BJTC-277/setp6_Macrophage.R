rm(list = ls())
options(stringsAsFactors = F) 
gc()
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)

load('./00.data/01.single_data/scRNA_bacth_celltype.Rdata')
Idents(scRNA) <- 'celltype'
scRNA <- subset(scRNA, idents = c("Monocyte/Macrophage"))
cellinfo <- scRNA@meta.data
table(cellinfo$Patients_ID)
# Monocyte/Macrophage
# 1358
# Li1 Li2 LN1 LN2  P1 PT1 PT2 PT3 
# 52  66  17  79 205 147 723  69 
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
ElbowPlot(scRNA,ndims = 50)
pc.num = 1:40
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims = pc.num) %>% FindClusters()
###diff clusters
diff.mast = FindAllMarkers(scRNA,test.use ="MAST")
all.markers = diff.mast %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top10 <- all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(all.markers, "./05.macrophage/diff_genes_mast.csv", row.names = F)
write.csv(top10, "./05.macrophage/top10_diff_genes_mast.csv", row.names = F)

refdata <- celldex::HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,file = "./05.macrophage/type.csv")

genes_to_check <- c("S100A9","S100A8",##Monocyte
                    "CD1C","LAMP3",#Dendritic cell
                    "CD68","CD163",#Macrophage
                    "TNF",##M1
                    "MRC1"##M2
)
 DotPlot(scRNA, features = unique(genes_to_check),
                    assay='RNA')  + coord_flip()+
             scale_colour_viridis()+RotatedAxis()
celltype <- read.csv('./05.macrophage/type.csv')
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
SaveH5Seurat(scRNA,'./00.data/01.single_data/seurat_macrophage.h5seurat',overwrite=T)
Convert('./00.data/01.single_data/seurat_macrophage.h5seurat','./00.data/01.single_data/seurat_macrophage.h5ad',overwrite=T)

p1 <- DotPlot(scRNA, features = unique(genes_to_check),
              assay='RNA',group.by ="celltype" )  + coord_flip()+
  scale_colour_viridis()+RotatedAxis()
ggsave('./05.macrophage/cell_DotPlot.pdf',p1,width =12,height = 8,dpi = 1000 )
ggsave('./05.macrophage/cell_DotPlot.jpg',p1,width =12,height = 8,dpi = 1000 )

mydata<- FetchData(scRNA,vars = c("tSNE_1","tSNE_2","celltype"))
allcolour <- c(brewer.pal(5, "Set1"))
if(T){
  p2 <- ggplot(mydata,aes(x= tSNE_1 , y = tSNE_2 ,color = celltype)) +geom_point(size = 1 , alpha =1 )+
    scale_color_manual(values = allcolour)+
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线 
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white"))+
    theme(legend.title = element_blank(), #去掉legend.title 
          legend.key=element_rect(fill='white'), #
          legend.text = element_text(size=20), #设置legend标签的大小
          legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=3.9)))+ #设置legend中 点的大小 
    geom_segment(aes(x = min(mydata$tSNE_1) , y = min(mydata$tSNE_2) ,
                     xend = min(mydata$tSNE_1) +6.5, yend = min(mydata$tSNE_2) ),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
    geom_segment(aes(x = min(mydata$tSNE_1)  , y = min(mydata$tSNE_2)  ,
                     xend = min(mydata$tSNE_1) , yend = min(mydata$tSNE_2) + 6.5),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(mydata$tSNE_1) +3, y = min(mydata$tSNE_2) -1.5, label = "tSNE_1",
             color="black",size = 3, fontface="bold" ) + 
    annotate("text", x = min(mydata$tSNE_1) -1.5, y = min(mydata$tSNE_2) + 3, label = "tSNE_2",
             color="black",size = 3, fontface="bold" ,angle=90) 
  ggsave('./05.macrophage/cell_type.pdf',p2,width =12,height = 8,dpi = 1000 )
  ggsave('./05.macrophage/cell_type.jpg',p2,width =12,height = 8,dpi = 1000 )
}

p3 <- FeaturePlot(scRNA,features = unique(genes_to_check), coord.fixed = T, order = T,
                  reduction = "tsne",cols = viridis(10),ncol = 4)
ggsave('./05.macrophage/cell_Feature.pdf',p3,width =24,height = 8,dpi = 600 )
ggsave('./05.macrophage/cell_Feature.jpg',p3,width =24,height = 8,dpi = 600)
##批量feature
if(T){
  dir.create('./05.macrophage/feature/')
  for(i in genes_to_check){
    p4 <- FeaturePlot(scRNA,features = unique(i), coord.fixed = T, order = T,
                      reduction = "tsne",cols = viridis(10))
    filenames1 <- paste0('./05.macrophage/feature/',i,'.pdf')
    filenames2 <- paste0('./05.macrophage/feature/',i,'.jpg')
    ggsave(filename = filenames1,p4,width =8,height = 8,dpi = 1000)
    ggsave(filename = filenames2,p4,width =8,height = 8,dpi = 1000)
  }
}

system.time(save(scRNA, file = "./00.data/01.single_data/scRNA_macrophage.Rdata"))
