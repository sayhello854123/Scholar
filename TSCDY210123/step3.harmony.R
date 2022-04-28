rm(list = ls())
options(stringsAsFactors = F)
gc()
library(future)
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

load('./00.data/single_data/scRNA_qc.Rdata')
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters(resolution = 0.1)
p1 <- DimPlot(scRNA, group.by = "Sample",reduction = 'tsne')
p2 <- DimPlot(scRNA, group.by = "Sample",reduction = 'umap')
p3<- DimPlot(scRNA, group.by = "Sample", reduction = 'tsne',split.by = "Sample", ncol = 2)
p4<- DimPlot(scRNA, group.by = "Sample", reduction = 'umap',split.by = "Sample", ncol = 2)
cellinfo <- scRNA@meta.data
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
p5 <- DimPlot(scRNA, group.by = "Sample",reduction = 'tsne')
p6 <- DimPlot(scRNA, group.by = "Sample",reduction = 'umap')
p7 <- p1+p5
p8 <- p2+p6
p9 <- p5/p6+plot_layout(guides = 'collect')
P10 <- DimPlot(scRNA, group.by = "Sample", reduction = 'tsne',split.by = "Sample", ncol = 2)
P11 <- DimPlot(scRNA, group.by = "Sample", reduction = 'umap',split.by = "Sample", ncol = 2)
p12 <- (p3|p4)/(P10|P11)
dir.create('./03.Harmony')
ggsave('./03.Harmony/harmony_1.pdf',p9,width =8,height = 8,dpi = 1000 )
ggsave('./03.Harmony/harmony_2.pdf',p12,width =12,height = 8,dpi = 1000 )
system.time(save(scRNA, file = "./00.data/single_data/scRNA_harmony.Rdata"))  

mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","Sample"))
allcolour=brewer.pal(11,'Paired')
allcolour[1:2] <- c('#FF6A6A','darkolivegreen2')
p13 <- ggplot(mydata,aes(x= UMAP_1 , y = UMAP_2 ,color = Sample)) +geom_point(size = 1 , alpha =1 )+
  scale_color_manual(values = allcolour)+theme(panel.grid.major = element_blank(), #主网格线
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
  geom_segment(aes(x = min(mydata$UMAP_1) , y = min(mydata$UMAP_2) ,
                   xend = min(mydata$UMAP_1) +2.5, yend = min(mydata$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(mydata$UMAP_1)  , y = min(mydata$UMAP_2)  ,
                   xend = min(mydata$UMAP_1) , yend = min(mydata$UMAP_2) + 2.5),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(mydata$UMAP_1) +1, y = min(mydata$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(mydata$UMAP_1) -1, y = min(mydata$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
ggsave('./03.Harmony/harmony_xx.pdf',p13,width =12,height = 8,dpi = 1000 )
ggsave('./03.Harmony/harmony_xx.jpg',p13,width =12,height = 8,dpi = 1000 )
