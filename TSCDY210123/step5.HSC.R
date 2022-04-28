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
scRNA <- subset(scRNA, idents = c("Hepatic stellate cell "))
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
#DoHeatmap(scRNA,features = as.character(unique(top10$gene)),group.by = "Sample",)+NoLegend()+scale_fill_viridis()

dir.create('./05.hsc')
celltype_marker=c(
"Gata3","Hmgb1",#项目需要的基因
'Atg5','Pten'#自噬
)
p1 <- FeaturePlot(scRNA,features = c("Gata3","Hmgb1"), coord.fixed = T, order = T,
            reduction = "umap",cols = viridis(10),split.by ='Sample')

ggsave('./05.hsc/first.pdf',p1,width = 8,height = 8,dpi = 1000)
ggsave('./05.hsc/first.jpg',p1,width = 8,height = 8,dpi = 1000)

p2 <- FeaturePlot(scRNA,features = c('Atg5','Pten'), coord.fixed = T, order = T,
                  reduction = "umap",cols = viridis(10),split.by ='Sample')
ggsave('./05.hsc/second.pdf',p2,width = 8,height = 8,dpi = 1000)
ggsave('./05.hsc/second.jpg',p2,width = 8,height = 8,dpi = 1000)

p3 <- FeaturePlot(scRNA,features = "Gata3", coord.fixed = T, order = T,
                  reduction = "umap",cols = viridis(10),split.by ='Sample')
ggsave('./05.hsc/Gata3.pdf',p3,width = 8,height = 6,dpi = 1000)
ggsave('./05.hsc/Gata3.jpg',p3,width = 8,height = 6,dpi = 1000)

p4 <- FeaturePlot(scRNA,features = "Hmgb1", coord.fixed = T, order = T,
                  reduction = "umap",cols = viridis(10),split.by ='Sample')
ggsave('./05.hsc/Hmgb1.pdf',p4,width = 8,height = 6,dpi = 1000)
ggsave('./05.hsc/Hmgb1.jpg',p4,width = 8,height = 6,dpi = 1000)

mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","Sample"))
allcolour=c("#FF0066",'#6633FF')

p1 <- ggplot(mydata,aes(x= UMAP_1 , y = UMAP_2 ,color = Sample)) +geom_point(size = 1 , alpha =1 )+
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
ggsave('./05.hsc/cell.pdf',p1,width =12,height = 8,dpi = 1000 )
ggsave('./05.hsc/cell.jpg',p1,width =12,height = 8,dpi = 1000 )
