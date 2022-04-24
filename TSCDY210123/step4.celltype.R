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

load('./00.data/single_data/scRNA_harmony.Rdata')

##鉴定簇
scRNA <- FindNeighbors(scRNA, reduction = "harmony",
              dims = 1:40) 

scRNA <- FindClusters(scRNA,resolution = 0.1)
Idents(scRNA) <- "seurat_clusters"
scRNA <- subset(scRNA, idents = c(0:8))
diff.mast = FindAllMarkers(scRNA,test.use ="MAST" )
top10 <- diff.mast%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)

##singleR鉴定
refdata <- MouseRNAseqData() 
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
DimPlot(scRNA, label = T)
write.csv(celltype,file = "type.csv")

##人工鉴定
celltype_marker=c(
"Dcn","Lrat",'Pdgfrb','Des',##HSC标志基因
'Gata3','Hmgb1',"Plrg1",#项目需要的基因
'Atg5','Pten','Tsc1',#自噬
"Cd68",##macrophage 
"Cd79a",##B细胞
"Ms4a4b",##T memory cells
"Alb","Apoa1",#Hepatocytes
"Pecam1"

)
FeaturePlot(scRNA,features = unique(celltype_marker), coord.fixed = T, order = T,
            reduction = "umap",cols = viridis(10))
genes_to_check <- c("Apoa1",#Hepatocytes
                    "Pecam1",	#Endothelial
                   "Dcn","Lrat",'Pdgfrb','Des',##HSC标志基因
                    "Cd68",##macrophage 
                   "Cd79a",##B细胞
                   "Ms4a4b"##T memory cells
                    )
celltype <- read.csv('./type.csv')
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

#DimPlot(scRNA,group.by = 'celltype',label=T, label.size=5, reduction='umap')
mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","celltype"))
allcolour=c("#FF6633","#663399","#FF3300","#663366","#FF9900","#66CC00")

p1 <- ggplot(mydata,aes(x= UMAP_1 , y = UMAP_2 ,color = celltype)) +geom_point(size = 1 , alpha =1 )+
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
ggsave('./04.cell_type/cell_type.pdf',p1,width =12,height = 8,dpi = 1000 )
ggsave('./04.cell_type/cell_type.jpg',p1,width =12,height = 8,dpi = 1000 )
p2 <- DotPlot(scRNA, features = unique(genes_to_check),
        assay='RNA',group.by ="celltype" )  + coord_flip()+
  scale_color_viridis()+RotatedAxis()
ggsave('./04.cell_type/cell_DotPlot.pdf',p2,width =12,height = 8,dpi = 1000 )
ggsave('./04.cell_type/cell_DotPlot.jpg',p2,width =12,height = 8,dpi = 1000 )



p3 <- FeaturePlot(scRNA,features = unique(genes_to_check), coord.fixed = T, order = T,
            reduction = "umap",cols = viridis(10),split.by ='Sample')
ggsave('./04.cell_type/cell_sample.pdf',p3,width =10,height = 24,dpi = 1000 )
p4 <- FeaturePlot(scRNA,features = unique(genes_to_check), coord.fixed = T, order = T,
            reduction = "umap",cols = viridis(10),ncol = 3)
ggsave('./04.cell_type/cell_marker.pdf',p4,width =16,height = 12,dpi = 1000 )
ggsave('./04.cell_type/cell_marker.jpg',p4,width =16,height = 12,dpi = 1000 )

p5=scRNA@meta.data %>% ggplot(aes(x = Sample))+geom_bar(aes(fill = celltype))+
scale_x_discrete("")+scale_y_continuous('cell number',expand = c(0.02,0))+
scale_fill_brewer(palette = "Set1")+
theme_bw()+theme(
panel.grid = element_blank(),
axis.ticks.x  = element_blank(),
#axis.text.x = element_blank(),
#legend.position = "none"
)

ggsave('./04.cell_type/cell_bar.jpg',p5,width =8,height = 6,dpi = 1000 )
ggsave('./04.cell_type/cell_bar.pdf',p5,width =8,height = 6,dpi = 1000 )
system.time(save(scRNA, file = "./00.data/single_data/scRNA_celltype.Rdata"))



