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

load("./00.data/01.single_data/scRNA_bath.Rdata")
dir.create('./03.celltype/02.bath/')
Idents(scRNA) <- 'seurat_clusters'
###diff clusters
diff.mast = FindAllMarkers(scRNA,test.use ="MAST")
all.markers = diff.mast %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top10 <- all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.csv(all.markers, "./03.celltype/02.bath/diff_genes_mast.csv", row.names = F)
write.csv(top10, "./03.celltype/02.bath/top10_diff_genes_mast.csv", row.names = F)

###singleR
refdata <- celldex::HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,file = "./03.celltype/02.bath/type.csv")

genes_to_check <- c("EPCAM", #epithelial 
                    "PECAM1",	#stromal
                    "FN1",#fibroblast
                    "MKI67",#proliferative
                    "CD3D","CD2",##T
                    "MS4A1","CD79A",##B
                    "MZB1",#浆细胞
                    "KLRD1",##NK
                    "KIT",#mast cells 
                    "CSF3R","FCGR3A",#Neutrophils
                    "LILRA4",#pDC
                    "CD163","CD68"#monocytes/macrophage
                    
)#myeloid
# # # marker ann
# DotPlot(scRNA, features = unique(genes_to_check),
#         assay='RNA')  + coord_flip()+
#      scale_colour_viridis()+RotatedAxis()
celltype <- read.csv('./03.celltype/02.bath/type.csv')
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

p1 <- DotPlot(scRNA, features = unique(genes_to_check),
              assay='RNA',group.by ="celltype" )  + coord_flip()+
  scale_colour_viridis()+RotatedAxis()
ggsave('03.celltype/02.bath/cell_DotPlot.pdf',p1,width =12,height = 8,dpi = 1000 )
ggsave('03.celltype/02.bath/cell_DotPlot.jpg',p1,width =12,height = 8,dpi = 1000 )


mydata<- FetchData(scRNA,vars = c("tSNE_1","tSNE_2","celltype"))
allcolour <- c("#FF0066", "#33CC99","#669900",
               "#99FF66", "#386CB0", "#6699FF",
               "#0066CC", "#FF6600", "#FDC086",
               "#FF66FF", "#E7298A", "#FF0033")
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
  ggsave('03.celltype/02.bath/cell_type.pdf',p2,width =12,height = 8,dpi = 1000 )
  ggsave('03.celltype/02.bath/cell_type.jpg',p2,width =12,height = 8,dpi = 1000 )
}
p3 <- FeaturePlot(scRNA,features = unique(genes_to_check), coord.fixed = T, order = T,
                  reduction = "tsne",cols = viridis(10),ncol = 4)
ggsave('03.celltype/02.bath/cell_Feature.pdf',p3,width =24,height = 16,dpi = 600 )
ggsave('03.celltype/02.bath/cell_Feature.jpg',p3,width =24,height = 16,dpi = 600)
##批量feature
if(T){
  dir.create('./03.celltype/02.bath/feature/')
  for(i in genes_to_check){
    p4 <- FeaturePlot(scRNA,features = unique(i), coord.fixed = T, order = T,
                      reduction = "tsne",cols = viridis(10))
    filenames1 <- paste0('./03.celltype/02.bath/feature/',i,'.pdf')
    filenames2 <- paste0('./03.celltype/02.bath/feature/',i,'.jpg')
    ggsave(filename = filenames1,p4,width =8,height = 8,dpi = 1000)
    ggsave(filename = filenames2,p4,width =8,height = 8,dpi = 1000)
  }
}

p5 <- scRNA@meta.data %>% ggplot(aes(x = Patients_ID))+geom_bar(aes(fill = celltype))+
  scale_x_discrete("")+scale_y_continuous('cell number',expand = c(0.02,0))+
  scale_fill_manual(values = allcolour)+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.ticks.x  = element_blank(),
    #axis.text.x = element_blank(),
    #legend.position = "none"
  )
ggsave('03.celltype/02.bath/cell_barplot.pdf',p5,width =12,height = 8,dpi = 1000 )
ggsave('03.celltype/02.bath/cell_barplot.jpg',p5,width =12,height = 8,dpi = 1000)

p6 <- scRNA@meta.data %>% ggplot(aes(x = Organization))+geom_bar(aes(fill = celltype))+
  scale_x_discrete("")+scale_y_continuous('cell number',expand = c(0.02,0))+
  scale_fill_manual(values = allcolour)+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.ticks.x  = element_blank(),
    #axis.text.x = element_blank(),
    #legend.position = "none"
  )
ggsave('03.celltype/02.bath/cell_barplot1.pdf',p6,width =12,height = 8,dpi = 1000 )
ggsave('03.celltype/02.bath/cell_barplot1.jpg',p6,width =12,height = 8,dpi = 1000)

p7 <- scRNA@meta.data %>% ggplot(aes(x = Type))+geom_bar(aes(fill = celltype))+
  scale_x_discrete("")+scale_y_continuous('cell number',expand = c(0.02,0))+
  scale_fill_manual(values = allcolour)+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.ticks.x  = element_blank(),
    #axis.text.x = element_blank(),
    #legend.position = "none"
  )
ggsave('03.celltype/02.bath/cell_barplot2.pdf',p7,width =8,height = 8,dpi = 1000 )
ggsave('03.celltype/02.bath/cell_barplot2.jpg',p7,width =8,height = 8,dpi = 1000)

system.time(save(scRNA, file = "./00.data/01.single_data/scRNA_bacth_celltype.Rdata"))



A <- scRNA@meta.data
match_celltype_levels <- c("B cell", "Endothelial", "Epithelial", "Fibroblast",
                           "mast cell", "Monocyte/Macrophage", "Neutrophils","NK cell",
                           "pDC","Plasma cell","Proliferative","T cell")
a <- A%>%
  group_by(Patients_ID) %>%
  mutate(celltype = factor(celltype, levels = match_celltype_levels)) %>%
  arrange(celltype)
ggplot() +
  geom_bar(data = a, aes(x = Patients_ID, fill = factor(celltype)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = allcolour) +
  labs(fill = "cell type", y = "fraction of cells")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线 
        panel.border = element_blank(), #边框
        #axis.title = element_blank(),  #轴标题
        #axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))

b <- A%>%
  group_by(Organization) %>%
  mutate(celltype = factor(celltype, levels = match_celltype_levels)) %>%
  arrange(celltype)
p15 <- ggplot() +
  geom_bar(data = a, aes(x = Organization, fill = factor(celltype)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = allcolour) +
  labs(fill = "cell type", y = "fraction of cells")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white',colour = "black",size = 1), #背景色
        plot.background=element_rect(fill="white")
)
ggsave('03.celltype/02.bath/cell_barplot2.pdf',p15,width = 12,height = 8,dpi = 1000)

c <- A%>%
  group_by(Type) %>%
  mutate(celltype = factor(celltype, levels = match_celltype_levels)) %>%
  arrange(celltype)
p16 <- ggplot() +
  geom_bar(data = c, aes(x = Type, fill = factor(celltype)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = allcolour) +
  labs(fill = "cell type", y = "fraction of cells")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white',colour = "black",size = 1), #背景色
        plot.background=element_rect(fill="white")
  )
ggsave('03.celltype/02.bath/cell_barplot3.pdf',p16,width = 12,height = 8,dpi = 1000)

p17 <- p15/p16+plot_layout(guides = 'collect')
ggsave('03.celltype/02.bath/cell_barplot4.pdf',p17,width = 12,height = 16,dpi = 1000)
ggsave('03.celltype/02.bath/cell_barplot4.jpg',p17,width = 12,height = 16,dpi = 1000)
