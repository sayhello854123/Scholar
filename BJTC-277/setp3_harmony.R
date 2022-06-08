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

load('./00.data/01.single_data/scRNA_qc.Rdata')
scRNA <- SCTransform(scRNA,vars.to.regress = 'percent.mt',variable.features.n = 3000)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
#Harmony整合
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
pc.num=1:36
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num)
scRNA <- FindNeighbors(scRNA, reduction = "harmony",
                       dims=pc.num) 
scRNA <- FindClusters(scRNA,resolution = 0.5)
dir.create('./02.harmony')
##查看分类
Idents(scRNA) <- 'seurat_clusters'
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","seurat_clusters"))
cell_type_med <- mydata %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
p2 <- ggplot(mydata, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values=col)+
  guides(colour=guide_legend(override.aes=list(size=10)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white")
  )+geom_label_repel(aes(label=seurat_clusters), fontface="bold",data = cell_type_med,
                     point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none")
ggsave('./02.harmony/1.UMAP/cluster.pdf',p2,width = 8,height = 6,dpi = 1000)
ggsave('./02.harmony/1.UMAP/cluster.jpg',p2,width = 8,height = 6,dpi = 1000)

##查看患者
Idents(scRNA) <- 'Patients_ID'
col <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
         "#FF9900","#FF0066","#A65628","#F781BF")
mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","Patients_ID"))
cell_type_med <- mydata %>%
  group_by(Patients_ID) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
p3 <- ggplot(mydata, aes(x = UMAP_1, y = UMAP_2, color = Patients_ID)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values=col)+
  guides(colour=guide_legend(override.aes=list(size=10)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white")
  )+geom_label_repel(aes(label=Patients_ID), fontface="bold",data = cell_type_med,
                     point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")
ggsave('./02.harmony/1.UMAP/Patients.pdf',p3,width = 8,height = 6,dpi = 1000)
ggsave('./02.harmony/1.UMAP/Patients.jpg',p3,width = 8,height = 6,dpi = 1000)

##查看原发与转移
Idents(scRNA) <- 'Type'
col <- c('#FF6600',"#FF0000")
mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2",'Type'))
cell_type_med <- mydata %>%
  group_by(Type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
p4 <- ggplot(mydata, aes(x = UMAP_1, y = UMAP_2, color = Type)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values=col)+
  guides(colour=guide_legend(override.aes=list(size=10)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white")
  )+geom_label_repel(aes(label=Type), fontface="bold",data = cell_type_med,
                     point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")
ggsave('./02.harmony/1.UMAP/Type.pdf',p4,width = 8,height = 6,dpi = 1000)
ggsave('./02.harmony/1.UMAP/Type.jpg',p4,width = 8,height = 6,dpi = 1000)

P3 <- p2|p3|p4
ggsave('./02.harmony/1.UMAP/all1.pdf',P3,width = 24,height = 8,dpi = 1000)
ggsave('./02.harmony/1.UMAP/all1.jpg',P3,width = 24,height = 8,dpi = 1000)

load('./01.QC/all1.RData')
P4 <- P1/P3
ggsave('./02.harmony/1.UMAP/all.pdf',P4,width = 24,height = 16,dpi = 1000)
ggsave('./02.harmony/1.UMAP/all.jpg',P4,width = 24,height = 16,dpi = 1000)


#查看分簇
Idents(scRNA) <- 'seurat_clusters'
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mydata<- FetchData(scRNA,vars = c("tSNE_1","tSNE_2","seurat_clusters"))
cell_type_med <- mydata %>%
  group_by(seurat_clusters) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
p5 <- ggplot(mydata, aes(x = tSNE_1, y = tSNE_2, color = seurat_clusters)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values=col)+
  guides(colour=guide_legend(override.aes=list(size=10)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white")
  )+geom_label_repel(aes(label=seurat_clusters), fontface="bold",data = cell_type_med,
                     point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none")
ggsave('./02.harmony/2.TSNE/cluster.pdf',p5,width = 8,height = 6,dpi = 1000)
ggsave('./02.harmony/2.TSNE/cluster.jpg',p5,width = 8,height = 6,dpi = 1000)

##查看患者
Idents(scRNA) <- 'Patients_ID'
col <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
         "#FF9900","#FF0066","#A65628","#F781BF")
mydata<- FetchData(scRNA,vars = c("tSNE_1","tSNE_2","Patients_ID"))
cell_type_med <- mydata %>%
  group_by(Patients_ID) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
p6 <- ggplot(mydata, aes(x = tSNE_1, y = tSNE_2, color = Patients_ID)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values=col)+
  guides(colour=guide_legend(override.aes=list(size=10)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white")
  )+geom_label_repel(aes(label=Patients_ID), fontface="bold",data = cell_type_med,
                     point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")
ggsave('./02.harmony/2.TSNE/Patients.pdf',p6,width = 8,height = 6,dpi = 1000)
ggsave('./02.harmony/2.TSNE/Patients.jpg',p6,width = 8,height = 6,dpi = 1000)


##查看原发与转移
Idents(scRNA) <- 'Type'
col <- c('#FF6600',"#FF0000")
mydata<- FetchData(scRNA,vars = c("tSNE_1","tSNE_2",'Type'))
cell_type_med <- mydata %>%
  group_by(Type) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
p7 <- ggplot(mydata, aes(x = tSNE_1, y = tSNE_2, color = Type)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values=col)+
  guides(colour=guide_legend(override.aes=list(size=10)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white")
  )+geom_label_repel(aes(label=Type), fontface="bold",data = cell_type_med,
                     point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")
ggsave('./02.harmony/2.TSNE/Type.pdf',p7,width = 8,height = 6,dpi = 1000)
ggsave('./02.harmony/2.TSNE/Type.jpg',p7,width = 8,height = 6,dpi = 1000)

P5 <- p5|p6|p7
ggsave('./02.harmony/2.TSNE/all2.pdf',P5,width = 24,height = 8,dpi = 1000)
ggsave('./02.harmony/2.TSNE/all2.jpg',P5,width = 24,height = 8,dpi = 1000)

load('./01.QC/all2.RData')
P6 <- P2/P5
ggsave('./02.harmony/2.TSNE/all.pdf',P6,width = 24,height = 16,dpi = 1000)
ggsave('./02.harmony/2.TSNE/all.jpg',P6,width = 24,height = 16,dpi = 1000)

system.time(save(scRNA, file = "./00.data/01.single_data/scRNA_harmony.Rdata"))
