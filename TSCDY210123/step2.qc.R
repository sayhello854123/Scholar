rm(list = ls())
options(stringsAsFactors = F) 
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)

load('./00.data/single_data/scRNA_orig.Rdata')
length(colnames(scRNA))#7727
length(rownames(scRNA))#17224
theme.set2 = theme(axis.title.x=element_blank())
# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb")
group = "Sample"
# 质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 4)    
dir.create("./01.QC")
ggsave("./01.QC/vlnplot_before_qc.pdf", plot = violin, width = 24, height = 12)
ggsave("./01.QC/vlnplot_before_qc.jpg", plot = violin, width = 24, height = 12,dpi = 1000)

### 设置质控标准
### 设置质控标准
minGene=500
maxGene=6000
maxUMI=30000
pctMT=20
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT )


### 数据质控并绘制小提琴图
length(colnames(scRNA))#4071
length(rownames(scRNA))#17224
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 4)     
ggsave("./01.QC/vlnplot_after_qc.pdf", plot = violin, width = 24, height = 12)
ggsave("./01.QC/vlnplot_after_qc.jpg", plot = violin, width = 24, height = 12,dpi = 1000)
system.time(save(scRNA, file = "./00.data/single_data/scRNA_qc.Rdata"))  
