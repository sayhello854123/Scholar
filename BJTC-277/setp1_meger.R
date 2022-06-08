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
library(SeuratDisk)

file <- list.dirs('./00.data/single_data')
file <- file[2:9]
metadata <- read.table('./00.data/metadata.txt',header = T,sep = '\t',
                       check.names = F,row.names = 1)
scRNAlist <- list()
for(i in 1:length(file)){
  project <- str_split((str_split(file[i],'_')[[1]][2]),'/')[[1]][2]
  scRNAlist[[i]] <- CreateSeuratObject(counts = Read10X(file[i]), project=project,
                     min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project) 
  ##添加样本信息
  meta <- metadata[i,]
  #scRNAlist[[i]] <- AddMetaData(object = scRNAlist[[i]], metadata = meta)
  a <- length(colnames(scRNAlist[[i]]))
  Organization <- rep(x =meta[1,2],a)
  Patients_ID <- rep(x =meta[1,1],a)
  Type <- rep(x =meta[1,3],a)
  scRNAlist[[i]] $Organization <- Organization
  scRNAlist[[i]] $Patients_ID <- Patients_ID
  scRNAlist[[i]] $Type <- Type
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  }
  #计算核糖体基因比例
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  #计算红细胞基因比例
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
}
### 给列表命名并保存数据
dir.create("./00.data/01.single_data")
files <- sapply(str_split(file,'_'),'[',2)
samples_name <- sapply(str_split(files,'/'),'[',2)
names(scRNAlist) <- samples_name
system.time(save(scRNAlist, file = "./00.data/01.single_data/scRNAlist.Rdata")) 

##数据合并
scRNA=merge(x=scRNAlist[[1]],
            scRNAlist[2:length(scRNAlist)])
table(scRNA$Patients_ID)          
system.time(save(scRNA, file = "./00.data/01.single_data/scRNA_orig.Rdata"))  

##转换scanpy
SaveH5Seurat(scRNA,'./00.data/01.single_data/seurat_raw.h5seurat',overwrite=T)
Convert('./00.data/01.single_data/seurat_raw.h5seurat','./00.data/01.single_data/seurat_raw.h5ad',overwrite=T)

##10x style
file <- list.dirs('./00.data/single_data')
file <- file[2:10]
metadata <- read.table('./00.data/metadata.txt',header = T,sep = '\t',
                       check.names = F,row.names = 1)
scRNAlist <- list()
for(i in 1:length(file)){
  project <- str_split((str_split(file[i],'_')[[1]][2]),'/')[[1]][2]
  scRNAlist[[i]] <- CreateSeuratObject(counts = Read10X(file[i]), project=project
  )
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project) 
  ##添加样本信息
  meta <- metadata[i,]
  #scRNAlist[[i]] <- AddMetaData(object = scRNAlist[[i]], metadata = meta)
  a <- length(colnames(scRNAlist[[i]]))
  Organization <- rep(x =meta[1,2],a)
  Patients_ID <- rep(x =meta[1,1],a)
  Type <- rep(x =meta[1,3],a)
  scRNAlist[[i]] $Organization <- Organization
  scRNAlist[[i]] $Patients_ID <- Patients_ID
  scRNAlist[[i]] $Type <- Type
}

##数据合并
scRNA=merge(x=scRNAlist[[1]],
            scRNAlist[2:length(scRNAlist)])
table(scRNA$Patients_ID)          
##To 10x
dir.10x = './00.data/seurat_10x/'
mtx = scRNA@assays$RNA@counts
Matrix::writeMM(mtx, paste0(dir.10x, 'matrix.mtx'))
write.table(rownames(scRNA@assays$RNA@counts), paste0(dir.10x, 'genes.tsv'), sep='\t', row.names = TRUE, col.names = FALSE)
write.table(colnames(scRNA@assays$RNA@counts), paste0(dir.10x, 'barcodes.tsv'), sep='\t', row.names = TRUE, col.names = FALSE)
metadata <- scRNA@meta.data
write.csv(metadata,file = 'GSE163558_meta.csv')


load('00.data/01.single_data/scRNA_orig.Rdata')
length(colnames(scRNA))#45388
length(rownames(scRNA))#25366
theme.set2 = theme(axis.title.x=element_blank())
# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
group = "Patients_ID"
# 质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)    
dir.create("./01.QC")
ggsave("./01.QC/vlnplot_before_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("./01.QC/vlnplot_before_qc.jpg", plot = violin, width = 18, height = 6,dpi = 1000)


### 设置质控标准
minGene=quantile(scRNA$nFeature_RNA,.02)#316
maxGene=quantile(scRNA$nFeature_RNA,.98)#5952
maxUMI=quantile(scRNA$nCount_RNA,.95)#20752
pctMT=20


### 数据质控并绘制小提琴图
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT )
length(colnames(scRNA))#38148
length(rownames(scRNA))#25366
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)      
ggsave("./01.QC/vlnplot_after_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("./01.QC/vlnplot_after_qc.jpg", plot = violin, width = 18, height = 6,dpi = 1000)
system.time(save(scRNA, file = "./00.data/01.single_data/scRNA_qc.Rdata"))  

