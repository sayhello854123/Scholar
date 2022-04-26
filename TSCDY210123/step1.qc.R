rm(list = ls())
options(stringsAsFactors = F)
library(devtools)
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(hdf5r)
library(data.table)
library(SeuratObject)

dir <- './00.data/GSE175939_RAW/'
samples=list.files(dir)
sampledir <- paste0(dir,samples)
metadata <- read.table('./00.data/metadata1.txt',header = T,sep = '\t',
                       check.names = F)
scRNAlist <- list()
for(i in 1:length(samples)){
  counts <- Read10X_h5(sampledir[i])
  project <- str_split(samples[i],'_')[[1]][1]
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=project,
                                       min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project)
  ##添加样本信息
  meta <- metadata[i,]
  scRNAlist[[i]] <- AddMetaData(object = scRNAlist[[i]], metadata = meta) 
  a <- length(colnames(scRNAlist[[i]]))
  #Patients_ID <- rep(x =meta[1,2],a)
  Sample <- rep(x =meta[1,2],a)
  #scRNAlist[[i]]$Patients_ID <- Patients_ID
  scRNAlist[[i]]$Sample <- Sample
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-") 
  }
  #计算核糖体基因比例
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rp[sl]")
  }
}
### 给列表命名并保存数据
samples_name <- sapply(strsplit(samples,"_"),'[',1)
names(scRNAlist) <- samples_name
system.time(save(scRNAlist, file = "./00.data/single_data/scRNAlist.Rdata"))  

##数据合并
scRNA=merge(x=scRNAlist[[1]],
            scRNAlist[2:length(scRNAlist)])
table(scRNA$Sample)          
system.time(save(scRNA, file = "./00.data/single_data/scRNA_orig.Rdata"))  
