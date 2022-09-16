library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
##step1
sce <- Read10X(data.dir = "./00.data/00.single Cell data/")
dim(sce)
metadata <- read.table('./00.data/00.single Cell data/GSE125449_Set1_samples.txt/GSE125449_Set1_samples.txt',
                       sep = '\t',header = T)
same <- intersect(metadata[,2],colnames(sce))
sce <- sce[,same]
rownames(metadata) <- metadata[,2]
sce.all <- CreateSeuratObject(counts = sce,  #和文章一样，初步质控
                             min.cells = 3, #在不少于三个细胞中表达
                             min.features = 200) #基因数不少于200
sce.all  <- AddMetaData(object = sce.all, metadata = metadata)
dim(sce.all)
#计算线粒体基因比例
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] 
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-.")
sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mt")
fivenum(sce.all@meta.data$percent_mt)
#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)

##可视化
VlnPlot(object = sce.all,features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)


#过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 700)
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]
sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all.filt) 

#过滤指标2:线粒体
selected_mito <- WhichCells(sce.all.filt, expression = percent_mt < 20)
sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
dim(sce.all.filt)
sce.all.filt<- NormalizeData(sce.all.filt, normalization.method = "LogNormalize", scale.factor = 10000)


sce.all.filt <- FindVariableFeatures(sce.all.filt, selection.method = "vst", nfeatures = 2000) 
##如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(sce.all.filt)
sce.all.filt <- ScaleData(sce.all.filt, features = scale.genes)

#细胞周期评分
##标准化数据
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(sce.all.filt))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(sce.all.filt))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(sce.all.filt))
sce.all.filt <- CellCycleScoring(object=sce.all.filt,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAa <- RunPCA(sce.all.filt, features = c(s_genes, g2m_genes))
DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
###### 检测doublets 
sce.all.filt = RunPCA(sce.all.filt, npcs = 20)
#nExp <- round(ncol(sce.all.filt) * 0.04)
saveRDS(sce.all.filt, file="./00.data/01.single cell/scRNA1.rds")

