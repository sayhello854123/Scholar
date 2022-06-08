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
library(ggrepel)
library(msigdbr)
library(tidyverse)
library(patchwork)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
library(limma)

load('./00.data/01.single_data/scRNA_macrophage.Rdata')
##Ecotyper输入数据准备
# A <- scRNA@assays$RNA@counts
# a <- data.frame(A)
# metadata <- scRNA@meta.data
# write.table(a,file = './06.EcoTyper/macrophage.txt',sep = '\t',quote = F)
# write.table(metadata,file = './06.EcoTyper/macrophage_annotation.txt',sep = '\t',quote = F)
# data1 = read.delim("./06.EcoTyper/ecotyper-master/ecotyper-master/example_data/scRNA_CRC_data.txt", nrow = 5)
# data2 = read.delim("./06.EcoTyper/ecotyper-master/ecotyper-master/example_data/scRNA_CRC_annotation.txt", nrow = 5)

###提取生态型
metadata <- scRNA@meta.data
Eco <- read.table('./07.Ecological_diff/Ecotyper.txt',header = T,sep = '\t',row.names = 1,check.names = F)
same <- intersect(row.names(metadata),row.names(Eco))
metadata <- metadata[same,]
Eco <- Eco[same,]
meatdata1 <- cbind(metadata,Eco)
scRNA_macrophage <- subset(scRNA, cells=row.names(meatdata1))
rm(scRNA)
scRNA <- CreateSeuratObject(scRNA_macrophage@assays$RNA@counts, meta.data = meatdata1 )
rm(scRNA_macrophage)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
ElbowPlot(scRNA,ndims = 50)
pc.num = 1:30
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()
SaveH5Seurat(scRNA,'./00.data/01.single_data/seurat_EcoTyper.h5seurat',overwrite=T)
Convert('./00.data/01.single_data/seurat_EcoTyper.h5seurat','./00.data/01.single_data/seurat_EcoTyper.h5ad',overwrite=T)
system.time(save(scRNA, file = "./00.data/01.single_data/scRNA_macrophageM2.Rdata"))
##计算差异基因
Idents(scRNA) <- 'EcoTyper'
diff.mast <- FindAllMarkers(scRNA,test.use = 'MAST')
all.markers = diff.mast %>% subset(p_val_adj<0.05)
write.csv(all.markers,file = './07.Ecological_diff/alldiff.csv')
all_E1 <- all.markers[all.markers$cluster=='E1',]
all_E3 <- all.markers[all.markers$cluster=='E3',]
EcoTyper_gene <- read.table('./07.Ecological_diff/Ecotyper_geni_nfo.txt',header = T,sep = '\t',row.names = 1,check.names = F)
EcoTyper_E1 <- EcoTyper_gene[EcoTyper_gene$EcoTyper=='E1',]
EcoTyper_E3 <- EcoTyper_gene[EcoTyper_gene$EcoTyper=='E3',]
save(EcoTyper_E1,EcoTyper_E3,file = './07.Ecological_diff/eco.RData')
same1 <- intersect(row.names(all_E1),row.names(EcoTyper_E1))
same2 <- intersect(all_E3[,7],row.names(EcoTyper_E3))
same3 <- c(same1,same2)
alldiff <- read.csv('./07.Ecological_diff/alldiff.csv',header = T,row.names = 1)
EcoTyper_diff <- alldiff[same3,]
write.csv(EcoTyper_diff,file = './07.Ecological_diff/result_diff.csv')
if(T){
###火山图
allDiff <- EcoTyper_diff
allDiff$gene <- rownames(allDiff)
adjustP = 0.05
logFoldChange = 1
threshold <- as.factor(ifelse(allDiff$p_val_adj < adjustP &abs(allDiff$avg_log2FC) >logFoldChange,
                              ifelse(allDiff$avg_log2FC > 1 ,'Up','Down'),'Stable'))
p1 <- ggplot(
  #设置数据
  allDiff, 
  aes(x = avg_log2FC, 
      y=-log10(p_val_adj),
      colour=threshold)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#00FF00", "#d2dae2","#DF0029"))+
  
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(adjustP),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
allDiff$label = ifelse(allDiff$p_val_adj < 0.001 & abs(allDiff$avg_log2FC) >2.5, 
                       as.character(allDiff$gene),"")
P<- p1+geom_text_repel(data = allDiff, aes(x = avg_log2FC, 
                                          y = -log10(p_val_adj), 
                                          label = label),
                      size = 3,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE)
ggsave('./07.Ecological_diff/Ecotyper_diff.pdf',P,width = 8,height = 6,dpi = 1000)
ggsave('./07.Ecological_diff/Ecotyper_diff.jpg',P,width = 8,height = 6,dpi = 1000)
}

###相关性分析
same4 <- intersect(row.names(scRNA),row.names(EcoTyper_diff))
scRNA <- scRNA[same4,]
A <- scRNA@assays$RNA@counts
data1 <- data.frame(A)
save(data1,file = './07.Ecological_diff/cor_input.RData')
##批量画图
outTab=data.frame()
dir <- './07.Ecological_diff/cor/'
dir.create(dir)
for(i in 1:(nrow(data1)-1)){
  x=log2(as.numeric(data1[i,])+1)
  gene1 = rownames(data1)[i]
  for(j in (i+1):nrow(data1)){
    y=log2(as.numeric(data1[j,])+1)
    gene2 = rownames(data1)[j]
    corT=cor.test(x,y)
    z=lm(y~x)
    cor=corT$estimate
    cor=round(cor,3)
    pvalue=corT$p.value
    if (is.na(pvalue)) {
      next
    }
    print(gene1)
    print(gene2)
    if(pvalue<0.001){
      pval=signif(pvalue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=round(pvalue,3)}
    if((abs(cor)>0.8) & (pvalue<0.001)){
      picDir <- paste0(dir,gene1)
      dir.create(picDir)
      tiffFile=paste(gene1,"_",gene2,".cor.pdf",sep="")
      outTiff=paste(picDir,tiffFile,sep="/")
      pdf(file=outTiff,width =6,height = 6)
      plot(x,y, type="p",pch=16,col="blue",main=paste("Cor=",cor," (p-value=",pval,")",sep=""),
           cex=1, cex.lab=1, cex.main=1,cex.axis=1,
           xlab=paste(gene1,"expression"),
           ylab=paste(gene2,"expression") )
      lines(x,fitted(z),col=2)
      dev.off()
    }
    outTab=rbind(outTab,cbind(gene1,gene2,cor,pvalue))
  } 
}
write.csv(outTab,file = './07.Ecological_diff/cor/all_cor.csv')


###
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H')
gs=split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs = lapply(gs, unique)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
geneset <- gs
X= as.matrix(data1) 
gsva_es <- gsva(X, geneset, 
                mx.diff=FALSE, verbose=FALSE, 
                parallel.sz=8)
write.csv(gsva_es, "./07.Ecological_diff/gsva/gsva_output.csv", quote = F)

cellinfo <- scRNA@meta.data
group_list <- data.frame(rownames(cellinfo ), cellinfo$EcoTyper)
design <- model.matrix(~ 0 + factor(group_list$cellinfo.EcoTyper))
colnames(design) <- levels(factor(group_list$cellinfo.EcoTyper))
rownames(design) <- colnames(gsva_es)

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(E3-E1, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)
write.csv(x, "./07.Ecological_diff/gsva/gsva_limma.csv", quote = F)

#输出t值，用做的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, "./07.Ecological_diff/gsva/easy_input2_for39bar.csv", quote = F, row.names = F)

##作图
rm(list = ls())
options(stringsAsFactors = F) 
gc()
df <- read.csv("./07.Ecological_diff/gsva/easy_input2_for39bar.csv")
head(df)

#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)


p1 <- ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "inward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score\n  E3 Macrophage M2 versus E1 Macrophage M1")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave('./07.Ecological_diff/gsva/gsva.pdf',width = 7.5, height = 8,p1,dpi = 1000)
ggsave('./07.Ecological_diff/gsva/gsva.jpg',width = 7.5, height = 8,p1,dpi = 1000)
