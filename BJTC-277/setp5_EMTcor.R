rm(list = ls())
options(stringsAsFactors = F) 
gc()
library(GSVA)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)
library(grid)
library(ggplot2)
library(patchwork)

##Reading marker genes
cellMarker <- data.table::fread('./04.EMT_cor/angiogenesis.txt',data.table = F)
colnames(cellMarker)[2] <- "celltype"
cellMarker <- lapply(split(cellMarker,cellMarker$celltype), function(x){
  dd = x$gene
  unique(dd)
})
##
load('./00.data/01.single_data/scRNA_qc.Rdata')
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
##如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
a <- as.data.frame(scRNA@assays$RNA@scale.data) 
expr <- as.matrix(a )

### 3.使用ssGSEA量化免疫浸润
gene1="ANGIOGENESIS"           
gene2="EMT"
gene3 <- "STEMNES"
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
write.csv(gsva_data,file = './04.EMT_cor/gsva_result.csv')

x=as.numeric(gsva_data[gene1,])
y=as.numeric(gsva_data[gene3,])
#相关性分析
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
  xlab('ANGIOGENESIS(hallmark)')+ylab("STEMNES(hallmark)")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p1=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
ggsave("./04.EMT_cor/cor1.pdf", plot = p1, width = 8, height = 8)
ggsave("./04.EMT_cor/cor1.jpg", plot = p1, width = 8, height = 8,dpi = 600)

x=as.numeric(gsva_data[gene1,])
y=as.numeric(gsva_data[gene2,])

#相关性分析
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p2=ggplot(df1, aes(x, y)) + 
  xlab('ANGIOGENESIS(hallmark)')+ylab("EMT(hallmark)")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p2, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
ggsave("./04.EMT_cor/cor2.pdf", plot = p2, width = 8, height = 8)
ggsave("./04.EMT_cor/cor2.jpg", plot = p2, width = 8, height = 8,dpi = 600)

x=as.numeric(gsva_data[gene3,])
y=as.numeric(gsva_data[gene2,])

#相关性分析
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p3=ggplot(df1, aes(x, y)) + 
  xlab("STEMNES(hallmark)")+ylab("EMT(hallmark)")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p3=ggMarginal(p3, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
ggsave("./04.EMT_cor/cor3.pdf", plot = p3, width = 8, height = 8)
ggsave("./04.EMT_cor/cor3.jpg", plot = p3, width = 8, height = 8,dpi = 600)

##组合图片
library("jpeg")
library(ggpubr)
library(ggplot2)
library(cowplot)
pir <- './04.EMT_cor/'
files=dir('./04.EMT_cor/')                        #获取目录下所有文件
files=grep("jpg$",files,value=T)   
files <- paste0(pir,files)
plotp <- list()
for(i in files){
  A1 <- readJPEG(i)
  gene1 <- gsub("(.*?)\\/(.*?)\\/(.*?)", "\\3", i)
  gene <- strsplit(gene1, ".",fixed = TRUE)[[1]][1]
  p<-ggplot()+
    background_image(A1)+
    theme_void()
  plotp[[gene]] <- p
}
p4 <- plot_grid(plotlist=plotp, ncol=3,labels = "AUTO")
ggsave(filename = '04.EMT_cor/all.pdf',p4,width = 24,height = 8,dpi = 1000)
ggsave(filename = '04.EMT_cor/all.jpg',p4,width = 24,height = 8,dpi = 1000)

p5 <- plot_grid(plotlist=plotp, ncol=1,labels = "AUTO")
ggsave(filename = '04.EMT_cor/all1.pdf',p5,width = 8,height = 24,dpi = 1000)
ggsave(filename = '04.EMT_cor/all1.jpg',p5,width = 8,height = 24,dpi = 1000)
