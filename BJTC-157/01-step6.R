library(GSVA)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)
library(ggimage)
library(grid)
library(ggplot2)
library(patchwork)
library(EBImage)
library(imager)
library("jpeg")
library(ggpubr)
cellMarker <- data.table::fread('./00.data/02.gene get/angiogenesis.txt',data.table = F)
colnames(cellMarker)[2] <- "celltype"
cellMarker <- lapply(split(cellMarker,cellMarker$celltype), function(x){
  dd = x$gene
  unique(dd)
})
scRNA1 <- readRDS("./00.data/01.single cell/scRNA1.rds")
a <- as.data.frame(scRNA1@assays$RNA@scale.data) 

expr <- as.matrix(a )

### 3.使用ssGSEA量化免疫浸润
gene1="ANGIOGENESIS"           
gene2="EMT"
gene3 <- "STEMNES"
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")

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
ggsave("./01.single cell/cor1.pdf", plot = p1, width = 5, height = 5)

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
ggsave("./01.single cell/cor2.pdf", plot = p2, width = 5, height = 5)

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
ggsave("./01.single cell/cor3.pdf", plot = p3, width = 5, height = 5) 

A1 <- readJPEG('./01.single cell/cor1.jpg')
p0<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./01.single cell/cor2.jpg')
p1<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./01.single cell/cor3.jpg')
p2<-ggplot()+
  background_image(A3)+
  theme_void()
P <- (p0|p1|p2)+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/cor.png", plot =P , width = 15, height = 5) 
