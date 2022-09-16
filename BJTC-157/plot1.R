rt <- read.table('1.txt',sep = '\t',header = T,check.names = F,row.names = 1)
type <- read.table('2.txt',sep = '\t',header = T,check.names = F,row.names = 1)
rownames(rt)[6] <- 'Naive')
table(type$Patient)
group <- type$Patient
H18 <- rt[,group=="H18"]
H21 <- rt[,group=="H21"]
H28 <- rt[,group=="H28"]
H30 <- rt[,group=="H30"]
H37 <- rt[,group=="H37"]
H38<- rt[,group=="H38"]
rt <- cbind(H18,H21,H28,H30,H37,H38)

patient <- c("H18" = "red","H21" = "green","H28" = "blue",
             "H30" = "hotpink","H37" = "mediumorchid1","H38" = "yellow")
column_ha = HeatmapAnnotation(Type=type$Type,Patient =type$Patien,
                              col = list(Patient =patient),
                              Type=c('Type'='mediumblue'))
pdf(file='./01.single cell/heat1_T.pdf',width = 6, height = 6)
Heatmap(rt,show_column_names = F,top_annotation = column_ha,
        cluster_rows = F,color = colorRampPalette(c("blue", "white", "red"))(50),
        row_names_gp = gpar( fontsize=12),
        heatmap_legend_param = list(
          title= '', title_position = "topcenter", 
          legend_height=unit(8,"cm"), legend_direction="vertical")
        )
dev.off()



library("jpeg")
library(ggpubr)
library(ggplot2)
A1 <- readJPEG('./01.single cell/01.Seurat/tSNE.jpg')
p0<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./01.single cell/01.Seurat/cor1.jpg')
p1<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./01.single cell/01.Seurat/cor2.jpg')
p2<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./01.single cell/01.Seurat/cor3.jpg')
p3<-ggplot()+
  background_image(A4)+
  theme_void()
P <- p0/(p1|p2|p3)+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/01.Seurat/picture1.png", plot =P , width =20, height = 16) 

A1 <- readJPEG('./01.single cell/01.Seurat/celltype.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./01.single cell/01.Seurat/immune_group.jpg')
p1<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./01.single cell/01.Seurat/immune_cell.jpg')
p2<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./01.single cell/01.Seurat/immune_heatmap.jpg')
p3<-ggplot()+
  background_image(A4)+
  theme_void()

P <- (p1|p2)/(p3|p4)+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/01.Seurat/picture1.png", plot =P , width =20, height = 16) 


