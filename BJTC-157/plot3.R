library(jpeg)
library(ggpubr)
library(patchwork)
A1 <- readJPEG('./01.single cell/celltype.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./01.single cell/immune_group.jpg')
p2<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./01.single cell/immune_cell.jpg')
p3<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./01.single cell/immune_heatmap.jpg')
p4<-ggplot()+
  background_image(A4)+
  theme_void()

P <- (p1|p2)/(p3|p4)+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/01.Seurat/picture1.png", plot =P , width =20, height = 16) 


####bulk###
A1 <- readJPEG('./02.bulk/plot1/1.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./02.bulk/plot1/2.jpg')
p2<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./02.bulk/plot1/3.jpg')
p3<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./02.bulk/plot1/4.jpg')
p4<-ggplot()+
  background_image(A4)+
  theme_void()





A1 <- readJPEG('./02.bulk/plot1/1.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./02.bulk/plot1/2.jpg')
p2<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./02.bulk/plot1/3.jpg')
p3<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./02.bulk/plot1/4.jpg')
p4<-ggplot()+
  background_image(A4)+
  theme_void()
A5 <- readJPEG('./02.bulk/plot1/5.jpg')
p5<-ggplot()+
  background_image(A5)+
  theme_void()
A6 <- readJPEG('./02.bulk/plot1/6.jpg')
p6<-ggplot()+
  background_image(A6)+
  theme_void()
A7 <- readJPEG('./02.bulk/plot1/7.jpg')
p7<-ggplot()+
  background_image(A7)+
  theme_void()
A8 <- readJPEG('./02.bulk/plot1/8.jpg')
p8<-ggplot()+
  background_image(A8)+
  theme_void()

P <- (p1|p2)/(p3|p4|p5)/(p6|p7|p8)+plot_annotation(tag_levels = 'A')
P<- P* theme_minimal()
ggsave("./02.bulk/plot1/picture1.png", plot =P , width =14, height = 14) 



A1 <- readJPEG('./02.bulk/plot2/1.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./02.bulk/plot2/2.jpg')
p2<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./02.bulk/plot2/3.jpg')
p3<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./02.bulk/plot2/4.jpg')
p4<-ggplot()+
  background_image(A4)+
  theme_void()
A5 <- readJPEG('./02.bulk/plot2/5.jpg')
p5<-ggplot()+
  background_image(A5)+
  theme_void()
A6 <- readJPEG('./02.bulk/plot2/6.jpg')
p6<-ggplot()+
  background_image(A6)+
  theme_void()
A7 <- readJPEG('./02.bulk/plot2/7.jpg')
p7<-ggplot()+
  background_image(A7)+
  theme_void()
A8 <- readJPEG('./02.bulk/plot2/8.jpg')
p8<-ggplot()+
  background_image(A8)+
  theme_void()
A9 <- readJPEG('./02.bulk/plot2/9.jpg')
p9<-ggplot()+
  background_image(A9)+
  theme_void()
P <- (p1|p2|p3)/(p4|p5)/(p6|p7)/(p8|p9)+plot_annotation(tag_levels = 'A')
P<- P* theme_minimal()
ggsave("./02.bulk/plot2/picture1.png", plot =P , width =14, height = 14) 
