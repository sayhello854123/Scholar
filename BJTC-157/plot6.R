library(jpeg)
library(ggpubr)
A0 <- readJPEG('./01.single cell/02.cellchat/heatmap1.jpg')
p0<-ggplot()+
  background_image(A0)+
  theme_void()
A5 <- readJPEG('./01.single cell/02.cellchat/heatmap2.jpg')
p5<-ggplot()+
  background_image(A5)+
  theme_void()
A6 <- readJPEG('./01.single cell/02.cellchat/heatmap.jpg')
p6<-ggplot()+
  background_image(A6)+
  theme_void()
A1 <- readJPEG('./01.single cell/02.cellchat/picture1.jpg')
p1<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('./01.single cell/02.cellchat/picture2.jpg')
p2<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('./01.single cell/02.cellchat/picture3.jpg')
p3<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('./01.single cell/02.cellchat/picture4.jpg')
p4<-ggplot()+
  background_image(A4)+
  theme_void()
P <- p1/p6/p2/p3/p4+plot_annotation(tag_levels = 'A')
ggsave("./01.single cell/02.cellchat/picture.png", plot =P , width =12, height = 30) 
