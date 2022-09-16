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
ggsave("./01.single cell/01.Seurat/picture1.png", plot =P , width =20, height = 16,dpi = 1000) 
