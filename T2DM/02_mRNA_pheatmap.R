rm(list = ls())  
options(stringsAsFactors = F)#清空环境
load('./01.data/02pheatmap_input.RData')##加载热图数据
dir.create('02.mRNA_pheatmap')

library(pheatmap)
Type=c(rep("non−diabetic condition",10),rep("diabetic condition_donor",10))#设置分型
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file='./02.mRNA_pheatmap/mRNA_pheatmap.pdf',width=8,height=10)
pheatmap(rt,
         annotation=Type,
         cluster_cols = T,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_colnames = T,
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

