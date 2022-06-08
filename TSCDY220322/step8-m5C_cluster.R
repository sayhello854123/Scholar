library(ggpubr)

rt <- read.table('./03.survial/merge.txt',header = T,sep ='\t',row.names = 1,check.names = F)
rt1 <- t(rt)
cluster=read.table("./04.cluster/Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
same <- intersect(row.names(rt1),row.names(cluster))

data <- cbind(rt1[same, ,drop=F],cluster[same, ,drop=F])
data1 <- data[,c('ERCC2','cluster')]
data1 <- read.csv('./2588.csv',header = T,row.names = 1)
x=colnames(data1)[2]
y=colnames(data1)[1]

#设置比较租
group=levels(factor(data1$cluster))
data1$group=factor(data1$cluster, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}



var="ERCC2"
p <- wilcox.test(data1[which(data1$cluster == "A"),var],
                 data1[which(data1$cluster == "B"),var])$p.value

jco <- c("#0066FF","#FF0066")
p15 <- ggplot(data = data1,aes(x = cluster, y = ERCC2, fill = cluster))+
  scale_fill_manual(values = jco[1:2]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("XPD gene expression")) +
  xlab("m5C cluster")  +
  annotate(geom="text", cex=6,
           x=1.5, y=8, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3)))), # 添加P值
           color="black") + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))
ggsave('./04.cluster/violin.pdf',p15,width = 8,height = 6,dpi = 1000)
ggsave('./04.cluster/violin.jpg',p15,width = 8,height = 6,dpi = 1000)

library(corrplot)
library(survival)
library(survminer)
rt <- read.table('./04.cluster/COR.txt',header = T,row.names = 1,check.names = F,sep = '\t')
M=cor(rt)
res1=cor.mtest(data, conf.level = 0.95)

size=c()
for(i in 1:ncol(M)){
  size=c(size,M[1:i,i])
}
size=as.numeric(size)
col1 <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow",
                           "#FF7F00","red","#7F0000"))
pdf(file="./04.cluster/corpot.pdf",width=8,height=8)
 corrplot(M, method = "circle",
         type = "upper",
         number.cex = abs(size),
         col=col1(200),
         tl.col="black",
         addCoef.col="black"
)
dev.off()
