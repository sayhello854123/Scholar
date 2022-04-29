#引用包
library(plyr)
library(ggplot2)
library(ggpubr)

MSIplot <- function(scoreFile,cliFile,x,y,z,a){
#读取输入文件
trait="MSI" 
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(score) <- gsub("(.*?)\\_(.*?)\\.*", "\\2", rownames(score))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(cli), gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score)))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])
rt <- na.omit(rt)
rt$MSI=factor(rt$MSI, levels=c("MSS", "MSI-L", "MSI-H"))
rt$group=factor(rt$group, levels=c("Low", "High"))
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

#统计高低评分组病人数目
rt1=rt[,c(trait, "group")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
#计算高低评分组的百分率
df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
#百分比位置
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
#绘制百分率图
p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("riskScore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()
ggsave(filename = x,plot =p,
       width=4, height=5,dpi = 1000)
ggsave(filename = y,plot =p,
       width=4, height=5,dpi = 1000)



#设置比较组
rt2=rt[,c(trait, "riskScore")]
colnames(rt2)=c("trait", "riskScore")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制箱线图
boxplot=ggboxplot(rt2, x="trait", y="riskScore", fill="trait",
                  xlab="",
                  ylab="riskScore",
                  legend.title=trait,
                  palette=bioCol
)+ 
  stat_compare_means(comparisons=my_comparisons)

ggsave(filename = z,plot =boxplot,
       width=4, height=5,dpi = 1000)

ggsave(filename = a,plot =boxplot,
       width=4, height=5,dpi = 1000)
}

scoreFile="./00.data/READ/11.Model/m6Ascore_train_group.txt"    
cliFile="./00.data/READ/12.MSI/MSI.txt"                 
x <- "./05.step5/11.MSI/barplot.pdf"
y <- "./05.step5/11.MSI/barplot.jpg"
z <- "./05.step5/11.MSI/boxplot.pdf"
a <- "./05.step5/11.MSI/boxplot.jpg"
MSIplot(scoreFile,cliFile,x,y,z,a)
