library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(pheatmap)
library(limma)
library(VennDiagram)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

load('./09.ConsensusClusterPlus/wgcna_input.RData')
dim(data)
data1 <- t(data)
group=sapply(strsplit(rownames(data1),"\\-"), "[", 4)
table(group)
group=sapply(strsplit(group,""), "[", 1)
table(group)
rt1=data1[group==0,]
dim(rt1)
rt2 <- t(rt1)
gene <- read.table('./09.ConsensusClusterPlus/intersect.txt',sep = '\t')
same <- intersect(row.names(rt2),gene[,1])
data <- rt2[same,]
data=as.matrix(data)

#聚类
maxK=8
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=1000,
                             pItem=0.8,
                             pFeature=1,
                             title='./09.ConsensusClusterPlus/',
                             clusterAlg="pam",
                             distance="pearson",
                             innerLinkage="complete",
                             seed=123456,
                             plot="png")
#输出分型结果
clusterNum=2       #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="./09.ConsensusClusterPlus/Cluster.txt", sep="\t", quote=F, col.names=F)


cli <-data.table::fread('./09.ConsensusClusterPlus/TCGA-STAD.survival.tsv',data.table = F)
row.names(cli) <- cli[,1]
cli$futime <- cli$OS.time/365
cli$fustat <- cli$OS
cli <- cli[,-c(1:4)]
cluster <- read.table('./09.ConsensusClusterPlus/Cluster.txt',header = T,row.names = 1,sep = '\t',check.names = F)
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#生存差异统计
length=length(levels(factor(rt$cluster)))
diff=survdiff(Surv(futime, fustat) ~ cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)

bioCol=c("#0066FF","#FF0066","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="cluster",
                   legend.labs=levels(factor(rt[,"cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.28)
pdf(file='./09.ConsensusClusterPlus/survival.pdf',onefile = FALSE,width=6.5,height=5)
print(surPlot)
dev.off()

###热图

exp=t(data)
cluster <- read.table('./09.ConsensusClusterPlus/Cluster1.txt',header = T,row.names = 1,sep = '\t',check.names = F)
#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)

data <- expCluster
data=data[order(data$cluster),]
Type=as.data.frame(data[,((ncol(exp)+1):ncol(data))])
colnames(Type)[1] <- c('cluster')
rownames(Type) <- rownames(data)
data=t(data[,1:ncol(exp)])
#聚类颜色
bioCol=c("#0066FF","#FF0066","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(Type$cluster)))]
names(CluCol)=levels(factor(Type$cluster))
ann_colors[["cluster"]]=CluCol
pdf("./09.ConsensusClusterPlus/heatmap.pdf", height=5, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames =F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

##PCA
exp=t(data)
cluster <- read.table('./09.ConsensusClusterPlus/Cluster1.txt',header = T,row.names = 1,sep = '\t',check.names = F)
same <- intersect(row.names(cluster),row.names(exp))
#PCA分析
data.pca=prcomp(exp, scale. = TRUE)
pcaPredict=predict(data.pca)
cluster1=as.vector(cluster[,1])
m6aCluCol=bioCol[1:length(levels(factor(cluster1)))]
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2],cluster=cluster1)
PCA.mean=aggregate(PCA[,1:2], list(cluster=PCA$cluster), mean)
p1 <- ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = cluster)) +
  scale_colour_manual(name="cluster", values =CluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave('./09.ConsensusClusterPlus/PCA.jpg',p1,width = 7,height = 6,dpi = 1000)
ggsave('./09.ConsensusClusterPlus/PCA.pdf',p1,width = 7,height = 6,dpi = 1000)


##TAM评分
