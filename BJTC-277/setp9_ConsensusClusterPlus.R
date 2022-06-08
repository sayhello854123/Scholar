library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(pheatmap)
library(scatterplot3d)
library(limma)
library(VennDiagram)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

load('./08.WGCNA/wgcna_input.RData')
dim(data)
data1 <- t(data)
group=sapply(strsplit(rownames(data1),"\\-"), "[", 4)
table(group)
group=sapply(strsplit(group,""), "[", 1)
table(group)
rt1=data1[group==0,]
dim(rt1)
rt2 <- t(rt1)
gene <- read.table('./08.WGCNA/venn/intersect.txt',sep = '\t')
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
write.table(clusterOut, file="./09.ConsensusClusterPlus/Cluster1.txt", sep="\t", quote=F, col.names=F)

cli <-data.table::fread('./08.WGCNA/TCGA-STAD.survival.tsv',data.table = F)
row.names(cli) <- cli[,1]
cli$futime <- cli$OS.time/365
cli$fustat <- cli$OS
cli <- cli[,-c(1:4)]
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
rt <- read.csv('./09')
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
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
ggsurvplot(fit, 
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
