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
if(T){
data=read.table('02.bulk/cluster/cluster.txt', header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
table(group)
#聚类
maxK=8
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=1000,
                             pItem=0.8,
                             pFeature=1,
                             title='02.bulk/cluster/',
                             clusterAlg="pam",
                             distance="pearson",
                             innerLinkage="complete",
                             seed=123456,
                             plot="png")
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK = Kvec[which.min(PAC)] 
icl = calcICL(results,
              title='02.bulk/cluster/',
              plot="pdf")

#输出分型结果
clusterNum=4       #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="./02.bulk/cluster/Cluster1.txt", sep="\t", quote=F, col.names=F)
}
if(T){
##合并生存数据
cli <-data.table::fread('./00.data/03.LIHC/TCGA-LIHC.survival.tsv',data.table = F)
row.names(cli) <- cli[,1]
cli$futime <- cli$OS.time/365
cli$fustat <- cli$OS
cli <- cli[,-c(1:4)]
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

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
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
pdf(file='./02.bulk/cluster/survival.pdf',onefile = FALSE,width=7,height=6.5)
print(surPlot)
dev.off()
}
###热图
if(T){
###热图
data=read.table('02.bulk/cluster/cluster.txt', header=T, sep="\t", check.names=F, row.names=1)
exp=t(data)
cluster <-read.table('02.bulk/cluster/Cluster1.txt', header=T, sep="\t", check.names=F, row.names=1) 
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
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(Type$cluster)))]
names(CluCol)=levels(factor(Type$cluster))
ann_colors[["cluster"]]=CluCol
pdf("./02.bulk/cluster/heatmap.pdf", height=5, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()
}
####M2相关性
if(T){
####M2相关性
cluster <-read.table('02.bulk/cluster/Cluster1.txt', header=T, sep="\t", check.names=F, row.names=1) 
data=read.table('02.bulk/cluster/cluster.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
same <- intersect(row.names(cluster),row.names(data))
rt <- cbind(cluster[same,],data[same,])
colnames(rt)[1] <- c('cluster')
rt <- data.frame(rt)
data=rt[,c(2:ncol(rt))]
Type=rt[,1]  #提取分组信息
var=colnames(rt)[1]
#颜色
group=levels(factor(Type))
bioCol=c("red","blue","green","yellow")
col= bioCol[match(Type,group)]
#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
pdf(file='02.bulk/cluster/pca.pdf', height=5, width=6)
par(oma=c(0.5,0.5,0.5,0.5))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=col)
legend("top", legend =group,pch = 16, inset = -0.2, box.col="white",xpd = TRUE, horiz = TRUE,col=bioCol[1:length(group)])
dev.off()
}
###富集文件构建###
if(T){
###富集文件构建###
data1 <- read.table('./02.bulk/survival/data.txt',sep = '\t',header = T,check.names = F)
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample,]
data <- cbind(gene_id = rownames(data), description = "na", data)
## 创建gct文件
gct_file = paste0("my.gct")
sink(gct_file)
cat("#1.2\n")
cat(paste0(nrow(data), "\t", (length(colnames(data)) - 2), "\n"))
sink()
write.table(data, gct_file, append = T, quote = F, row.names = F, col.names = T,
            sep = "\t")
cluster <-read.table('02.bulk/cluster/Cluster1.txt', header=T, sep="\t", check.names=F, row.names=1)
group_list <- cluster[,1]
cls_file = paste0("my.cls")
sink(cls_file)
cat(paste0(length(group_list), "\t", length(unique(group_list)), "\t1\n"))
cat(paste0("#", paste(unique(group_list), collapse = "\t"), "\n"))
cat(paste(group_list, collapse = "\t"), "\n")
sink()
}
###富集图##
if(T){
gseaplot <- function(dir,files,x,outfile){
  files <- paste0(dir,files)
  data=lapply(files, read.delim)                             #读取每个文件
  names(data)=sapply(strsplit(files,'/'),'[',5)
  dataSet=ldply(data, data.frame)
  dataSet$pathway = gsub(".tsv", "", dataSet$.id)
  #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #处理后有73种差异还比较明显的颜色，基本够用
  #gseaCol = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
  
  gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
    geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
    labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
    geom_hline(yintercept = 0) + 
    guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
  pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
    scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
    labs(x = x, y = "", title = "") + 
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)
  gGsea = ggplot_gtable(ggplot_build(pGsea))
  gGene = ggplot_gtable(ggplot_build(pGene))
  maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
  gGsea$widths = as.list(maxWidth)
  gGene$widths = as.list(maxWidth)
  dev.off()
  pdf(file=outfile,    #输出图片的文件
      width=9,                    #设置输出图片高度
      height=5.5)                 #设置输出图片高度
  par(mar=c(5,5,2,5))
  grid.arrange(arrangeGrob(gGsea, gGene, nrow=2, heights=c(.8,.25)))
  dev.off()
}
##A##
dir <- './02.bulk/cluster/A/'
files=grep(".tsv", dir('./02.bulk/cluster/A'), value=T) #获取目录下的所有tsv文件
x <- "A cluster score<----------->Rest cluster score"
outfile <- "./02.bulk/cluster/A_multipleGSEA.pdf"
gseaplot(dir=dir,files=files,x=x,outfile=outfile)

##B##
dir <- './02.bulk/cluster/B/'
files=grep(".tsv", dir('./02.bulk/cluster/B'), value=T) #获取目录下的所有tsv文件
x <- "B cluster score<----------->Rest cluster score"
outfile <- "./02.bulk/cluster/B_multipleGSEA.pdf"
gseaplot(dir=dir,files=files,x=x,outfile=outfile)

##C##
dir <- './02.bulk/cluster/C/'
files=grep(".tsv", dir('./02.bulk/cluster/C'), value=T) #获取目录下的所有tsv文件
x <- "C cluster score<----------->Rest cluster score"
outfile <- "./02.bulk/cluster/C_multipleGSEA.pdf"
gseaplot(dir=dir,files=files,x=x,outfile=outfile)


##D##
dir <- './02.bulk/cluster/D/'
files=grep(".tsv", dir('./02.bulk/cluster/D'), value=T) #获取目录下的所有tsv文件
x <- "D cluster score<----------->Rest cluster score"
outfile <- "./02.bulk/cluster/D_multipleGSEA.pdf"
gseaplot(dir=dir,files=files,x=x,outfile=outfile)
}
