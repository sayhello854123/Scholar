library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
if(T){
  rt=read.table('./00.data/READ/04.meger/merge.txt', header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp), colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  data=avereps(data)
  #ssGSEA分析
  cellMarker <- data.table::fread("./00.data/gene.csv",data.table = F)
  colnames(cellMarker)[2] <- "type"
  cellMarker <- lapply(split(cellMarker,cellMarker$type), function(x){
    dd = x$gene
    unique(dd)
  })
  ssgseaScore=gsva(data, cellMarker, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  ssgseaScore1=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(ssgseaScore1), ssgseaScore1)
  write.table(ssgseaOut,file="./00.data/READ/07.GSVA/ssGSEA.result.txt",sep="\t",quote=F,col.names=F)
  
  
  #读取分型文件
  cluster=read.table('./00.data/READ/06.cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
  row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
  #数据合并
  ssgseaScore1=t(ssgseaScore1)
  sameSample=intersect(row.names(ssgseaScore1), row.names(cluster))
  ssgseaScore1=ssgseaScore1[sameSample,,drop=F]
  cluster=cluster[sameSample,,drop=F]
  scoreCluster=cbind(ssgseaScore1, cluster)
  
  #把数据转换成ggplot2输入文件
  data=melt(scoreCluster, id.vars=c("cluster"))
  colnames(data)=c("m6Acluster", "metabolize", "Fraction")
  
  #绘制箱线图
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length(levels(factor(data[,"m6Acluster"])))]
  p=ggboxplot(data, x="metabolize", y="Fraction", color="m6Acluster", 
              ylab="metabolize infiltration",
              xlab="",
              legend.title="m6Acluster",
              palette=bioCol)
  p=p+rotate_x_text(50)
  
  plot2 <- p+stat_compare_means(aes(group=m6Acluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                                       symbols = c("***", "**", "*", "ns")),label = "p.signif")
  ggsave(filename = "./05.step5/05.ssGSEA/metabolize_boxplot.pdf",plot =plot2,
         width=8, height=6.5,dpi = 1000)
  ggsave(filename = "./05.step5/05.ssGSEA/metabolize_boxplot.jpg",plot =plot2,
         width=8, height=6.5,dpi = 1000)
  
}