library(limma)
library(ggplot2)
library(VennDiagram)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
if(T){
  #读取输入文件,并对输入文件进行整理
  rt=read.table('00.data/ESCA/04.meger/m6A_merge.txt', header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  data=t(data)
  cluster=read.table('./00.data/ESCA/06.cluster/m6aCluster1.txt', header=T, sep="\t", check.names=F, row.names=1)
  row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
  sameSample=intersect(row.names(data), row.names(cluster))
  data=data[sameSample,,drop=F]
  cluster=cluster[sameSample,,drop=F]
  #PCA分析
  data.pca=prcomp(data, scale. = TRUE)
  pcaPredict=predict(data.pca)
  write.table(pcaPredict, file="./00.data/ESCA/08.PCA/newTab.xls", quote=F, sep="\t")
  
  
  m6Acluster=as.vector(cluster[,1])
  
  #设置颜色
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  
  m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]
  
  #可视化
  PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], m6Acluster=m6Acluster)
  PCA.mean=aggregate(PCA[,1:2], list(m6Acluster=PCA$m6Acluster), mean)
  pdf(file="./03.step3/06.PCA/PCA.pdf", height=5, width=6.5)
  ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = m6Acluster)) +
    scale_colour_manual(name="m6Acluster", values =m6aCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$m6Acluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  dev.off()
  p1 <- ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = m6Acluster)) +
    scale_colour_manual(name="m6Acluster", values =m6aCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$m6Acluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggsave(filename = "./03.step3/06.PCA/PCA.jpg",plot =p1,
         height=5, width=6.5,dpi = 1000)
}
