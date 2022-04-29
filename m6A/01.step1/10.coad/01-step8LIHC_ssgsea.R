library(ConsensusClusterPlus) 
library(pheatmap)
library(limma)
library(GSEABase)
library(GSVA)
if(T){
  rt=read.table('./00.data/LIHC/megerData/merge.txt', header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp), colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  data=avereps(data)
  
  geneSets=getGmt('./00.data/c2.cp.kegg.v7.2.symbols.gmt', geneIdType=SymbolIdentifier())
  gsvaResult=gsva(data, 
                  geneSets, 
                  min.sz=10, 
                  max.sz=500, 
                  verbose=TRUE,
                  parallel.sz=1)
  gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
  write.table(gsvaOut, file="./00.data/LIHC/cluster/gsvaOut.txt", sep="\t", quote=F, col.names=F)
  
  cluster=read.table('./00.data/LIHC/cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
  row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
  #数据合并
  gsvaResult=t(gsvaResult)
  sameSample=intersect(row.names(gsvaResult), row.names(cluster))
  gsvaResult=gsvaResult[sameSample,,drop=F]
  cluster=cluster[sameSample,,drop=F]
  gsvaCluster=cbind(gsvaResult, cluster)
  Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
  gsvaCluster=cbind(gsvaCluster, Project)
  
  #差异分析
  dir <- './00.data/LIHC/cluster/'
  a <- './01.step1/06.m6ALIHC_cluster/'
  adj.P.Val.Filter=0.05
  allType=as.vector(gsvaCluster$cluster)
  comp=combn(levels(factor(allType)), 2)
  for(i in 1:ncol(comp)){
    #样品分组
    treat=gsvaCluster[gsvaCluster$cluster==comp[2,i],]
    con=gsvaCluster[gsvaCluster$cluster==comp[1,i],]
    data=rbind(con, treat)
    #差异分析
    Type=as.vector(data$cluster)
    ann=data[,c(ncol(data), (ncol(data)-1))]
    data=t(data[,-c((ncol(data)-1), ncol(data))])
    design=model.matrix(~0+factor(Type))
    colnames(design)=levels(factor(Type))
    fit=lmFit(data, design)
    contrast=paste0(comp[2,i], "-", comp[1,i])
    cont.matrix=makeContrasts(contrast, levels=design)
    fit2=contrasts.fit(fit, cont.matrix)
    fit2=eBayes(fit2)
    
    #输出所有通路的差异情况
    allDiff=topTable(fit2,adjust='fdr',number=200000)
    allDiffOut=rbind(id=colnames(allDiff),allDiff)
    write.table(allDiffOut, file=paste0(dir,contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
    
    #输出显著的差异
    diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
    diffSigOut=rbind(id=colnames(diffSig),diffSig)
    write.table(diffSigOut, file=paste0(dir,contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
    
    #聚类颜色
    bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    ann_colors=list()
    m6aCluCol=bioCol[1:length(levels(factor(allType)))]
    names(m6aCluCol)=levels(factor(allType))
    ann_colors[["cluster"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
    
    #绘制差异通路热图
    termNum=20
    diffTermName=as.vector(rownames(diffSig))
    diffLength=length(diffTermName)
    if(diffLength<termNum){termNum=diffLength}
    hmGene=diffTermName[1:termNum]
    hmExp=data[hmGene,]
    pdf(file=paste0(a,contrast,".heatmap.pdf"),height=6,width=10)
    pheatmap(hmExp, 
             annotation=ann,
             annotation_colors = ann_colors,
             color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
             cluster_cols =F,
             show_colnames = F,
             gaps_col=as.vector(cumsum(table(Type))),
             scale="row",
             fontsize = 10,
             fontsize_row=7,
             fontsize_col=10)
    dev.off()
  }
}
if(T){
cellMarker <- data.table::fread("./00.data/gene.csv",data.table = F)
colnames(cellMarker)[2] <- "type"
cellMarker <- lapply(split(cellMarker,cellMarker$type), function(x){
  dd = x$gene
  unique(dd)
})
rt=read.table('./00.data/LIHC/megerData/merge.txt', header=T, sep="\t", check.names=F)
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
rt=read.table('./00.data/LIHC/megerData/merge.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
ssgseaScore=gsva(data, cellMarker, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="./00.data/LIHC/cluster/ssGSEA.result.txt",sep="\t",quote=F,col.names=F)


#读取分型文件
cluster=read.table('./00.data/LIHC/cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
#数据合并
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

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
pdf(file="./01.step1/06.m6ALIHC_cluster/metabolize_boxplot.pdf", width=8, height=6.5)                          #输出图片文件
p+stat_compare_means(aes(group=m6Acluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()
plot2 <- p+stat_compare_means(aes(group=m6Acluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
ggsave(filename = "./01.step1/06.m6ALIHC_cluster/metabolize_boxplot.jpg",plot =plot2,
       width=8, height=6.5,dpi = 1000)
}
