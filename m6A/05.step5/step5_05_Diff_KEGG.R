library(ConsensusClusterPlus) 
library(pheatmap)
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
  
  geneSets=getGmt('./00.data/c2.cp.kegg.v7.2.symbols.gmt', geneIdType=SymbolIdentifier())
  gsvaResult=gsva(data, 
                  geneSets, 
                  min.sz=10, 
                  max.sz=500, 
                  verbose=TRUE,
                  parallel.sz=1)
  gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
  write.table(gsvaOut, file="./00.data/READ/07.GSVA/gsvaOut.txt", sep="\t", quote=F, col.names=F)
  
  cluster=read.table('./00.data/READ/06.cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
  row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
  #数据合并
  gsvaResult <- read.table('./00.data/READ/07.GSVA/gsvaOut.txt', header=T, sep="\t", check.names=F, row.names=1)
  gsvaResult=t(gsvaResult)
  sameSample=intersect(row.names(gsvaResult), row.names(cluster))
  gsvaResult=gsvaResult[sameSample,,drop=F]
  cluster=cluster[sameSample,,drop=F]
  gsvaCluster=cbind(gsvaResult, cluster)
  Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
  gsvaCluster=cbind(gsvaCluster, Project)
  
  #差异分析
  dir <- './00.data/READ/07.GSVA/'
  a <- './05.step5/04.kegg/'
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