library(limma)
library(ggplot2)
library(VennDiagram)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
if(T){
  adj.P.Val.Filter=0.05
  
  #读取输入文件，并对输入文件整理
  rt=read.table('./00.data/ESCA/04.meger/merge.txt', header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  
  #读取cluster文件
  cluster=read.table('./00.data/ESCA/06.cluster/m6aCluster1.txt', header=T, sep="\t", check.names=F, row.names=1)
  row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
  #提取交集文件
  sameSample=intersect(colnames(data), row.names(cluster))
  data=data[,sameSample]
  cluster=cluster[sameSample,]
  
  #差异分析
  logFCfilter=0
  geneList=list()
  Type=as.vector(cluster)
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  comp=combn(levels(factor(Type)), 2)
  allDiffGenes=c()
  dir.create('./00.data/ESCA/09.PAAD_VENN/')
  dir <- './00.data/ESCA/09.PAAD_VENN/'
  for(i in 1:ncol(comp)){
    fit=lmFit(data, design)
    contrast=paste0(comp[2,i], "-", comp[1,i])
    #print(contrast)
    cont.matrix=makeContrasts(contrast, levels=design)
    fit2=contrasts.fit(fit, cont.matrix)
    fit2=eBayes(fit2)
    
    #输出所有基因差异情况
    allDiff=topTable(fit2,adjust='fdr',number=200000)
    allDiffOut=rbind(id=colnames(allDiff),allDiff)
    write.table(allDiffOut, file=paste0(dir,contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
    
    #输出差异结果
    diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
    diffSigOut=rbind(id=colnames(diffSig),diffSig)
    write.table(diffSigOut, file=paste0(dir,contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
    geneList[[contrast]]=row.names(diffSig)
  }
  
  #绘制venn图
  venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
  pdf(file="./03.step3/07.VENN/venn.pdf", width=5, height=5)
  grid.draw(venn.plot)
  dev.off()
  
  #保存交集基因
  interGenes=Reduce(intersect,geneList)
  write.table(file="./00.data/ESCA/09.PAAD_VENN/interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)
  
  #保存交集基因的表达量
  interGeneExp=data[interGenes,]
  interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
  write.table(interGeneExp, file="./00.data/ESCA/09.PAAD_VENN/interGeneExp.txt", sep="\t", quote=F, col.names=F)
}
