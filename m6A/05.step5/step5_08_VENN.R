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
  rt=read.table('./00.data/READ/04.meger/merge.txt', header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  
  #读取cluster文件
  cluster=read.table('./00.data/READ/06.cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
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
  dir.create('./00.data/READ/09.VENN/')
  dir <- './00.data/READ/09.VENN/'
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
  color=c( "#3C5488B2","#00A087B2", 
           "#F39B7FB2","#91D1C2B2", 
           "#8491B4B2", "#DC0000B2", 
           "#7E6148B2","yellow", 
           "darkolivegreen1", "lightskyblue", 
           "darkgreen", "deeppink", "khaki2", 
           "firebrick", "brown1", "darkorange1", 
           "cyan1", "royalblue4", "darksalmon", 
           "darkgoldenrod1", "darkseagreen", "darkorchid")
  pdf(file="./05.step5/07.venn/venn.pdf", width=5, height=5)
  venn(geneList, zcolor = color[1:(length(geneList))])
  dev.off()
  
  #保存交集基因
  interGenes=Reduce(intersect,geneList)
  write.table(file="./00.data/READ/09.VENN/interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)
  
  #保存交集基因的表达量
  interGeneExp=data[interGenes,]
  interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
  write.table(interGeneExp, file="./00.data/READ/09.VENN/interGeneExp.txt", sep="\t", quote=F, col.names=F)
}
