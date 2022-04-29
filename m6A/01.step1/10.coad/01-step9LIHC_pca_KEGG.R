library(limma)
library(ggplot2)
library(VennDiagram)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
if(T){
#读取输入文件,并对输入文件进行整理
rt=read.table('00.data/LIHC/megerData/m6Amerge.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)
#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="./00.data/LIHC/PCA/newTab.xls", quote=F, sep="\t")

cluster=read.table('./00.data/LIHC/cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
m6Acluster=as.vector(cluster[,1])

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")

m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]

#可视化
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], m6Acluster=m6Acluster)
PCA.mean=aggregate(PCA[,1:2], list(m6Acluster=PCA$m6Acluster), mean)
pdf(file="./01.step1/07.m6ALIHC_PCA/PCA.pdf", height=5, width=6.5)
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
ggsave(filename = "./01.step1/07.m6ALIHC_PCA/PCA.jpg",plot =p1,
       height=5, width=6.5,dpi = 1000)
}
if(T){
adj.P.Val.Filter=0.001

#读取输入文件，并对输入文件整理
rt=read.table('./00.data/LIHC/megerData/merge.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取cluster文件
cluster=read.table('00.data/LIHC/cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
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
dir <- './00.data/LIHC/PCA/'
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
pdf(file="./01.step1/07.m6ALIHC_PCA/venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#保存交集基因
interGenes=Reduce(intersect,geneList)
write.table(file="./00.data/LIHC/PCA/interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

#保存交集基因的表达量
interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="./00.data/LIHC/PCA/interGeneExp.txt", sep="\t", quote=F, col.names=F)
}
if(T){
pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

rt=read.table("./00.data/LIHC/PCA/interGene.txt", header=F, sep="\t", check.names=F)     #读取输入文件
#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]    
gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="./00.data/LIHC/PCA/GO.txt", sep="\t", quote=F, row.names = F)

#定义显示Term数目
showNum=20
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="./01.step1/07.m6ALIHC_PCA/GO_barplot.pdf", width=10, height=13)
bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="./01.step1/07.m6ALIHC_PCA/GO_bubble.pdf",width=10, height=13)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG, file="./00.data/LIHC/PCA/KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="./01.step1/07.m6ALIHC_PCA/KEGG_barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()

#气泡图
pdf(file="./01.step1/07.m6ALIHC_PCA/KEGG_bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()
}
