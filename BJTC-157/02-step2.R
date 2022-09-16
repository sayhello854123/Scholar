library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)
load('./00.data/03.LIHC/wgcna_input.RData')
dge=DGEList(counts=rt)
cpmValue=cpm(dge)
keep=rowMeans(cpmValue)>1
y=dge[keep, , keep.lib.sizes=FALSE]
y=calcNormFactors(y)
logCPM=cpm(y, log=TRUE, prior.count=3)

#输出RPKM值
geneRT=read.table("./00.data/03.LIHC/geneLength.txt",sep="\t",header=T,check.names=F)
geneRT=geneRT[,c(1,2,2)]
geneRT=as.matrix(geneRT)
rownames(geneRT)=geneRT[,1]
exp1=geneRT[,2:ncol(geneRT)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)
same <- intersect(row.names(y$counts),row.names(data1))
length(same)
data1 <- data1[same,]
y$counts <- y$counts[same,]
data1=data1[row.names(y$counts),]
geneLen=data1[,1]
rpkm=rpkm(y,geneLen)
dim(rpkm)
save(rpkm,file = './00.data/03.LIHC/wgcna_input1.RData')
#rpkm=rbind(id=colnames(rpkm),rpkm)
#write.table(rpkm,file="TCGA_rpkm.txt",sep="\t",quote=F,col.names=F)