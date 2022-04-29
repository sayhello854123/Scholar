library(GEOquery)
library(limma)
library(tidyverse)
library(sva)
library(AnnoProbe)
library(tidyverse)

files=c("./00.data/STAD/01.TCGA/TCGA.TPM.txt", "./00.data/STAD/03.GEO/GSE15459.txt",
        "./00.data/STAD/03.GEO/GSE34942.txt")
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit((strsplit(inputFile, "\\/|\\-")[[1]][5]), "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)
#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit((strsplit(inputFile, "\\/|\\-")[[1]][5]), "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对TCGA删除正常样品
  if(header[1] == "TCGA"){
    group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
    group=sapply(strsplit(group,""), "[", 1)
    rt=rt[,group==0]
    rt=t(rt)
    #row.names(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(rt))
    rt=avereps(rt)
    rt=t(rt)
  }
  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] != "TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}
#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="./00.data/STAD/04.meger/merge.txt", sep="\t", quote=F, col.names=F)
gene=read.table('./00.data/m6A_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(outTab))
geneExp=outTab[sameGene,]
tpmOut=rbind(ID=colnames(geneExp), geneExp)
write.table(tpmOut, file="./00.data/STAD/04.meger/m6A_merge.txt", sep="\t", col.names=F, quote=F)
