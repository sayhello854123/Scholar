rm(list = ls())  
options(stringsAsFactors = F)#清空环境
gc()
library(limma)
library(data.table)
library(reshape2)
library(ggpubr)
gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
geneRT=read.table("./00.data/m5c_gene.txt", header=T, sep="\t", check.names=F)
diffplot <- function(x,y){
  rt<-data.table::fread(file = x,data.table = F)
  rownames(rt) <- rt[,1]
  rt <- rt[,-1]
  data <- 2^rt-1
  same <- intersect(row.names(gencode),row.names(data))
  length(same)
  data1 <- cbind(gencode[same,],data[same,])
  data1 <- data1[,-c(1,3:6)]
  rt=as.matrix(data1)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  
  #FPKM转换为TPM
  fpkmToTpm=function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  tpm=apply(rt, 2, fpkmToTpm)
  
  sameGene=intersect(as.vector(geneRT[,1]), row.names(tpm))
  data=tpm[sameGene,]
  ##
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  table(group)
  group_list=ifelse(as.numeric(substr(colnames(data),14,15)) < 10,'tumor','normal')
  table(group_list)
  Tumor <- data[,group_list == "tumor"]
  Normal <- data[,group_list == "normal"]
  data <- cbind(Normal,Tumor)
  data <- read.csv('./1.csv',header = T,row.names = 1)
  conNum=length(colnames(Normal))       
  treatNum=length(colnames(Tumor))     
  sampleType=c(rep(1,conNum), rep(2,treatNum))
  
  #把数据转换成ggplot2输入文件
  exp=log2(data+1)
  exp=as.data.frame(t(exp))
  exp=cbind(exp, Type=sampleType)
  exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
  data=melt(exp, id.vars=c("Type"))
  colnames(data)=c("Type", "Gene", "Expression")
  
  p=ggboxplot(data, x="Gene", y="Expression", color = "black", fill="Type",
              ylab="Gene expression",
              xlab="",
              legend.title="Type",
              palette="aaas",
              width=1)
  p=p+rotate_x_text(60)
  p1=p+stat_compare_means(aes(group=Type),
                          method="wilcox.test",
                          symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                          label = "p.signif")
  ggsave(filename = y,p1,width = 8,height = 6,dpi = 1000)
  ggsave(filename = y1,p1,width = 8,height = 6,dpi = 1000)
}
###LIHC###
x <- '00.data/TCGA-LIHC.htseq_fpkm.tsv'
y <- './02.m5Cdiff/LIHC_boxplot.pdf'
y1 <- './02.m5Cdiff/LIHC_boxplot.jpg'
diffplot(x,y)



