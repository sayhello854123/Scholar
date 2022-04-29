rm(list = ls())  
options(stringsAsFactors = F)#清空环境
library(limma)
library(data.table)
library(reshape2)
library(ggpubr)
gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F)
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

p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            xlab="",
            legend.title="Type",
            palette = c("blue", "red"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
pdf(file=y, width=7, height=5)
print(p1)
dev.off()
}

##COAD##
x <- '00.data/COAD/01.TCGA/TCGA-COAD.htseq_fpkm.tsv/TCGA-COAD.htseq_fpkm.tsv'
y <- './01.step1/04.m6Adiff/COAD_boxplot.pdf'
diffplot(x,y)
##ESCA##
x <- '00.data/ESCA/01.TCGA/TCGA-ESCA.htseq_fpkm.tsv/TCGA-ESCA.htseq_fpkm.tsv'
y <- './01.step1/04.m6Adiff/ESCA_boxplot.pdf'
diffplot(x,y)
###LIHC###
x <- '00.data/LIHC/TCGA-LIHC.htseq_fpkm.tsv/TCGA-LIHC.htseq_fpkm.tsv'
y <- './01.step1/04.m6Adiff/LIHC_boxplot.pdf'
diffplot(x,y)
###PAAD###
x <- '00.data/PAAD/01.TCGA/TCGA-PAAD.htseq_fpkm.tsv/TCGA-PAAD.htseq_fpkm.tsv'
y <- './01.step1/04.m6Adiff/PAAD_boxplot.pdf'
diffplot(x,y)
##STAD###
x <- '00.data/STAD/01.TCGA/TCGA-STAD.htseq_fpkm.tsv'
y <- './01.step1/04.m6Adiff/STAD_boxplot.pdf'
diffplot(x,y)
##READ###
x <- '00.data/READ/01.TCGA/TCGA-READ.htseq_fpkm.tsv'
y <- './01.step1/04.m6Adiff/READ_boxplot.pdf'
diffplot(x,y)
