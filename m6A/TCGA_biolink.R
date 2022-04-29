library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
query <- GDCquery(project = 'TCGA-STAD', 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - FPKM")
GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query)
FPKM_matrix=assay(expdat)
gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- do.call(rbind,strsplit(gencode$id,'\\.'))[,1]
fpkmToTpm=function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(FPKM_matrix, 2, fpkmToTpm)
same <- intersect(row.names(gencode),row.names(tpm))
length(same)
data1 <- cbind(gencode[same,],tpm[same,])
head(data1)[1:5,1:5]
data1 <- data1[,-c(1,3:6)]
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(rt))
rt1=rt[sameGene,]
group=sapply(strsplit(colnames(rt1),"\\-"), "[", 4)
table(group)
group=sapply(strsplit(group,""), "[", 1)
table(group)
Tumor <- rt1[,group == 0]
Normal <- rt1[,group == 1]
data <- cbind(Normal,Tumor)
write.csv(data,file = './00.data/STAD/01.TCGA/33.csv')
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])
sampleType=c(rep(1,conNum), rep(2,treatNum))


exp=log2(data+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
write.csv(exp,file = './00.data/STAD/01.TCGA/33.csv')
exp <- read.csv('./00.data/STAD/01.TCGA/33.csv',header = T,row.names = 1)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
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
ggsave(filename = './06.step6/diff_plot.pdf',plot = p1,width=7, height=5,dpi = 1000)
ggsave(filename = './06.step6/diff_plot.jpg',plot = p1,width=7, height=5,dpi = 1000)