
library(tidyverse)
library(devtools)
library(AnnoProbe)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 172157430)#GSE102238
gset <- getGEO('GSE5081', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
dat1=exprs(a)
dim(dat1)
#metadata=pData(a)
gpl='GPL570'
ids=idmap(gpl)
colnames(ids)=c('probe_id','symbol')
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat1),]
dat1=dat1[ids$probe_id,] 
ids$median=apply(dat1,1,median) 
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat1=dat1[ids$probe_id,]
rownames(dat1)=ids$symbol
rt1 <- dat1
rt1=normalizeBetweenArrays(rt1)
metadata=pData(a)
clinical = data.frame(gsm=metadata[,2],
                     
            tissue=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2]))
)                     
outTab=rbind(geneNames=colnames(rt1), rt1)
write.table(outTab, file="./00.data/STAD/02.GSE5081/GSE5081.txt", sep="\t", quote=F, col.names=F)
outTab1=rbind(id_name=colnames(clinical), clinical)
write.table(outTab1, file="./00.data/STAD/02.GSE5081/clinical.txt", sep="\t", quote=F, col.names=F)

if(T){
rt=read.table("./00.data/STAD/02.GSE5081/GSE5081.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

gene=read.table('./00.data/m6A_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), row.names(data))
data=data[sameGene,]

erosion <- data[,c(1:16)]
sampleType=c(rep(1,8), rep(2,8))
exp=log2(erosion+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "HP-", "HP+")
erosion=melt(exp, id.vars=c("Type"))
colnames(erosion)=c("Type", "Gene", "Expression")

p=ggboxplot(erosion, x="Gene", y="Expression", color = "Type", 
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
p1 <- p1+ggtitle("Site of erosion") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "./00.data/STAD/02.GSE5081/erosion.pdf",plot =p1,
       width=8, height=6.5,dpi = 1000)
ggsave(filename = "./00.data/STAD/02.GSE5081/erosion.jpg",plot =p1,
       width=8, height=6.5,dpi = 1000)
}
if(T){
rt=read.table("./00.data/STAD/02.GSE5081/GSE5081.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

gene=read.table('./00.data/m6A_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), row.names(data))
data=data[sameGene,]

non_erosive <- data[,c(17:32)]
sampleType=c(rep(1,8), rep(2,8))
exp=log2(non_erosive+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "HP+", "HP-")
non_erosive=melt(exp, id.vars=c("Type"))
colnames(non_erosive)=c("Type", "Gene", "Expression")

p=ggboxplot(non_erosive, x="Gene", y="Expression", color = "Type", 
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
p1 <- p1+ggtitle("Site of non_erosive") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "./00.data/STAD/02.GSE5081/non_erosive.pdf",plot =p1,
       width=8, height=6.5,dpi = 1000)
ggsave(filename = "./00.data/STAD/02.GSE5081/non_erosive.jpg",plot =p1,
       width=8, height=6.5,dpi = 1000)
}
if(T){
rt=read.table("./00.data/STAD/02.GSE5081/GSE5081.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

gene=read.table('./00.data/m6A_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), row.names(data))
data=data[sameGene,]
HP1 <-data[,c(1:8,25:32)] 
HP2 <-data[,c(9:24)]  
HP <- cbind(HP1,HP2)

sampleType=c(rep(1,16), rep(2,16))
exp=log2(HP+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "HP-", "HP+")
HP=melt(exp, id.vars=c("Type"))
colnames(HP)=c("Type", "Gene", "Expression")

p=ggboxplot(HP, x="Gene", y="Expression", color = "Type", 
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
p1 <- p1+ggtitle("HP- vs HP+ ") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "./00.data/STAD/02.GSE5081/HP.pdf",plot =p1,
       width=8, height=6.5,dpi = 1000)
ggsave(filename = "./00.data/STAD/02.GSE5081/HP.jpg",plot =p1,
       width=8, height=6.5,dpi = 1000)
}
if(T){
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

  rt=read.table('./00.data/STAD/02.GSE5081/GSE5081.txt', header=T, sep="\t", check.names=F)
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
  ssgseaScore=gsva(data, cellMarker, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  ssgseaScore1=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(ssgseaScore1), ssgseaScore1)
  write.table(ssgseaOut,file="./00.data/STAD/02.GSE5081/ssGSEA.result.txt",sep="\t",quote=F,col.names=F)
  ann=read.table('./00.data/STAD/02.GSE5081/group.txt',header=T,sep="\t",row.names=1,check.names=F)
  pdf(file="./00.data/STAD/02.GSE5081/SSgsea.pdf",width=8,height=5.5)
   pheatmap(ssgseaScore1, scale = "row", show_colnames = F, annotation=ann,cluster_cols =F,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
   dev.off()
}

rt=read.table("./00.data/STAD/02.GSE5081/GSE5081.txt", header=T, sep="\t", check.names=F)
riskScore=predict(multiCox,type="risk",newdata=rt)
