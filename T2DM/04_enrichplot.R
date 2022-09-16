rm(list = ls())  
options(stringsAsFactors = F)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
dir.create('04.enrichplot')

pvalueFilter=0.05         #p值过滤条件
qvalueFilter=0.05

load('01.data/04.enrichplot_input.RData')

rt=rt[is.na(rt[,"entrezID"])==F,]                            #去除基因id为NA的基因
gene=rt$entrezID
geneFC=2^rt$logFC
names(geneFC)=gene

if(T){
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

#GO富集分析
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<0.2),]
#保存富集结果
write.table(GO,file="./04.enrichplot/GO.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="./04.enrichplot/barplot.pdf",width = 9,height = 7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="./04.enrichplot/bubble.pdf",width = 9,height = 7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#圈图
pdf(file="./04.enrichplot/circos.pdf",width = 10,height = 8)
cnet=cnetplot(kk, foldChange=geneFC, showCategory = 5, circular = TRUE, colorEdge = TRUE)
print(cnet)
dev.off()


colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<0.56),]
#保存富集结果
write.table(KEGG,file="./04.enrichplot/KEGG.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="./04.enrichplot/KEGG_barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

#气泡图
pdf(file="./04.enrichplot/KEGG_bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

#圈图
pdf(file="./04.enrichplot/KEGG_circos.pdf",width = 11,height = 7)
kkx=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, foldChange=geneFC,showCategory = 5, circular = TRUE, colorEdge = TRUE,node_label="all")
dev.off()
}
