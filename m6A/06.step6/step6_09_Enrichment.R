library(limma)
library(ggplot2)
library(VennDiagram)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
if(T){
  pvalueFilter=0.05       #p值过滤条件
  qvalueFilter=0.05       #矫正后的p值过滤条件
  #定义颜色
  colorSel="qvalue"
  if(qvalueFilter>0.05){
    colorSel="pvalue"
  }
  rt=read.table("./00.data/STAD/09.VENN/interGene.txt", header=F, sep="\t", check.names=F)     #读取输入文件
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
  write.table(GO, file="./00.data/STAD/10.Enrichment_analysis/GO.txt", sep="\t", quote=F, row.names = F)
  
  #定义显示Term数目
  showNum=10
  if(nrow(GO)<30){
    showNum=nrow(GO)
  }
  
  #柱状图
  pdf(file="./06.step6/08.Enrichment/GO_barplot.pdf", width=8, height=7)
  bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
  print(bar)
  dev.off()
  
  #气泡图
  pdf(file="./06.step6/08.Enrichment/GO_bubble.pdf",width=8, height=7)
  bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
  print(bub)
  dev.off()
  
  
  #kegg富集分析
  kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
  KEGG=as.data.frame(kk)
  KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
  KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
  #保存富集结果
  write.table(KEGG, file="./00.data/STAD/10.Enrichment_analysis/KEGG.txt", sep="\t", quote=F, row.names = F)
  
  #定义显示通路的数目
  showNum=30
  if(nrow(KEGG)<showNum){
    showNum=nrow(KEGG)
  }
  
  #柱状图
  pdf(file="./06.step6/08.Enrichment/KEGG_barplot.pdf", width=8, height=6)
  barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
  dev.off()
  
  #气泡图
  pdf(file="./06.step6/08.Enrichment/KEGG_bubble.pdf", width = 8, height = 6)
  dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
  dev.off()
}
