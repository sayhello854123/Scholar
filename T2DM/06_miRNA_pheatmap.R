rm(list = ls())  
options(stringsAsFactors = F)#清空环境
setwd('../')
library(limma)
library(pheatmap)
load('./01.data/06.miRNA_input.RData')
dir.create('09.miRNA_Pheatmap')

#整理miRNA数据
if(T){
  logFoldChange=1
  adjustP=0.05
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data <- data[,-c(11:17)]
  ##进行差异分析
  rt <- data
  modType=c(rep("normal",10),rep("tumor",9))
  design <- model.matrix(~0+factor(modType))
  colnames(design) <- c("con","treat")
  fit <- lmFit(rt,design)
  cont.matrix<-makeContrasts(treat-con,levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  ##提取差异
  diffSig <- allDiff[with(allDiff, (abs(logFC)>=logFoldChange & P.Value < adjustP )), ]
  write.table(diffSig,file="./09.miRNA_Pheatmap/miRNA_diff.xls",sep="\t",quote=F)
  diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]
  write.table(diffUp,file="./09.miRNA_Pheatmap/up.xls",sep="\t",quote=F)
  diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]
  write.table(diffDown,file="./09.miRNA_Pheatmap/down.xls",sep="\t",quote=F)
  same <- intersect(row.names(diffSig),row.names(rt))
  diff_data <- rt[same,]
}

#进行差异热图绘制
if(T){
  Type=c(rep("non−diabetic condition",10),rep("diabetic condition_donor",9))
  names(Type)=colnames(rt)
  Type=as.data.frame(Type)
  pdf(file='./09.miRNA_Pheatmap/miRNA_pheatmap.pdf',width=8,height=10)
  pheatmap(diff_data,
           annotation=Type,
           cluster_cols = T,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           show_colnames = T,
           scale="row",  #矫正
           #border_color ="NA",
           fontsize = 8,
           fontsize_row=6,
           fontsize_col=6)
  dev.off()
}
