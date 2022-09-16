rm(list = ls())  
options(stringsAsFactors = F)#清空环境
load('./01.data/01GSE20966.RData')

logFoldChange=1
adjustP=0.05
rt <- dat

##limma差异分析
modType=c(rep("non−diabetic condition",10),rep("diabetic condition_donor",10))##设置对照组
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)#提取差异信息
diffSig <- allDiff[with(allDiff, (abs(logFC)>=logFoldChange & P.Value < adjustP )), ]#筛选符合的基因2916个
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]#筛选上调的基因1569个
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]#筛选下调的基因1347个

diff_data <- rt[diffSig[,1],]
dim(diff_data)
save(diff_data,file = '02pheatmap_input.RData')

