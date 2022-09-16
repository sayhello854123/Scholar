rm(list = ls())  
options(stringsAsFactors = F)
dir.create('./03.wgcna')
load('./01.data/03wgcna_input.RData')

library(WGCNA)
#step1: 输入数据的准备和表型文件准备
WGCNA_matrix = t(rt[order(apply(rt,1,mad), decreasing = T),])#矩阵转置并排序
datExpr0 <- WGCNA_matrix 
datExpr <- datExpr0
ann <- list()
ann$id <-rownames(datExpr)
ann$type <- c(rep("non-diabetic condition",10),rep("diabetic condition_donor",10))
ann <- data.frame(ann)
rownames(ann) <- ann[,1]
datTraits <- ann#制作表型文件 
sampleNames = rownames(datExpr)
traitRows = match(sampleNames, datTraits$id) #匹配表型文件 
rownames(datTraits) = datTraits[traitRows, 1] 

##step2确定最佳beta值
if(T){
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)##设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
pdf(file='./03.wgcna/01.pdf',width=5,height=5)
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
}

#step3：构建共表达矩阵
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,#power是最佳beta值10
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
table(net$colors) 

#step4: 模块可视化
if(T){
mergedColors = labels2colors(net$colors)
table(mergedColors)#查看模块的基因数目
pdf(file='./03.wgcna/02.pdf',width=5,height=5)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
}
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#首先针对样本做个系统聚类树
if(T){
datExpr_tree<-hclust(dist(datExpr), method = "average")
pdf(file='./03.wgcna/03.pdf',width=6,height=5)
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
sample_colors <- numbers2colors(as.numeric(factor(datTraits$type)), 
                                colors = c("blue","red","yellow","green"),signed = FALSE)

plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")
dev.off()
}
table(datTraits$type)

#step5:模块和性状的关系
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design=model.matrix(~0+ datTraits$type)
  colnames(design)=levels(as.factor(datTraits$type))
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  pdf(file='./03.wgcna/04.Module-trait-relationships.pdf',width=12,height=7)
  par(mar = c(8, 8.5, 3, 3))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()  
}


###step6:性状的模块的具体基因分析

#首先计算模块与基因的相关性矩阵
if(T){
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#再计算性状与基因的相关性矩阵
T2D = as.data.frame(design[,2])
names(T2D) = "T2D"
geneTraitSignificance = as.data.frame(cor(datExpr, T2D, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(T2D), sep="")
names(GSPvalue) = paste("p.GS.", names(T2D), sep="")
}

#最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
if(T){
module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;
pdf(file='./03.wgcna/05.Gene significance for T2DM.pdf',width=5,height=5)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for T2D",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

####step7:提取指定模块的基因名
module = "black";
# Select module probes
probes = colnames(datExpr) 
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
head(modProbes)
which.module="black";
dat=datExpr[,moduleColors==which.module ] 
}
