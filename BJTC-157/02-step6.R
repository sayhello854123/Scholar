library(glmnet)
library(survival)
library(survivalROC)
library(survminer)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
pvalueFilter=0.05      
qvalueFilter=0.05 
load('./02.bulk/survival/cox_input.RData')
if(T){
pFilter=0.05
outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
      sigGenes=c(sigGenes,gene)
      outTab=rbind(outTab,
                   cbind(gene=gene,
                         #KM=pValue,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         coxPvalue=coxP) )
    }
  }
surSigExp=rt[,sigGenes]
rt <- surSigExp
rt$futime[rt$futime<=0]=0.003

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)

trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore))
predictTime=3      #1年的ROC曲线，需要做3年或5年改成相应的数值
roc=survivalROC(Stime=outTab$futime, status=outTab$fustat, marker = outTab$riskScore, predict.time =predictTime, method="KM")
sum=roc$TP-roc$FP
cutOp=roc$cut.values[which.max(sum)]
cutTP=roc$TP[which.max(sum)]
cutFP=roc$FP[which.max(sum)]
pdf(file="./02.bulk/survival/ROC.pdf",width=5.5,height=5.5)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="black", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.2, cex.lab=1.2, cex.axis=1.2, font=1.2)
points(cutFP,cutTP, pch=20, col="red",cex=1.5)
text(cutFP+0.1,cutTP-0.05,paste0("Cutoff:",sprintf("%0.3f",cutOp)))
dev.off()
risk=as.vector(ifelse(trainScore>cutOp,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="./02.bulk/survival/tcgaRisk.txt",sep="\t",quote=F,row.names=F)

#根据公式计算test组风险值,输出test组风险值结果
load('./02.bulk/survival/GEO_input.RData')
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>cutOp,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="./02.bulk/survival/geoRisk.txt",sep="\t",quote=F,row.names=F)
}


bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile,header=T,sep="\t")
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  HR = (diff$obs[2]/diff$exp[2])/(diff$obs[1]/diff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%0.3f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     title='TCGA cohort',
                     #conf.int=TRUE,
                     pval=pValue,
                     pval.size=6,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
  
  #保存输出图片
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="./02.bulk/survival/tcgaRisk.txt",outFile="./02.bulk/survival/tcgaRisk.pdf")
bioSurvival(inputFile="./02.bulk/survival/geoRisk.txt",outFile="./02.bulk/survival/geoRisk.pdf")

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
rt=read.table("./02.bulk/survival/gene.txt", header=F, sep="\t", check.names=F)  
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="./02.bulk/survival/GO_barplot.pdf", width=20, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="./02.bulk/survival/GO_bubble.pdf", width=20, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]


#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="./02.bulk/survival/KEGG_barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()

#气泡图
pdf(file="./02.bulk/survival/KEGG_bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()