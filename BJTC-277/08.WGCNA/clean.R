rm(list = ls())
options(stringsAsFactors = F) 
gc()
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(WGCNA)
library(survival)
library(survminer)
library(limma)
library(edgeR)

##data clean
stad<-data.table::fread(file = './08.WGCNA/TCGA-STAD.htseq_fpkm.tsv',data.table = F)
rownames(stad) <- stad[,1]
stad <- stad[,-1]
data <- 2^stad-1
gencode <- data.table::fread('./08.WGCNA/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
same <- intersect(row.names(gencode),row.names(data))
length(same)
data1 <- cbind(gencode[same,],data[same,])
data1 <- data1[,-c(1,3:6)]
write.table(data1,file = './08.WGCNA/data.txt',sep = '\t',quote = F)
data1 <- read.table('./08.WGCNA/data.txt',sep = '\t',header = T,check.names = F)
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
save(data,file = './08.WGCNA/wgcna_input.RData')




###
cli <-data.table::fread('./08.WGCNA/TCGA-STAD.survival.tsv',data.table = F)
row.names(cli) <- cli[,1]
cli$futime <- cli$OS.time/365
cli$fustat <- cli$OS
cli <- cli[,-c(1:4)]
pFilter=0.05 
rt <- read.table('./08.WGCNA/CIBERSOFT.txt',sep = '\t',header = T,row.names = 1)
data=rt[rt[,"P.value"]<pFilter,]
data=data[,1:(ncol(rt)-4)]
group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
table(group)
group=sapply(strsplit(group,""), "[", 1)
table(group)
rt1=data[group==0,]
sameSample=intersect(row.names(rt1), row.names(cli))
rt=cbind(cli[sameSample,], data[sameSample,])
var="Macrophages.M2"
rt=rt[,c("futime","fustat",var)]
colnames(rt)=c("futime","fustat","var")

#获取最优cutoff
res.cut=surv_cutpoint(rt, time = "futime", event = "fustat",variables =c("var"))
res.cut#0.2307188
res.cat=surv_categorize(res.cut)
fit=survfit(Surv(futime, fustat) ~var, data = res.cat)

#比较高低表达生存差异
diff=survdiff(Surv(futime, fustat) ~var,data =res.cat)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
surPlot=ggsurvplot(fit, 
                   data=res.cat,
                   #conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf(file='./08.WGCNA/M2.pdf',onefile = FALSE,width = 8,height =6)
print(surPlot)
dev.off()


####样本划分
load('./08.WGCNA/wgcna_input.RData')
table(res.cat$var)
same <- intersect(row.names(res.cat),colnames(data))
data <- data[,same]
dim(data)
meta <- data.frame(colnames(data),res.cat[same,])
rt2 <- data[,res.cat$var=="high"]#47
dim(rt2)
rt3 <- data[,res.cat$var=="low"]#302
dim(rt3)
save(data,file = './08.WGCNA/wgcna_input.RData')
save(meta,file = './08.WGCNA/wgcna_ann.RData')



