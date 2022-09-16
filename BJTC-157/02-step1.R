rm(list = ls())  
options(stringsAsFactors = F)#清空环境
library(limma)
luad<-data.table::fread(file = './00.data/03.LIHC/TCGA-LIHC.htseq_counts.tsv',data.table = F)
rownames(luad) <- luad[,1]
luad <- luad[,-1]
data <- 2^luad-1
gencode <- data.table::fread('./00.data/03.LIHC/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
same <- intersect(row.names(gencode),row.names(data))
length(same)
data1 <- cbind(gencode[same,],data[same,])
data1 <- data1[,-c(1,3:6)]
write.table(data1,file = './00.data/03.LIHC/data.txt',sep = '\t',quote = F)

library(data.table)
library(survival)
library(survminer)
cli <-data.table::fread('./00.data/03.LIHC/TCGA-LIHC.survival.tsv',data.table = F)
row.names(cli) <- cli[,1]
cli$futime <- cli$OS.time/365
cli$fustat <- cli$OS
cli <- cli[,-c(1:4)]
pFilter=0.05 
rt <- read.table('./00.data/03.LIHC/CIBERSOFT.txt',sep = '\t',header = T,row.names = 1)
data=rt[rt[,"P.value"]<pFilter,]
data=data[,1:(ncol(rt)-4)]
#rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
table(group)
group=sapply(strsplit(group,""), "[", 1)
table(group)
rt1=data[group==0,]
sameSample=intersect(row.names(rt1), row.names(cli))
rt=cbind(cli[sameSample,], data[sameSample,])
#####生存分析####
picDir='./02.bulk'
dir.create(picDir)
var="Macrophages.M2"
rt=rt[,c("futime","fustat",var)]
colnames(rt)=c("futime","fustat","var")

#获取最优cutoff
res.cut=surv_cutpoint(rt, time = "futime", event = "fustat",variables =c("var"))
res.cut
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

#绘制
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
pdf(file='02.bulk/M2.pdf',onefile = FALSE,width = 6,height =5)
print(surPlot)
dev.off()
data1 <- read.table('./00.data/03.LIHC/data.txt',sep = '\t',header = T,check.names = F)
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
table(res.cat$var)
same <- intersect(row.names(res.cat),colnames(data))
rt <- data[,same]
rt2 <- rt[,res.cat$var=="high"]#138
dim(rt2)
rt3 <- rt[,res.cat$var=="low"]#202
dim(rt3)
rt <- cbind(rt3,rt2)
dim(rt)
save(rt,file = './00.data/03.LIHC/wgcna_input.RData')
