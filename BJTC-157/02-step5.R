rm(list = ls())  
options(stringsAsFactors = F)#清空环境
library(limma)
luad<-data.table::fread(file = '02.bulk/survival/TCGA-LIHC.htseq_fpkm.tsv',data.table = F)
rownames(luad) <- luad[,1]
luad <- luad[,-1]
data <- 2^luad-1
gencode <- data.table::fread('./00.data/03.LIHC/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
same <- intersect(row.names(gencode),row.names(data))
length(same)
data1 <- cbind(gencode[same,],data[same,])
data1 <- data1[,-c(1,3:6)]
write.table(data1,file = '02.bulk/survival/data.txt',sep = '\t',quote = F)

data1 <- read.table('./02.bulk/survival/data.txt',sep = '\t',header = T,check.names = F)
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
dim(data)
gene <- read.table('./02.bulk/survival/intersect.txt',sep = '\t',header = T,check.names = F)
same <- intersect(row.names(data),gene[,1])
data2 <- data[same,]
group_list=ifelse(as.numeric(substr(colnames(data2),14,15)) < 10,'tumor','normal')
table(group_list)
rt1 <- t(data2[,group_list == "tumor" ])

cli <-data.table::fread('./00.data/03.LIHC/TCGA-LIHC.survival.tsv',data.table = F)
row.names(cli) <- cli[,1]
cli$futime <- cli$OS.time/365
cli$fustat <- cli$OS
cli <- cli[,-c(1:4)]
sameSample=intersect(row.names(rt1), row.names(cli))
rt=cbind(cli[sameSample,], rt1[sameSample,])
save(rt,file = './02.bulk/survival/cox_input.RData')
