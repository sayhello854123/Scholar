library(tidyverse)
library(devtools)
library(AnnoProbe)
library(GEOquery)
if(T){
Sys.setenv("VROOM_CONNECTION_SIZE" = 172157430)#GSE102238
gset <- getGEO('GSE15459', destdir=".",
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
outTab=rbind(geneNames=colnames(rt1), rt1)
write.table(outTab, file="./00.data/STAD/03.GEO/GSE15459.txt", sep="\t", quote=F, col.names=F)
}
if(T){
Sys.setenv("VROOM_CONNECTION_SIZE" = 172157430)#GSE102238
gset <- getGEO('GSE34942', destdir=".",
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
outTab=rbind(geneNames=colnames(rt1), rt1)
write.table(outTab, file="./00.data/STAD/03.GEO/GSE34942.txt", sep="\t", quote=F, col.names=F)
}