library(tidyverse)
library(devtools)
library(AnnoProbe)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 172157430)#GSE102238
gset <- getGEO('GSE133057', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
dat1=exprs(a)
dim(dat1)
gpl='GPL6102'
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
metadata=pData(a)
clinical = data.frame(gsm=metadata[,2],
                      futime=trimws(sapply(as.character(metadata$characteristics_ch1.3),function(x) strsplit(x,":")[[1]][2])),
                      fustat=trimws(sapply(as.character(metadata$characteristics_ch1.4),function(x) strsplit(x,":")[[1]][2])),
                      Age = trimws(sapply(as.character(metadata$characteristics_ch1.5),function(x) strsplit(x,":")[[1]][2])),
                      Gender = trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2]))
                      
                      
                      
)
rownames(clinical) <- clinical[,1]
gene=read.table('./00.data/m6A_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(rt1))
outTab=rbind(geneNames=colnames(rt1), rt1)
write.table(outTab, file="./00.data/READ/02.GSE17536/GSE133057.txt", sep="\t", quote=F, col.names=F)
outTab1=rbind(id_name=colnames(clinical), clinical)
write.table(outTab1, file="./00.data/READ/02.GSE17536/clinical1.txt", sep="\t", quote=F, col.names=F)
