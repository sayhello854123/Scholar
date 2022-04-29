#install_github("jmzeng1314/AnnoProbe")
library(tidyverse)
library(devtools)
library(AnnoProbe)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 524288 * 2)#GSE102238
gset <- getGEO('GSE53625', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
dat1=exprs(a)
dim(dat1)
gpl='GPL18109'
ids=idmap(gpl,type = 'pipe')
colnames(ids)=c('probe_id','symbol')
ids<- ids[!duplicated(ids$symbol),]
ids=ids[ids$symbol != '',]
rownames(ids) <- ids[,1]
same <- intersect(row.names(ids),row.names(dat1))
data <- cbind(ids[same,],dat1[same,])
data <- data[,-1]
rt=as.matrix(data)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
tpmOut=rbind(ID=colnames(data), data)
#write.table(tpmOut, file="./00.data/ESCA/02.GSE53625/GSE53625.txt", sep="\t", col.names=F, quote=F)

metadata=pData(a)
clinical = data.frame(gsm=metadata[,2],
                      #id = metadata[,1],
                      tissue=metadata[,8],
                      futime=trimws(sapply(as.character(metadata$characteristics_ch1.15),function(x) strsplit(x,":")[[1]][2])),
                      fustat=trimws(sapply(as.character(metadata$characteristics_ch1.14),function(x) strsplit(x,":")[[1]][2])),
                      Age=trimws(sapply(as.character(metadata$characteristics_ch1.1),function(x) strsplit(x,":")[[1]][2])),
                      Gender=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2])),
                      Grade = trimws(sapply(as.character(metadata$characteristics_ch1.6),function(x) strsplit(x,":")[[1]][2])),
                      Stage = trimws(sapply(as.character(metadata$characteristics_ch1.9),function(x) strsplit(x,":")[[1]][2])),
                      T = trimws(sapply(as.character(metadata$characteristics_ch1.7),function(x) strsplit(x,":")[[1]][2])),
                      N= trimws(sapply(as.character(metadata$characteristics_ch1.8),function(x) strsplit(x,":")[[1]][2]))
                      
)
group <- grep('cancer.*?',clinical$tissue)
clinical <- clinical[group,]
table(clinical$Grade)
clinical$Grade <- gsub('moderately','G2',clinical$Grade)
clinical$Grade <- gsub('poorly','G3',clinical$Grade)
clinical$Grade <- gsub('well','G1',clinical$Grade)
clinical$fustat <- ifelse(clinical$fustat=="yes",1,0)
rownames(clinical) <- clinical[,1]
rt <- read.table('./00.data/ESCA/02.GSE53625/GSE53625.txt',sep = '\t',header = T)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
same <- intersect(row.names(clinical),colnames(data))
rt <- data[,same]
tpmOut=rbind(ID=colnames(rt), rt)
write.table(tpmOut, file="./00.data/ESCA/02.GSE53625/GSE53625.txt", sep="\t", col.names=F, quote=F)
tpmOut1=rbind(ID=colnames(clinical), clinical)
write.table(tpmOut1, file="./00.data/ESCA/02.GSE53625/clinical.txt", sep="\t", col.names=F, quote=F)

rt <- read.table('./00.data/ESCA/02.GSE53625/GSE53625.txt',sep = '\t',header = T,row.names = 1,check.names = F)
rt1 <- rt[,c(1:60)]
colnames(rt1)
rt2 <- rt[,c(61:179)]
tpmOut=rbind(ID=colnames(rt1), rt1)
write.table(tpmOut, file="./00.data/ESCA/02.GSE53625/GSE53622.txt", sep="\t", col.names=F, quote=F)
tpmOut1=rbind(ID=colnames(rt2), rt2)
write.table(tpmOut1, file="./00.data/ESCA/02.GSE53625/GSE53624.txt", sep="\t", col.names=F, quote=F)
