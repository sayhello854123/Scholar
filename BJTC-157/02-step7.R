library(GEOquery)
library(limma)
gset <- getGEO('GSE14520', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
dat1=exprs(a)
dim(dat1)

#转换ID#GPL3921	
library(hthgu133a.db)
ids=toTable(hthgu133aSYMBOL)
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
##GPL571
library(hgu133a2.db)
b=gset[[2]]
dat2=exprs(b)
id2=toTable(hgu133a2SYMBOL)
colnames(id2)=c('probe_id','symbol')  
id2=id2[id2$symbol != '',]
id2=id2[id2$probe_id %in%  rownames(dat2),]
dat2=dat2[id2$probe_id,] 
id2$median=apply(dat2,1,median) 
id2=id2[order(id2$symbol,id2$median,decreasing = T),]
id2=id2[!duplicated(id2$symbol),]
dat2=dat2[id2$probe_id,]
rownames(dat2)=id2$symbol
rt2 <- dat2
same <- intersect(row.names(rt2),row.names(rt1))
length(same)
rt <- cbind(rt2[same,],rt1[same,])
###临床数据整理##
cli <- read.table('./00.data/03.GEO/GSE14520_Extra_Supplement.txt',sep = '\t',header = T,row.names = 1,check.names = F)
clinical <- cli[,c(1,8,9)]
table(clinical$`Tissue Type`)
clinical <- clinical[clinical$`Tissue Type`=='Tumor',]
clinical <- na.omit(clinical)
clinical$futime <- clinical$`Survival months`
clinical$fustat <- clinical$`Survival status`  
clinical <- clinical[,-c(1:3)]
clinical$futime <- clinical$futime/12
#rt1 <- t(rt)
same1 <- intersect(row.names(clinical),colnames(rt))
rt <- rt[,same1]




gset <- getGEO('GSE76427', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
c=gset[[1]]
dat3=exprs(c)
dim(dat3)
library(illuminaHumanv4.db)
ids=toTable(illuminaHumanv4SYMBOL)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat3),]
dat3=dat3[ids$probe_id,] 
ids$median=apply(dat3,1,median) 
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat3=dat3[ids$probe_id,]
rownames(dat3)=ids$symbol

metadata=pData(c)
cli = data.frame(gsm=metadata[,2],
                 Type =trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2])),
                 futime=trimws(sapply(as.character(metadata$characteristics_ch1.4),function(x) strsplit(x,":")[[1]][2])),
                 fustat=trimws(sapply(as.character(metadata$characteristics_ch1.3),function(x) strsplit(x,":")[[1]][2]))
)
rownames(cli) <- cli[,1]
table(cli$Type)
cli <- cli[cli$Type=='primary hepatocellular carcinoma tumor',]
same2 <- intersect(row.names(cli),colnames(dat3))
dat3 <- dat3[,same2]
same4 <- intersect(row.names(dat3),row.names(rt))
length(same4)
rt3 <- cbind(dat3[same4,],rt[same4,])
rt3=normalizeBetweenArrays(rt3)

##整合临床数据
cli <- cli[,-c(1,2)]
write.table(cli,file = 'cli1.txt',sep = '\t',quote = F)
write.table(clinical,file = 'cli2.txt',sep = '\t',quote = F)
##合并
cli2 <- read.table('./cli1.txt',sep = '\t',header = T,row.names = 1,check.names = F)
rt4 <- t(rt3)
same5 <- intersect(row.names(cli2),row.names(rt4))
rt <- cbind(cli2[same5,],rt4[same5,])
save(rt,file = './00.data/03.GEO/GEO_input.RData')