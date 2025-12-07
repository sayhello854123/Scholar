library(tidyverse)
library(devtools)
library(AnnoProbe)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 172157430)#GSE102238
gset <- getGEO('GSE39582', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
dat1=a@assayData$exprs
dim(dat1)
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
metadata=pData(a)
clinical = data.frame(gsm=metadata[,2],
                      Type = trimws(sapply(as.character(metadata$characteristics_ch1.1),function(x) strsplit(x,":")[[1]][2])),
                      futime= trimws(sapply(as.character(metadata$characteristics_ch1.14),function(x) strsplit(x,":")[[1]][2])),
                      fustat=trimws(sapply(as.character(metadata$characteristics_ch1.13),function(x) strsplit(x,":")[[1]][2])),
                      Sex=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2])),
                      Age= trimws(sapply(as.character(metadata$characteristics_ch1.3),function(x) strsplit(x,":")[[1]][2])),
                      tnm_stage = trimws(sapply(as.character(metadata$characteristics_ch1.4),function(x) strsplit(x,":")[[1]][2])),
                      tnm_T =trimws(sapply(as.character(metadata$characteristics_ch1.5),function(x) strsplit(x,":")[[1]][2])),
                      tnm_N =trimws(sapply(as.character(metadata$characteristics_ch1.6),function(x) strsplit(x,":")[[1]][2])),
                      tnm_M =trimws(sapply(as.character(metadata$characteristics_ch1.7),function(x) strsplit(x,":")[[1]][2]))
)

write.csv(clinical,file = './step00-data/GSE39582_cli.csv')
table(clinical$Type)
clinical_uninfected <- subset(clinical, Type %in% c( "Non Tumoral"))
rownames(clinical_uninfected) = clinical_uninfected$gsm
same1 =intersect(row.names(clinical_uninfected),colnames(rt1))
normal =rt1[,same1]


clinical_Sepsis <- subset(clinical, Type %in% c( "discovery",'validation'))
rownames(clinical_Sepsis) = clinical_Sepsis$gsm
same2 =intersect(row.names(clinical_Sepsis),colnames(rt1))
Sepsis= rt1[,same2]
data2 = cbind(normal,Sepsis)#19 566
write.csv(data2,file = './step00-data/GSE39582.csv')



rt = read.csv('./step00-data/GSE39582.csv',header = T,row.names = 1)
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)
conNum=19
treatNum=566
Type=c(rep("Control",conNum),rep("CRC",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0('./step01-normalize/GSE39582',".normalize.txt"),sep="\t",quote=F,col.names=F)