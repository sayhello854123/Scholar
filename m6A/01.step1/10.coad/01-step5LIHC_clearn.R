library(GEOquery)
library(limma)
library(tidyverse)
library(sva)

Sys.setenv("VROOM_CONNECTION_SIZE" = 524288 * 2)
gset <- getGEO('GSE76427', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
dat1=exprs(a)
dim(dat1)
metadata=pData(a)
clinical = data.frame(gsm=metadata[,2],
                      tissue=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2])),
                      fustat=trimws(sapply(as.character(metadata$characteristics_ch1.3),function(x) strsplit(x,":")[[1]][2])),
                      futime=trimws(sapply(as.character(metadata$characteristics_ch1.4),function(x) strsplit(x,":")[[1]][2])),
                      rfs_fustat = trimws(sapply(as.character(metadata$characteristics_ch1.5),function(x) strsplit(x,":")[[1]][2])),
                      rfs_futime = trimws(sapply(as.character(metadata$characteristics_ch1.6),function(x) strsplit(x,":")[[1]][2]))
                      )
clinical <- clinical[clinical$tissue=="primary hepatocellular carcinoma tumor",]
rownames(clinical) <- clinical[,1]
same <- intersect(row.names(clinical),colnames(rt1))
rt1 <- rt1[,same]
save(clinical,file = './00.data/LIHC/clinical/GSE76427.RData')
#转换ID#GPL3921	
library(illuminaHumanv4.db)
ids=toTable(illuminaHumanv4SYMBOL)
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
if(T){
geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)
rt1 <- read.table('./00.data/LIHC/megerData/GSE76427.txt',sep = '\t',header = T,check.names = F,row.names = 1)
same <- intersect(row.names(geneRT),row.names(rt1))
#write.table(rt1,file = './00.data/LIHC/GSE76427.txt',sep = '\t',quote = F)
}

if(T){
rt<-data.table::fread(file ='00.data/LIHC/TCGA-LIHC.htseq_fpkm.tsv/TCGA-LIHC.htseq_fpkm.tsv',data.table = F)
gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
rownames(rt) <- rt[,1]
rt <- rt[,-1]
data <- 2^rt-1
same <- intersect(row.names(gencode),row.names(data))
length(same)
data1 <- cbind(gencode[same,],data[same,])
data1 <- data1[,-c(1,3:6)]
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)
dir.create('./00.data/LIHC/megerData')
#FPKM转换为TPM
fpkmToTpm=function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(rt, 2, fpkmToTpm)
tpmOut=rbind(ID=colnames(tpm), tpm)
write.table(tpmOut, file="./00.data/LIHC/megerData/TCGA_LIHC.TPM.txt", sep="\t", col.names=F, quote=F)
}
if(T){
specimen_file <- read_delim(file="./00.data/LIHC/specimen.LIRI-JP.tsv/specimen.LIRI-JP.tsv",delim="\t",col_names = TRUE)
specimen_file$specimen_type=gsub("(.*?)\\-.*", "\\1", specimen_file$specimen_type)
table(specimen_file$specimen_type)
sepcimen_simplify_file <-   specimen_file %>% select(1,5,7) 
table(sepcimen_simplify_file$specimen_type)
sepcimen_simplify_file$specimen_type <- ifelse(sepcimen_simplify_file$specimen_type=='Normal ','N','T')

## 2.expression matrix 读入
exp_file <- read_delim(file = "./00.data/LIHC/exp_seq.LIRI-JP.tsv/exp_seq.LIRI-JP.tsv",delim = "\t",col_names = TRUE)

####3.将两文件合并，形成合并的表达矩阵
merge_df <- merge(exp_file,sepcimen_simplify_file,by="icgc_specimen_id",all.x=TRUE)

## gather 格式的表达数据
tmp <- merge_df[,c(1,23,24,8,9)] %>% tbl_df() %>% 
  unite(sample_id,icgc_specimen_id,icgc_donor_id.y,specimen_type,sep="-") 

## 查看捐赠者数目
length(unique(exp_file$icgc_donor_id))    # [1] 232
## 查看sample 数目
length(unique(exp_file$icgc_specimen_id)) # [1] 445

exp_matrix <- tmp %>% group_by(sample_id) %>% mutate(id=1:n()) %>% 
  spread(sample_id,normalized_read_count) %>% select(-"id")
}

ICGC <- data.table::fread('./00.data/LIHC/megerData/ICGC_LIHC.txt',data.table = F)
rownames(ICGC) <- ICGC[,1]
ICGC <- ICGC[,-1]
group <- grep('.*T$',colnames(ICGC))
tumor_icgc <- ICGC[,group]
outTab=rbind(geneNames=colnames(tumor_icgc), tumor_icgc)
write.table(outTab, file = './00.data/LIHC/megerData/ICGC_LIHC.txt', sep="\t", quote=F, col.names=F)


dir <- './00.data/LIHC/megerData/'
files=dir('./00.data/LIHC/megerData/') 
files=grep("txt$",files,value=T) 
files <- paste0(dir,files)

#获取交集基因
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit((strsplit(inputFile, "\\/|\\-")[[1]][5]), "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)


#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit((strsplit(inputFile, "\\/|\\-")[[1]][5]), "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对TCGA删除正常样品
  if(header[1] == "TCGA"){
    group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
    group=sapply(strsplit(group,""), "[", 1)
    rt=rt[,group==0]
    rt=t(rt)
    #row.names(rt)=substr(rownames(rt),1,20)
    rt=avereps(rt)
    rt=t(rt)
  }
  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] == "GSE76427"){
    rt=normalizeBetweenArrays(rt)
  }
  if(header[1] =="ICGC"){
    group <- grep('.*T$',colnames(rt))
    rt <- rt[,group]
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
      rt[rt<0]=0
      rt=log2(rt+1)}
    
  }
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="./00.data/LIHC/megerData/merge.txt", sep="\t", quote=F, col.names=F)
gene=read.table('./00.data/m6A_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(outTab))
geneExp=outTab[sameGene,]
save(geneExp,file = './00.data/LIHC/megerData/LIHC.RData')
