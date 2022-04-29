library(GEOquery)
library(limma)
library(stringr)
library(AnnoProbe)
library(tinyarray)
rm(list = ls())
options(stringsAsFactors = F)
if(T){
Sys.setenv("VROOM_CONNECTION_SIZE" = 524288 * 2)
gset <- getGEO('GSE79668', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
a=gset[[1]]
metadata=pData(a)
clinical = data.frame(gsm=metadata[,2],
                      id = metadata[,1],
                      tissue=metadata[,8],
                      futime=trimws(sapply(as.character(metadata$characteristics_ch1.1),function(x) strsplit(x,":")[[1]][2])),
                      fustat=trimws(sapply(as.character(metadata$characteristics_ch1),function(x) strsplit(x,":")[[1]][2])),
                      Age = trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2])),
                      Gender = trimws(sapply(as.character(metadata$characteristics_ch1.4),function(x) strsplit(x,":")[[1]][2])),
                      T = trimws(sapply(as.character(metadata$characteristics_ch1.5),function(x) strsplit(x,":")[[1]][2])),
                      N = trimws(sapply(as.character(metadata$characteristics_ch1.6),function(x) strsplit(x,":")[[1]][2])),
                      M = trimws(sapply(as.character(metadata$characteristics_ch1.7),function(x) strsplit(x,":")[[1]][2]))
                      
)
rownames(clinical) <- clinical[,2]
rownames(clinical) <- substr(rownames(clinical),1,13)
outTab1=rbind(id_name=colnames(clinical), clinical)
write.table(outTab1, file="./00.data/PAAD/03.GSE79668/clinical.txt", sep="\t", quote=F, col.names=F)
clinical<- read.table('./00.data/PAAD/03.GSE79668/clinical.txt',sep = '\t',header = T,row.names = 1,check.names = F)
}
#转换TPM
if(T){
rt <- read.table('./00.data/PAAD/03.GSE79668/GSE79668_51_tumors_sharedgenecounts.txt',sep = '\t',header = T)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)
expMatrix <- rt
eff_length2 <-read.csv("./00.data/PAAD/03.GSE79668/eff_length.csv", row.names = 1, header = T)
gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
same <- intersect(row.names(gencode),row.names(eff_length2))
length(same)
eff_length2 <- cbind(eff_length2[same,],gencode[same,])
eff_length2$gene_id <- rownames(eff_length2)
feature_ids <- rownames(expMatrix)
# 检查gtf文件和表达量输入文件里基因名的一致性
if (! all(feature_ids %in% eff_length2[,3])){
  tbl <- table(feature_ids %in% eff_length2[,3])
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
  warning(msg1)
  
} 

if (! identical(feature_ids, eff_length2[,3])){
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}
# trim the expression matrix and effetive gene length
expMatrix <- expMatrix[feature_ids %in% eff_length2[,3],]
mm <- match(rownames(expMatrix), eff_length2[,3])
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), eff_length2[,3])){
  print("GTF and expression matix now have the same gene and gene in same order")
}
colnames(eff_length2)[1]='eff_length'
x <- expMatrix / eff_length2$eff_length
expMatrix_tpm <- t( t(x) / colSums(x) ) * 1e6 
#查看前三个基因的TPM值
expMatrix_tpm[1:3,]
}
if(T){
rt <- expMatrix_tpm
# rownames(rt) <- gsub("(.*?)\\_(.*?)", "\\2",rownames(rt))
rt1 <- t(rt)
same1 <- intersect(row.names(clinical),row.names(rt1))
rt2 <- cbind(clinical[same1,],rt1[same1,])
rownames(rt2) <- rt2[,1]
rt2 <- rt2[,-c(1:3)]
outTab2=rbind(id_name=colnames(rt2), rt2)
write.table(outTab2, file="./00.data/PAAD/03.GSE79668/merge_clinical.txt", sep="\t", quote=F, col.names=F)
rt <- t(rt2[,-c(1:7)])
outTab=rbind(geneNames=colnames(rt), rt)
write.table(outTab, file="./00.data/PAAD/03.GSE79668/GSE79668.txt", sep="\t", quote=F, col.names=F)
}