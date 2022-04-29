library(GEOquery)
library(limma)
library(tidyverse)
library(sva)
library(AnnoProbe)
library(tidyverse)

rt<-data.table::fread(file ='00.data/COAD/01.TCGA/TCGA-COAD.htseq_fpkm.tsv/TCGA-COAD.htseq_fpkm.tsv',data.table = F)
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
#FPKM转换为TPM
fpkmToTpm=function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(rt, 2, fpkmToTpm)
tpmOut=rbind(ID=colnames(tpm), tpm)
write.table(tpmOut, file="./00.data/COAD/01.TCGA/TCGA.TPM.txt", sep="\t", col.names=F, quote=F)

##生存数据的整理
cli <- data.table::fread(file ='./00.data/COAD/01.TCGA/TCGA-COAD.GDC_phenotype.tsv',data.table = F)
rownames(cli) <- cli[,1]
cli1 <- data.table::fread(file ='00.data/COAD/01.TCGA/TCGA-COAD.survival.tsv',data.table = F)
rownames(cli1) <- cli1[,1]
same <- intersect(row.names(cli1),row.names(cli))
cli3 <- cbind(cli1[same,],cli[same,])
colnames(cli3)
clinical <- cli3[,c(4,2,6,77,97,96,45,44,43,79)]
cliOut=rbind(ID=colnames(clinical), clinical)
write.table(cliOut, file="./00.data/COAD/01.TCGA/TCGA_cli.txt", sep="\t", col.names=F, quote=F)
