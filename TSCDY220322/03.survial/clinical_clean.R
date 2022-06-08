library(limma)
library(survival)
library(survminer)
library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")

##ICGC
if(T){
  ICGC_cli <- data.table::fread('./00.data/donor.LIRI-JP.tsv',data.table = F)
  cli <- ICGC_cli[,c(1,17,6,9,5)]
  cli$fustat <- ifelse(cli$donor_vital_status=='alive','0','1')
  cli$futime <- cli$donor_survival_time
  ICGC_cli <- cli[,c(1,7,6,4,5)] 
  colnames(ICGC_cli) <- c('id','futime','fustat','age','gender')
  rownames(ICGC_cli) <- ICGC_cli[,1]
  ICGC_cli$futime=ICGC_cli$futime/365
  load("./03.survial/LIHC.RData")
  colnames(geneExp)
  group1 <- grep('ICGC.*',colnames(geneExp))
  length(group1)
  ICGC <- t(geneExp[,group1])
  ICGC1 <- ICGC 
  id <- row.names(ICGC)
  ICGC1 <- cbind(id,ICGC1)
  rownames(ICGC1) <- sapply(strsplit(row.names(ICGC),"\\-"), "[", 2)
  ICGC1 <- data.frame(ICGC1)
  same <- intersect(row.names(ICGC_cli),row.names(ICGC1))
  ICGC2 <- cbind(ICGC_cli[same,],ICGC1[same,])
  rownames(ICGC2) <- ICGC2[,6]
  ICGC <- ICGC2[,-c(1,4,5,6)]
}

#GSE76427
if(T){
  load('./00.data/GSE76427.RData')
  load("./03.survial/LIHC.RData")
  group2 <- grep('GSE76427.*',colnames(geneExp))
  GSE76427 <- t(geneExp[,group2])
  GSE76427L <- data.frame(GSE76427) 
  id <- row.names(GSE76427)
  GSE76427L  <- cbind(id,GSE76427L)
  rownames(GSE76427L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE76427L))
  cli <- clinical[,c(4,3)]
  same <- intersect(row.names(cli),row.names(GSE76427L))
  GSE76427 <- cbind(cli[same,],GSE76427L[same,])
  rownames(GSE76427) <- GSE76427[,3]
  GSE76427 <- GSE76427[,-3]
}

###TCGA
if(T){
  tcga_cli <- data.table::fread('./00.data/TCGA-LIHC.survival.tsv',data.table = F)
  tcga_cli$futime <- tcga_cli$OS.time
  tcga_cli$fustat <- tcga_cli$OS
  cli <-tcga_cli[,c(1,5,6)] 
  rownames(cli) <- cli[,1]
  cli$futime=cli$futime/365
  load("./03.survial/LIHC.RData")
  colnames(geneExp)
  group3 <- grep('TCGA.*',colnames(geneExp))
  TCGA <- t(geneExp[,group3])
  TCGA1 <- TCGA
  id <- row.names(TCGA)
  TCGA1 <- cbind(id,TCGA1)
  rownames(TCGA1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(TCGA1))
  same <- intersect(row.names(cli),row.names(TCGA1))
  TCGA<- cbind(cli[same,],TCGA1[same,])
  rownames(TCGA) <- TCGA[,4]
  TCGA <- TCGA[,-c(1,4)]
}
LIHC_m5C <- rbind(GSE76427,ICGC,TCGA)
outTab=rbind(geneNames=colnames(LIHC_m5C),LIHC_m5C)
write.table(outTab, file="./03.survial/LIHC_m5C.txt", sep="\t", quote=F, col.names=F)
