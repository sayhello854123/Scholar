rm(list = ls())
options(stringsAsFactors = F) 
gc()
library(GEOquery)
library(limma)
library(tidyverse)
library(sva)

rt <- read.table('./03.survial/merge.txt',sep = '\t',header = T,row.names = 1,check.names = F)
gene=read.table('./00.data/m5c_gene.txt', header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(rt))
geneExp=rt[sameGene,]
save(geneExp,file = './03.survial/LIHC.RData')
