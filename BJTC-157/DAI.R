scRNA <- readRDS("C:/Users/16385/Desktop/BJTC-157/00.data/01.single cell/scRNA2.rds")
library(limma)
library(ser)
a <- as.data.frame(scRNA@assays$RNA@counts)
colnames(a) <-scRNA@meta.data$Type
table(colnames(a))

B <- a[,colnames(a)=="B cell"]
Tcell <- a[,colnames(a)=="T cell"]
TAM <- a[,colnames(a)=="TAM"]
rt <- cbind(B,Tcell,TAM)
write.table(rt,file = './00.data/01.single cell/cibersoft_input.txt',sep = '\t',quote = F)

#if(T){
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("./00.data/01.single cell/Homo_sapiens.GRCh38.104.chr.gtf/Homo_sapiens.GRCh38.104.chr.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

exons_gene_2_lens=data.frame(t(data.frame(exons_gene_lens)))

names(exons_gene_2_lens)="length"

head(exons_gene_2_lens)[1]

write.csv(exons_gene_2_lens,"./00.data/01.single cell/Homo_sapiens.GRCh38.104.chr.gtf/eff2_length.csv",row.names = T)
}

eff_length2 <-read.csv("./00.data/01.single cell/Homo_sapiens.GRCh38.104.chr.gtf/eff2_length.csv", row.names = 1, header = T)
eff_length2$gene_id <- rownames(eff_length2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
id = bitr(eff_length2$gene_id,
          fromType = "ENSEMBL",
          toType = "SYMBOL",
          OrgDb = "org.Hs.eg.db")

group <- eff_length2$gene_id%in%id$ENSEMBL
table(group)
eff_length <- eff_length2[group==TRUE,]
id1 <- id[!duplicated(id$ENSEMBL),]
rownames(id1) <- id1[,1]
same <- intersect(row.names(id1),row.names(eff_length))
genelenth <- cbind(id1[same,],eff_length[same,])

expMatrix <- read.table('./00.data/01.single cell/cibersoft_input.txt',sep = '\t',header = T)
rt=as.matrix(expMatrix)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
expMatrix=data[rowMeans(data)>0,]

# 从输入数据里提取基因名
feature_ids <- rownames(expMatrix)


# 检查gtf文件和表达量输入文件里基因名的一致性
if (! all(feature_ids %in% rownames(eff_length2))){
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
  warning(msg1)
} 

if (! identical(feature_ids, rownames(eff_length2))){
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

# trim the expression matrix and effetive gene length
expMatrix <- expMatrix[feature_ids %in% genelenth$SYMBOL,]
mm <- match(rownames(expMatrix), genelenth$SYMBOL)
eff_length2 <- genelenth[mm, ]


#identical : The safe and reliable way to test two objects for being exactly  equal. 
#It returns TRUE in this case, FALSE in every other case.--------

if (identical(rownames(eff_length2), rownames(expMatrix))){
  print("GTF and expression matix now have the same gene and gene in same order")
}


countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}


fpkms <- apply(expMatrix, 2, countToFpkm, effLen = eff_length2$length)

fpkms.m<-data.frame(fpkms)
colnames(fpkms.m)<-colnames(expMatrix)
dim(fpkms.m)

write.csv(fpkms.m,"./00.data/01.single cell/fpkm.csv",row.names = T)
