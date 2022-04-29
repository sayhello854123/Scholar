gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F)
#FPKM转换为TPM
fpkmToTpm=function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

COAD<-data.table::fread(file = '00.data/COAD/TCGA-COAD.htseq_fpkm.tsv/TCGA-COAD.htseq_fpkm.tsv',data.table = F)#513
rownames(COAD) <- COAD[,1]
COAD <- COAD[,-1]
COAD <- 2^COAD-1
same <- intersect(row.names(gencode),row.names(COAD))
length(same)
COAD <- cbind(gencode[same,],COAD[same,])
COAD <- COAD[,-c(1,3:6)]
rt=as.matrix(COAD)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
COAD=avereps(data)
#COAD=apply(COAD, 2, fpkmToTpm)
group_list=ifelse(as.numeric(substr(colnames(COAD),14,15)) < 10,'tumor','normal')
table(group_list)
COAD <- COAD[,group_list == "tumor"]
same1 <- intersect(as.vector(geneRT[,1]),row.names(COAD))
COAD <- COAD[same1,]
coad<- t(COAD)
type <- rep('COAD',times=length(rownames(coad)))
coad <- cbind(coad,type)
COAD <- coad

ESCA<-data.table::fread(file = '00.data/ESCA/TCGA-ESCA.htseq_fpkm.tsv/TCGA-ESCA.htseq_fpkm.tsv',data.table = F)#174
rownames(ESCA) <- ESCA[,1]
ESCA<- ESCA[,-1]
ESCA <- 2^ESCA-1
same <- intersect(row.names(gencode),row.names(ESCA))
length(same)
ESCA <- cbind(gencode[same,],ESCA[same,])
ESCA <- ESCA[,-c(1,3:6)]
rt=as.matrix(ESCA)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
ESCA=avereps(data)
#ESCA=apply(ESCA, 2, fpkmToTpm)
group_list=ifelse(as.numeric(substr(colnames(ESCA),14,15)) < 10,'tumor','normal')
table(group_list)
ESCA <- ESCA[,group_list == "tumor"]
same1 <- intersect(as.vector(geneRT[,1]),row.names(ESCA))
ESCA <- ESCA[same1,]
ESCA<- t(ESCA)
type <- rep('ESCA',times=length(rownames(ESCA)))
ESCA <- cbind(ESCA,type)

LIHC<-data.table::fread(file = '00.data/LIHC/TCGA-LIHC.htseq_fpkm.tsv/TCGA-LIHC.htseq_fpkm.tsv',data.table = F)#425
rownames(LIHC) <- LIHC[,1]
LIHC <- LIHC[,-1]
LIHC <- 2^LIHC-1
same <- intersect(row.names(gencode),row.names(LIHC))
length(same)
LIHC <- cbind(gencode[same,],LIHC[same,])
LIHC <- LIHC[,-c(1,3:6)]
rt=as.matrix(LIHC)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
LIHC=avereps(data)
#LIHC=apply(LIHC, 2, fpkmToTpm)
group_list=ifelse(as.numeric(substr(colnames(LIHC),14,15)) < 10,'tumor','normal')
table(group_list)
LIHC <- LIHC[,group_list == "tumor"]
same1 <- intersect(as.vector(geneRT[,1]),row.names(LIHC))
LIHC <- LIHC[same1,]
LIHC<- t(LIHC)
type <- rep('LIHC',times=length(rownames(LIHC)))
LIHC <- cbind(LIHC,type)



PAAD<-data.table::fread(file = '00.data/PAAD/TCGA-PAAD.htseq_fpkm.tsv/TCGA-PAAD.htseq_fpkm.tsv',data.table = F)
rownames(PAAD) <- PAAD[,1]
PAAD<- PAAD[,-1]
PAAD <- 2^PAAD-1
same <- intersect(row.names(gencode),row.names(PAAD))
length(same)
PAAD <- cbind(gencode[same,],PAAD[same,])
PAAD <- PAAD[,-c(1,3:6)]
rt=as.matrix(PAAD)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
PAAD=avereps(data)
#PAAD=apply(PAAD, 2, fpkmToTpm)
group_list=ifelse(as.numeric(substr(colnames(PAAD),14,15)) < 10,'tumor','normal')
table(group_list)
PAAD <- PAAD[,group_list == "tumor"]
same1 <- intersect(as.vector(geneRT[,1]),row.names(PAAD))
PAAD <- PAAD[same1,]
PAAD<- t(PAAD)
type <- rep('PAAD',times=length(rownames(PAAD)))
PAAD <- cbind(PAAD,type)



STAD<-data.table::fread(file = '00.data/STAD/TCGA-STAD.htseq_fpkm.tsv/TCGA-STAD.htseq_fpkm.tsv',data.table = F)
rownames(STAD) <- STAD[,1]
STAD <- STAD[,-1]
STAD <- 2^STAD-1
same <- intersect(row.names(gencode),row.names(STAD))
length(same)
STAD <- cbind(gencode[same,],STAD[same,])
STAD <- STAD[,-c(1,3:6)]
rt=as.matrix(STAD)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
STAD=avereps(data)
#STAD=apply(STAD, 2, fpkmToTpm)
group_list=ifelse(as.numeric(substr(colnames(STAD),14,15)) < 10,'tumor','normal')
table(group_list)
STAD <- STAD[,group_list == "tumor"]
same1 <- intersect(as.vector(geneRT[,1]),row.names(STAD))
STAD <- STAD[same1,]
STAD<- t(STAD)
type <- rep('STAD',times=length(rownames(STAD)))
STAD <- cbind(STAD,type)



READ<-data.table::fread(file = '00.data/READ/TCGA-READ.htseq_fpkm.tsv/TCGA-READ.htseq_fpkm.tsv',data.table = F)
rownames(READ) <- READ[,1]
READ<- READ[,-1]
READ <- 2^READ-1
same <- intersect(row.names(gencode),row.names(READ))
length(same)
READ <- cbind(gencode[same,],READ[same,])
READ <- READ[,-c(1,3:6)]
rt=as.matrix(READ)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
READ=avereps(data)
#READ=apply(READ, 2, fpkmToTpm)
group_list=ifelse(as.numeric(substr(colnames(READ),14,15)) < 10,'tumor','normal')
table(group_list)
READ <- READ[,group_list == "tumor"]
same1 <- intersect(as.vector(geneRT[,1]),row.names(READ))
READ <- READ[same1,]
READ<- t(READ)
type <- rep('READ',times=length(rownames(READ)))
READ <- cbind(READ,type)



library(ggplot2)
rt <- read.table('./00.data/1.txt',sep = '\t',header = T,check.names = F,row.names = 1 )
com1 <- prcomp(rt[,1:25], center = TRUE,scale. = TRUE)
df1<-com1$x
df1<-data.frame(df1,rt$type)
summ<-summary(com1)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
p2<-ggplot(data = df1,aes(x=PC1,y=PC2,color=rt$type))+
  stat_ellipse(aes(fill=rt$type),
               type = "norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)
p2+scale_fill_manual(values = c("#DC143C","#FF00FF","#8A2BE2",'#0000FF','#00FF7F','#FFA500'))+
  scale_colour_manual(values = c("#DC143C","#FF00FF","#8A2BE2",'#0000FF','#00FF7F','#FFA500'))

