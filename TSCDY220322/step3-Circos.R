rm(list = ls())
options(stringsAsFactors = F) 
gc()
library("RCircos")
if(T){
  gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
  geneRT=read.table("./00.data/m5c_gene.txt", header=T, sep="\t", check.names=F, row.names=1)
  same <- intersect(row.names(geneRT),gencode[,2])
  data <- gencode[same,]
  group <- gencode[,2]%in%rownames(geneRT)
  data <- gencode[group==TRUE,]
  geneLabel<- data[,c(3:5)]
  geneLabel$Gene <- data[,2]
  cytoBandIdeogram=read.table("./00.data/refer.txt", header=T, sep="\t")
}
Circosplot <- function(x,y){
  tmp <- read.table(x,sep = '\t',
                    header = T,row.names = 1,check.names = F)
  #gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
  rownames(gencode) <- gencode[,1]
  same <- intersect(row.names(gencode),row.names(tmp))
  length(same)
  data1 <- cbind(gencode[same,],tmp[same,])
  data1 <- data1[,-c(1,3:6)]
  rt=as.matrix(data1)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  #geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F, row.names=1)
  same1 <- intersect(row.names(geneRT),row.names(data))
  rt <- data[same1,]
  GAIN=rowSums(rt> 0)
  LOSS=rowSums(rt< 0)
  group1 <- ifelse(GAIN>LOSS,'1','-1')
  data2 <- cbind(row.names(geneRT),group1)
  colnames(data2) <- c('Gene','seg.mean')
  scatter<- geneLabel
  scatter1 <-merge(scatter,data2,by ='Gene')
  scatter <- scatter1[,-1]
  #初始化圈图
  chr.exclude <- NULL
  cyto.info <- cytoBandIdeogram
  tracks.inside <- 5
  tracks.outside <- 0
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
  #设置圈图参数
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$text.size=1
  rcircos.params$point.size=5
  RCircos.Reset.Plot.Parameters(rcircos.params)
  #输出文件
  pdf(file=y, width=8, height=8)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  #散点图
  RCircos.Scatter.Data=scatter
  data.col <- 4
  track.num <- 1
  side <- "in"
  RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)
  #加上基因名称
  RCircos.Gene.Label.Data=geneLabel
  name.col <- 4
  side <- "in"
  track.num <- 2
  RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
  track.num <- 3
  RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
  dev.off()
}
###LIHC###
x <- './00.data/TCGA-LIHC.gistic.tsv'
y <- './01.maftools/LIHC_Circos.pdf'
Circosplot(x,y)
