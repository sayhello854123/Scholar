
cnvplot <- function(x,y){
tmp <- read.table(x,sep = '\t',
                  header = T,row.names = 1,check.names = F)

gencode <- data.table::fread('00.data/gencode.v22.annotation.gene.probeMap',data.table = F)
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
geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F, row.names=1)
same1 <- intersect(row.names(geneRT),row.names(data))
rt <- data[same1,]

GAIN=rowSums(rt> 0)       #拷贝数增加的样品数目
LOSS=rowSums(rt< 0)       #拷贝数缺失的样品数目
GAIN=GAIN/ncol(rt)*100      #拷贝数增加的百分率
LOSS=LOSS/ncol(rt)*100      #拷贝数缺失的百分率
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#绘制图形
data.max = apply(data, 1, max)
pdf(file=y, width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col='#DC143C', cex=3)
points(bar,data[,"LOSS"], pch=20, col='#00BFFF', cex=3)
legend("top", legend=c('GAIN','LOSS'), col=c('#DC143C','#00BFFF'), pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()
}
##COAD##
x <- './00.data/COAD/TCGA-COAD.gistic.tsv/TCGA-COAD.gistic.tsv'
y <- './01.step1/02.CNVfreq/COAD_CNVfreq.pdf'
cnvplot(x,y)
##ESCA##
x <- './00.data/ESCA/TCGA-ESCA.gistic.tsv/TCGA-ESCA.gistic.tsv'
y <- './01.step1/02.CNVfreq/ESCA_CNVfreq.pdf'
cnvplot(x,y)
###LIHC###
x <- './00.data/LIHC/TCGA-LIHC.gistic.tsv/TCGA-LIHC.gistic.tsv'
y <- './01.step1/02.CNVfreq/LIHC_CNVfreq.pdf'
cnvplot(x,y)
###PAAD###
x <- './00.data/PAAD/TCGA-PAAD.gistic.tsv/TCGA-PAAD.gistic.tsv'
y <- './01.step1/02.CNVfreq/PAAD_CNVfreq.pdf'
cnvplot(x,y)
##STAD###
x <- './00.data/STAD/TCGA-STAD.gistic.tsv/TCGA-STAD.gistic.tsv'
y <- './01.step1/02.CNVfreq/STAD_CNVfreq.pdf'
cnvplot(x,y)
##READ###
x <- './00.data/READ/TCGA-READ.gistic.tsv/TCGA-READ.gistic.tsv'
y <- './01.step1/02.CNVfreq/READ_CNVfreq.pdf'
cnvplot(x,y)
