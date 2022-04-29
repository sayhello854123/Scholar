library(limma)
library(survival)
library(survminer)
library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")
##ICGC
if(T){
ICGC_cli <- data.table::fread('./00.data/LIHC/clinical/donor.LIRI-JP.tsv',data.table = F)
cli <- ICGC_cli[,c(1,17,6,9,5)]
cli$fustat <- ifelse(cli$donor_vital_status=='alive','0','1')
cli$futime <- cli$donor_survival_time
ICGC_cli <- cli[,c(1,7,6,4,5)] 
colnames(ICGC_cli) <- c('id','futime','fustat','age','gender')
rownames(ICGC_cli) <- ICGC_cli[,1]
ICGC_cli$futime=ICGC_cli$futime/365
load("~/m6A/00.data/LIHC/megerData/LIHC.RData")
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
##GSE76427
if(T){
load('./00.data/LIHC/clinical/GSE76427.RData')
load("~/m6A/00.data/LIHC/megerData/LIHC.RData")
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
tcga_cli <- data.table::fread('./00.data/LIHC/clinical/TCGA-LIHC.survival.tsv',data.table = F)
tcga_cli$futime <- tcga_cli$OS.time
tcga_cli$fustat <- tcga_cli$OS
cli <-tcga_cli[,c(1,5,6)] 
rownames(cli) <- cli[,1]
cli$futime=cli$futime/365
load("~/m6A/00.data/LIHC/megerData/LIHC.RData")
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
LIHC_m6A <- rbind(GSE76427,ICGC,TCGA)
outTab=rbind(geneNames=colnames(LIHC_m6A),LIHC_m6A)
write.table(outTab, file="./00.data/LIHC/megerData/LIHC_m6A.txt", sep="\t", quote=F, col.names=F)
if(T){
#对基因进行循环，找出预后相关的基因
rt <- read.table('./00.data/LIHC/megerData/LIHC_m6A.txt', header=T, sep="\t", check.names=F,row.names = 1)
picDir <- './01.step1/05.m6A_survival/'
outTab=data.frame()
km=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  #km分析
  data=rt[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")
  #获取最优cutoff
  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
  #print(paste0(i, " ", res.cut$cutpoint[1]))
  #比较高低表达生存差异
  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  km=c(km, pValue)
  #对pvalue<0.05的基因绘制生存曲线
  if(pValue<0.05){
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    
    #绘制生存曲线
    surPlot=ggsurvplot(fit,
                       data=res.cat,
                       pval=pValue,
                       pval.size=6,
                       legend.title=i,
                       legend.labs=c("high","low"),
                       xlab="Time(years)",
                       ylab="Overall survival",
                       palette=c("red", "blue"),
                       break.time.by=1,
                       conf.int=T,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    filename <- paste0(i, ".pdf",sep='')
    outfile <- paste0(picDir,filename,sep="")
    pdf(file=outfile,onefile = FALSE,
        width = 6,         #图片的宽度
        height =5)         #图片的高度
    print(surPlot)
    dev.off()
  }
}
outTab=cbind(outTab, km)
write.table(outTab,file="./00.data/LIHC/survival_data/uniCox.txt",sep="\t",row.names=F,quote=F)
}
if(T){
load('./00.data/LIHC/megerData/LIHC.RData')
outTab=rbind(geneNames=colnames(geneExp), geneExp)
write.table(outTab, file="./00.data/LIHC/megerData/m6Amerge.txt", sep="\t", quote=F, col.names=F)
gene.exp  <-  read.table("./00.data/LIHC/megerData/m6Amerge.txt",header=T,sep="\t",check.names = F,row.names = 1)
gene.group <- read.table("./00.data/m6A_gene.txt",header=T,sep="\t")
gene.cox <- read.table("./00.data/LIHC/survival_data/uniCox.txt",header=T,sep="\t")

network_plot <- function(nodefile,edgefile,hua){
#基因取交集
colnames(gene.group) <- c('id','group')
genelist <- intersect(gene.group$id, gene.cox$id)
genelist <- intersect(genelist, rownames(gene.exp))
gene.group <- gene.group[match(genelist,gene.group$id),]
gene.group <- gene.group[order(gene.group$group),]
gene.exp <- gene.exp[match(gene.group$id,rownames(gene.exp)),]
gene.cox <- gene.cox[match(gene.group$id,gene.cox$id),]

#准备网络文件
gene.cor <- corr.test(t(gene.exp))
gene.cor.cor <- gene.cor$r
gene.cor.pvalue <- gene.cor$p
gene.cor.cor[upper.tri(gene.cor.cor)] = NA
gene.cor.pvalue[upper.tri(gene.cor.pvalue)] = NA
gene.cor.cor.melt <- melt(gene.cor.cor)   #gene1 \t gene2 \t cor
gene.cor.pvalue.melt <- melt(gene.cor.pvalue)
gene.melt <- data.frame(from = gene.cor.cor.melt$Var2,to=gene.cor.cor.melt$Var1,cor=gene.cor.cor.melt$value,pvalue=gene.cor.pvalue.melt$value)
gene.melt <- gene.melt[gene.melt$from!=gene.melt$to&!is.na(gene.melt$pvalue),,drop=F]
gene.edge <- gene.melt[gene.melt$pvalue<0.0001,,drop=F]
gene.edge$color <- ifelse(gene.edge$cor>0,'pink','#6495ED')
gene.edge$weight <- abs(gene.edge$cor)*6

#准备结果属性文件
gene.node <- gene.group
group.color <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(gene.node$group)))
gene.node$color <- group.color[as.numeric(as.factor(gene.node$group))]
gene.node$shape <- "circle"
gene.node$frame <- ifelse(gene.cox$HR>1,'purple',"green")
gene.node$pvalue <- gene.cox$pvalue
# pvalue size
pvalue.breaks <- c(0,0.0001,0.001,0.01,0.05,1)
pvalue.size <- c(16,14,12,10,8)
cutpvalue <- cut(gene.node$pvalue,breaks=pvalue.breaks)
gene.node$size <- pvalue.size[as.numeric(cutpvalue)]
nodefile <- nodefile
edgefile <- edgefile
write.table(gene.node, nodefile, sep="\t", col.names=T, row.names=F, quote=F)
write.table(gene.edge, edgefile, sep="\t", col.names=T, row.names=F, quote=F)


#绘制网络图
node = read.table(nodefile, header=T, sep="\t", comment.char="")
edge = read.table(edgefile, header=T, sep="\t", comment.char="")

g = graph.data.frame(edge,directed = FALSE)
node = node[match(names(components(g)$membership),node$id),]

if(!is.na(match('color',colnames(node)))) V(g)$color = node$color
if(!is.na(match('size',colnames(node)))) V(g)$size = node$size
if(!is.na(match('shape',colnames(node)))) V(g)$shape = node$shape
if(!is.na(match('frame',colnames(node)))) V(g)$frame = node$frame

#输出文件
pdf(file=hua, width=10, height=8)
par(mar=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nc=2),height=c(4,4,2),width=c(8,3))

#节点坐标 
coord = layout_in_circle(g)
degree.x = acos(coord[,1])
degree.y = asin(coord[,2])
degree.alpha = c()
for(i in 1:length(degree.x)){
  if(degree.y[i]<0) degree.alpha=c(degree.alpha,2*pi-degree.x[i]) else degree.alpha=c(degree.alpha,degree.x[i])
}
degree.cut.group = (0:8)/4*pi
degree.cut.group[1] = -0.0001
degree.cut = cut(degree.alpha,degree.cut.group)
degree.degree = c(-pi/4,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/4)
degree = degree.degree[as.numeric(degree.cut)]

#定义饼图,左半圆颜色代表m6A的类型，右半圆代表基因是高风险基因,还是低风险基因
values <- lapply(node$id,function(x)c(1,1))
V(g)$pie.color = lapply(1:nrow(node),function(x)c(node$color[x],node$frame[x]))
V(g)$frame = NA 

#绘制图形
plot(g,layout=layout_in_circle,vertex.shape="pie",vertex.pie=values,
     vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color=V(g)$color,vertex.frame.color=V(g)$frame,edge.color=E(g)$color,
     vertex.label.cex=2,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.4,
     vertex.color=V(g)$color,vertex.label.dist=1,vertex.label.degree=degree)
# label.degree : zero means to the right; and pi means to the left; up is -pi/2 and down is pi/2;  The default value is -pi/4
# label.dist If it is 0 then the label is centered on the vertex; If it is 1 then the label is displayed beside the vertex.

#绘制节点属性图例(m6A的类型)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
groupinfo = unique(data.frame(group=node$group,color=node$color))
legend("left",legend=groupinfo$group,col=groupinfo$color,pch=16,bty="n",cex=3)
#绘制风险图例(哪些基因是高风险的基因,哪些基因是低风险的基因)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
legend("left",legend=c('Risk factors','Favorable factors'),col=c('purple','green'),pch=16,bty="n",cex=2.5)
#绘制预后pvalue图例
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",axes=F,ylab="")
legend("top",legend=c('Postive correlation with P<0.0001','Negative correlation with P<0.0001'),lty=1,lwd=4,col=c('pink','#6495ED'),bty="n",cex=2.2)
legend('bottom',legend=c(0.0001,0.001,0.01,0.05,1),pch=16,pt.cex=c(1.6,1.4,1.2,1,0.8)*6,bty="n",ncol=5,cex=2.2,col="black",title="Cox test, pvalue")
dev.off()
}
hua <- "./01.step1/05.m6A_survival/network.pdf"
nodefile <- "./00.data/LIHC/survival_data/network.node.txt"
edgefile <- "./00.data/LIHC/survival_data/network.edge.txt"
network_plot(nodefile,edgefile,hua)
}