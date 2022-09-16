
library(impute)
library(limma)
library(ggplot2)
library(ggpubr)
if(T){
rt=read.table('./02.bulk/pRRophetic/drug.txt',sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
drug=rt[,2:ncol(rt)]
dimnames=list(rownames(drug),colnames(drug))
data=matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)

#对药物数据补缺
mat=impute.knn(data)
drug=mat$data
drug=avereps(drug)

#读取表达输入文件,并对输入文件整理
rt=read.table('./02.bulk/pRRophetic/geneExp.txt',sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
exp=avereps(data)

#提取特定基因表达
gene=read.table('./02.bulk/pRRophetic/gene.txt',sep="\t",header=F,check.names=F)
genelist=as.vector(gene[,1])
genelist=gsub(" ","",genelist)
genelist=intersect(genelist,row.names(exp))
exp=exp[genelist,]

#对基因循环
outTab=data.frame()
for(Gene in row.names(exp)){
  x=as.numeric(exp[Gene,])
  #对药物循环
  for(Drug in row.names(drug)){
    y=as.numeric(drug[Drug,])
    corT=cor.test(x,y,method="pearson")
    cor=corT$estimate
    pvalue=corT$p.value
    if(pvalue<0.05){
      outVector=cbind(Gene,Drug,cor,pvalue)
      outTab=rbind(outTab,outVector)
    }
  }
}

#输出相关性分析结果
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="./02.bulk/pRRophetic/drugCor.txt",sep="\t",row.names=F,quote=F)

#可视化
plotList=list()
for(i in rownames(outTab)){
  Gene=outTab[i,1]
  Drug=outTab[i,2]
  x=as.numeric(exp[Gene,])
  y=as.numeric(drug[Drug,])
  cor=sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=sprintf("%.03f",as.numeric(outTab[i,4]))
  if((abs(abs(as.numeric(outTab[i,3]))>0.5)&(as.numeric(outTab[i,4])<0.01))){
  df1=as.data.frame(cbind(x,y))
  if(as.numeric(outTab[i,4])<0.001){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  p1=ggplot(data = df1, aes(x = x, y = y))+
    geom_point(size=1)+
    stat_smooth(method="lm",se=FALSE, formula=y~x)+
    labs(x="",y="",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
    theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
    theme_bw()
  plotList[[i]]=p1
  }
}
}
#保存输出的图片
pdf(file="./02.bulk/pRRophetic/drugCor.pdf", width = 12,height =4)
ggarrange(plotlist=plotList,nrow=1,ncol=3)
dev.off()

