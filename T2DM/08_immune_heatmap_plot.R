rm(list = ls())  
options(stringsAsFactors = F)#清空环境
setwd('../')
dir.create('11.immune_Pheatmap')
library(pheatmap)
library(corrplot)

pFilter=0.05
load('./01.data/07.immune_heatmap.RData')
data=rt[rt[,"P-value"]<pFilter,]#删除p<0.05的细胞
data1=t(data[,1:(ncol(rt)-4)])
data2=data[,1:(ncol(rt)-4)]
#绘制热图
if(T){
    ##热图
    Type=c(rep("non−diabetic condition",10),rep("diabetic condition_donor",10))    #修改对照和处理组样品数目
    names(Type)=colnames(data1)
    Type=as.data.frame(Type)
  
    pdf("./11.immune_Pheatmap/immune_heatmap.pdf",height=5,width=10)
    pheatmap(data1, 
           annotation=Type, 
           color = colorRampPalette(c("green", "black", "red"))(50),
           cluster_cols =F,
           fontsize = 8,
           fontsize_row=7,
           fontsize_col=8)
     dev.off()
     ##柱状类图
     col=rainbow(nrow(data1),s=0.7,v=0.7)
     outpdf="./11.immune_Pheatmap/immune_barplot.pdf"
     pdf(outpdf,height=8,width=18)
     par(las=1,mar=c(8,5,4,16))
     a1 = barplot(data1,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
     a2=axis(2,tick=F,labels=F)
     axis(2,a2,paste0(a2*100,"%"))
     axis(1,a1,labels=F)
     par(srt=60,xpd=T);text(a1,-0.02,colnames(data1),adj=1,cex=1);par(srt=0)
     ytick2 = cumsum(data1[,ncol(data1)])
     ytick1 = c(0,ytick2[-length(ytick2)])
     legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data1),col=col,pch=15,bty="n",cex=1.3)
     dev.off()
     ###相关性热图
     pdf("./11.immune_Pheatmap/immune_corHeatmap.pdf",height=13,width=13)              #保存图片的文件名称
     corrplot(corr=cor(data2),
              method = "color",
              order = "hclust",
              tl.col="black",
              addCoef.col = "black",
              number.cex = 0.8,
              col=colorRampPalette(c("blue", "white", "red"))(50),
     )
     dev.off()
}




