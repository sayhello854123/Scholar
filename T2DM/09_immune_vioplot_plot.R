rm(list = ls())  
options(stringsAsFactors = F)#清空环境
library(vioplot) 
dir.create('12.immune_vioplot/')

load('./01.data/08.immune_vioplot_input.RData')

if(T){
    pdf("./12.immune_vioplot/vioplot.pdf",height=8,width=15)#保存图片的文件名称
    normal=10 #样本数目                                                           
    tumor=10
    par(las=1,mar=c(10,6,3,3))#上下坐下间距
    x=c(1:ncol(rt))
    y=c(1:ncol(rt))
    plot(x,y,
         xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
         main="",xlab="", ylab="Fraction",
         pch=21,
         col="white",
         xaxt="n")
    
    #对每个免疫细胞循环，绘制vioplot，正常用蓝色表示，不正常用红色表示
    for(i in 1:ncol(rt)){
      normalData=rt[1:normal,i]
      tumorData=rt[(normal+1):(normal+tumor),i]
      vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
      vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
      wilcoxTest=wilcox.test(normalData,tumorData)# wilcoxTest检验
      p=round(wilcoxTest$p.value,3)
      mx=max(c(normalData,tumorData))
      lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
      text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
      text(seq(1,64,3),-0.02,xpd = NA,labels=colnames(rt),cex = 0.8,srt = 40,pos=2)
    }
      dev.off()
}
