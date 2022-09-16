rm(list = ls())  
options(stringsAsFactors = F)#清空环境

library(pROC)             
library(stringr)      
dir.create('08.Hub_Roc') 

load('./01.data/05.ROC_input.RData')
setwd("./08.Hub_Roc")
y=colnames(rt)[1]

for(i in colnames(rt[,2:ncol(rt)])){
  #绘制
  rocobj1=roc(rt[,y], as.vector(rt[,i]))
  ROC_Value <- as.numeric(str_split(rocobj1$auc, " ")[[1]]) 
  if(ROC_Value>0.6){
    pdf(file=paste0(i,".pdf"),width=6,height=5)
    plot(rocobj1, print.auc=TRUE, col="black",auc.polygon=TRUE,auc.polygon.col="gray",
         print.auc.x=0.4,print.auc.y=0.4,print.thres.cex=0.6,print.auc.col = "black",
         print.thres = T,print.thres.pch = 12, print.thres.adj = 1.1, print.thres.col = "red", 
         lty=1,main=i,mfrow=c(1,1))
    dev.off()}
}
