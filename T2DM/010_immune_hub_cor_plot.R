rm(list = ls())  
options(stringsAsFactors = F)#清空环境
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)
dir.create('13.immunCell_hub_cor')
gene="T cells CD4 memory resting" 

load('./01.data/09.immune_cor_input.RData')
if(T){
  x=as.numeric(immune1[,1])
  p <- list()
  for(i in colnames(immune1[,2:ncol(immune1)])) {
    y=as.numeric(immune1[,i])
    if(sd(y)>0.001){
      df1=as.data.frame(cbind(x,y))
      corT=cor.test(x,y,method="spearman")
      cor=corT$estimate
      pValue=corT$p.value}
    p <- c(p, 
           list(ggplot(df1, aes(x, y)) + labs(title = paste(i, sep = "")) +
                  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
                  ylab(paste(i, sep = ""))+xlab(gene)+
                  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
                  stat_cor(method = 'spearman', aes(x =x, y =y))))
    }
     
  ###
  x=as.numeric(immune2[,1])
  #k <- list()
  for(i in colnames(immune2[,2:ncol(immune2)])) {
    y=as.numeric(immune2[,i])
    if(sd(y)>0.001){
      df1=as.data.frame(cbind(x,y))
      corT=cor.test(x,y,method="spearman")
      cor=corT$estimate
      pValue=corT$p.value
    }
    p <- c(p, 
           list(ggplot(df1, aes(x, y)) + labs(title = paste(i, sep = "")) +
                  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
                  ylab(paste(i, sep = ""))+xlab(gene)+
                  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
                  stat_cor(method = 'spearman', aes(x =x, y =y))))
  }
}

Rmisc::multiplot(plotlist = p, layout = matrix(1:6, nrow = 2))

