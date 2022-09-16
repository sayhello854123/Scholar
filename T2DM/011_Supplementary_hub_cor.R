library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(patchwork)

pFilter=0.05 
gene="T cells CD4 memory resting" 

load('./01.data/09.immune_cor_input.RData') 
x=as.numeric(immune1[,1])
y=as.numeric(immune1[,2])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
p1=ggplot(df1, aes(x, y)) + ggtitle("STY1") +
  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
  ylab(colnames(immune1)[2])+xlab(gene)+
  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))

y=as.numeric(immune1[,3])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
p2=ggplot(df1, aes(x, y)) + ggtitle("SPP1") +
  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
  ylab(colnames(immune1)[3])+xlab(gene)+
  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
###
y=as.numeric(immune1[,4])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
p3=ggplot(df1, aes(x, y)) + ggtitle("CDK5") +
  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
  ylab(colnames(immune1)[4])+xlab(gene)+
  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))

x=as.numeric(immune2[,1])
y=as.numeric(immune2[,2])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
p4=ggplot(df1, aes(x, y)) + ggtitle("STY1") +
  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
  ylab(colnames(immune2)[2])+xlab(gene)+
  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))

y=as.numeric(immune2[,3])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
p5=ggplot(df1, aes(x, y)) + ggtitle("SPP1") +
  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
  ylab(colnames(immune2)[3])+xlab(gene)+
  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
###
y=as.numeric(immune2[,4])
if(sd(y)>0.001){
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="spearman")
  cor=corT$estimate
  pValue=corT$p.value}
p6=ggplot(df1, aes(x, y)) + ggtitle("CDK5") +
  theme(plot.title = element_text(size=28,  family="serif",face = "bold",hjust = 0.5)) +
  ylab(colnames(immune2)[4])+xlab(gene)+
  geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))

P <- (p1|p2|p3)/(p4|p5|p6)
ggplot2::ggsave(filename = './13.immunCell_hub_cor/cor.pdf',plot = P, width =8, height =6) 

