library(ggpubr)                    
tciaFile="./02.bulk/IPS/tcia.txt"                
scoreFile="./02.bulk/IPS/tcgaRisk.txt"   


#读取免疫治疗打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#读取打分分组文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(ips), row.names(score))
ips=ips[sameSample, , drop=F]
score=score[sameSample,, drop=F]
data=cbind(ips, score)

#设置比较组
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
dir <- './02.bulk/IPS'
#绘制
for(i in colnames(data)[1:4]){
  rt=data[,c(i, "risk")]
  colnames(rt)=c("IPS", "group")
  gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
               xlab="TCGA-LIHC", ylab=i,
               legend.title="TCGA-LIHC",
               palette=c("#0066FF", "#FF0000"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  filename <- paste0(i, ".pdf",sep='')
  outfile <- paste(dir,filename,sep="/")
  pdf(outfile, width=6, height=5)
  print(gg1)
  dev.off()
}
