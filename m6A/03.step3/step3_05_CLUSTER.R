library(ConsensusClusterPlus) 
library(pheatmap)
library(limma)
library(GSEABase)
library(GSVA)
library(survival)
library(survminer)
options(stringsAsFactors = F)

data=read.table("./00.data/ESCA/04.meger/m6A_merge.txt",header=T,sep="\t",check.names = F,row.names = 1)
data=as.matrix(data)

  #聚类
  maxK=9
  results=ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='./03.step3/03.cluster_heatmap/',
                               clusterAlg="pam",
                               distance="euclidean",
                               seed=123456,
                              plot="png")
  clusterNum=3     #分几类，根据判断标准判断
  if(T){
  
  cluster=results[[clusterNum]][["consensusClass"]]
  cluster=as.data.frame(cluster)
  colnames(cluster)=c("cluster")
  letter=c("A","B","C","D","E","F","G")
  uniqClu=levels(factor(cluster$cluster))
  cluster$cluster=letter[match(cluster$cluster, uniqClu)]
  clusterOut=rbind(ID=colnames(cluster), cluster)
  write.table(clusterOut, file="./00.data/ESCA/06.cluster/m6aCluster1.txt", sep="\t", quote=F, col.names=F)
  
  cluster=read.table("./00.data/ESCA/06.cluster/m6aCluster.txt", header=T, sep="\t", check.names=F, row.names=1)
  row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
  cli=read.table("./00.data/ESCA/04.meger/ESCA_m6A.txt", header=T, sep="\t", check.names=F, row.names=1)
  cli <- cli[,c(1:2)]
  sameSample=intersect(row.names(cluster), row.names(cli))
  rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
  
  length=length(levels(factor(rt$cluster)))
  diff=survdiff(Surv(futime, fustat) ~ cluster, data = rt)
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="cluster",
                     legend.labs=levels(factor(rt[,"cluster"])),
                     legend = c(0.9, 0.8),
                     font.legend=10,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.28)
  surPlot
  pdf(file="./03.step3/03.cluster_heatmap/survival.pdf",onefile = FALSE,width=7,height=6.5)
  print(surPlot)
  dev.off()
  }
  if(T){
    exp=read.table("./00.data/ESCA/04.meger/m6A_merge.txt", header=T, sep="\t", check.names=F, row.names=1)
    exp=t(exp)
    cluster=read.table('./00.data/ESCA/06.cluster/m6aCluster.txt', header=T, sep="\t", check.names=F, row.names=1)
    row.names(cluster) <- gsub('\\.', '-',row.names(cluster))
    sameSample=intersect(row.names(exp), row.names(cluster))
    exp=exp[sameSample, , drop=F]
    cluster=cluster[sameSample, , drop=F]
    expCluster=cbind(exp, cluster)
    Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
    #rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
    expCluster=cbind(expCluster, Project)
    cli=read.csv("./00.data/ESCA/04.meger/ALL_cli.csv", row.names = 1, header = T)
    cli$Fustat <- ifelse(cli$fustat=='0','Alive','Dead')
    cli <- cli[,-c(1,2)]
    sameSample=intersect(row.names(expCluster), row.names(cli))
    expCluster=expCluster[sameSample,,drop=F]
    cli=cli[sameSample,,drop=F]
    data=cbind(expCluster, cli)
    #data <- data[,-c(28:30)]
    data=data[order(data$cluster),]
    Type=data[,((ncol(exp)+1):ncol(data))]
    data=t(data[,1:ncol(exp)])
    
    bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    ann_colors=list()
    m6aCluCol=bioCol[1:length(levels(factor(Type$cluster)))]
    names(m6aCluCol)=levels(factor(Type$cluster))
    ann_colors[["cluster"]]=m6aCluCol
    
    plot1=pheatmap(data,
                   annotation=Type,
                   annotation_colors = ann_colors,
                   color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
                   cluster_cols =F,
                   cluster_rows =F,
                   scale="row",
                   show_colnames=F,
                   fontsize=6,
                   fontsize_row=6,
                   fontsize_col=6)
    
    ggsave("./03.step3/03.cluster_heatmap/heatmap.jpg",plot =plot1, height=7, width=9,dpi = 1000)
    ggsave("./03.step3/03.cluster_heatmap/heatmap.pdf",plot =plot1, height=7, width=9,dpi = 1000)
  }
  