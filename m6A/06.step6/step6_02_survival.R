library(limma)
library(survival)
library(survminer)
library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")
##GSE15459
if(T){
  geneExp <- read.table('./00.data/STAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  group2 <- grep('GSE15459.*',colnames(geneExp))
  GSE15459 <- t(geneExp[,group2])
  GSE15459L <- data.frame(GSE15459) 
  id <- row.names(GSE15459)
  GSE15459L  <- cbind(id,GSE15459L)
  rownames(GSE15459L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE15459L))
  cli <- read.csv('./00.data/STAD/03.GEO/GSE15459_cli.csv',header = T,row.names = 1)
  #cli$fustat <- ifelse(cli$fustat=='no death','0','1')
  clinical <- cli
  clinical$futime=clinical$futime/12
  same <- intersect(row.names(clinical),row.names(GSE15459L))
  GSE15459 <- cbind(clinical[same,],GSE15459L[same,])
  rownames(GSE15459) <- GSE15459[,6]
  outTab1=rbind(geneNames=colnames(GSE15459), GSE15459)
  write.table(outTab1, file = './00.data/STAD/04.meger/GSE15459_ALL.txt', sep="\t", quote=F, col.names=F)
  GSE15459 <- GSE15459[,-c(3:6)]
}

##GSE34942
if(T){
  geneExp <- read.table('./00.data/STAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  group3 <- grep('GSE34942.*',colnames(geneExp))
  GSE34942 <- t(geneExp[,group3])
  GSE34942L <- data.frame(GSE34942) 
  id <- row.names(GSE34942)
  GSE34942L  <- cbind(id,GSE34942L)
  rownames(GSE34942L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE34942L))
  cli <- read.csv('./00.data/STAD/03.GEO/GSE39942_cli.csv',header = T,row.names = 1)
  clinical <- cli
  clinical$futime=clinical$futime/12
  same <- intersect(row.names(clinical),row.names(GSE34942L))
  GSE34942 <- cbind(clinical[same,],GSE34942L[same,])
  rownames(GSE34942) <- GSE34942[,6]
  outTab1=rbind(geneNames=colnames(GSE34942), GSE34942)
  write.table(outTab1, file = './00.data/STAD/04.meger/GSE34942_ALL.txt', sep="\t", quote=F, col.names=F)
  GSE34942 <- GSE34942[,-c(3:6)]
}
#TCGA
if(T){
  tcga_cli <- read.csv('./00.data/STAD/01.TCGA/TCGA_cli.csv',header = T,row.names = 1)
  cli<- tcga_cli
  cli$futime=cli$futime/365
  #colnames(cli) <- c('futime','fustat','age','gender','stage','grade','T','N','M')
  geneExp <- read.table('./00.data/STAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  colnames(geneExp)
  group3 <- grep('TCGA.*',colnames(geneExp))
  TCGA <- t(geneExp[,group3])
  TCGA1 <- TCGA
  id <- row.names(TCGA)
  TCGA1 <- cbind(id,TCGA1)
  rownames(TCGA1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(TCGA1))
  same <- intersect(row.names(cli),row.names(TCGA1))
  TCGA<- cbind(cli[same,],TCGA1[same,])
  rownames(TCGA) <- TCGA[,10]
  outTab2=rbind(geneNames=colnames(TCGA), TCGA)
  write.table(outTab2, file = './00.data/STAD/04.meger/TCGA_ALL.txt', sep="\t", quote=F, col.names=F)
  TCGA <- TCGA[,-c(3:10)]
}
PAAD_m6A <- rbind(GSE15459,GSE34942,TCGA)
outTab3=rbind(IDs=colnames(PAAD_m6A),PAAD_m6A)
write.table(outTab3, file="./00.data/STAD/04.meger/STAD_m6A.txt", sep="\t", quote=F, col.names=F)

if(T){
  #对基因进行循环，找出预后相关的基因
  rt <- read.table('./00.data/STAD/04.meger/STAD_m6A.txt', header=T, sep="\t", check.names=F,row.names = 1)
  picDir <- './06.step6/01.survial_m6A/'
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
  write.table(outTab,file="./00.data/STAD/05.survival/uniCox.txt",sep="\t",row.names=F,quote=F)
}