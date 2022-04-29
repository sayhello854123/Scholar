library(limma)
library(survival)
library(survminer)
library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")
##GSE79668
if(T){
  geneExp <- read.table('./00.data/PAAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  group2 <- grep('GSE79668.*',colnames(geneExp))
  GSE79668 <- t(geneExp[,group2])
  GSE79668L <- data.frame(GSE79668) 
  id <- row.names(GSE79668)
  GSE79668L  <- cbind(id,GSE79668L)
  rownames(GSE79668L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE79668L))
  cli <- read.table('./00.data/PAAD/03.GSE79668/clinical.txt',header = T,sep = '\t',,check.names = F,row.names = 1)
  rownames(cli) <- cli[,1]
  table(cli$fustat)
  cli$fustat <- ifelse(cli$fustat=='Alive','0','1')
  clinical <- cli[,c(1,4:10)]
  clinical$futime=clinical$futime/365
  same <- intersect(row.names(clinical),row.names(GSE79668L))
  GSE79668 <- cbind(clinical[same,],GSE79668L[same,])
  rownames(GSE79668) <- GSE79668[,9]
  outTab1=rbind(geneNames=colnames(GSE79668), GSE79668)
  write.table(outTab1, file = './00.data/PAAD/04.meger/GSE79668_ALL.txt', sep="\t", quote=F, col.names=F)
  GSE79668 <- GSE79668[,-c(1,4:9)]
}
##GSE62452
if(T){
  geneExp <- read.table('./00.data/PAAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  group2 <- grep('GSE62452.*',colnames(geneExp))
  GSE62452 <- t(geneExp[,group2])
  GSE62452L <- data.frame(GSE62452) 
  id <- row.names(GSE62452)
  GSE62452L  <- cbind(id,GSE62452L)
  rownames(GSE62452L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE62452L))
  cli <- read.table('./00.data/PAAD/02.GSE62452/clinical.txt',header = T,sep = '\t',,check.names = F,row.names = 1)
  #cli$fustat <- ifelse(cli$fustat=='alive','0','1')
  clinical <- cli[,c(1,3:6)]
  clinical$futime=clinical$futime/12
  same <- intersect(row.names(clinical),row.names(GSE62452L))
  GSE62452 <- cbind(clinical[same,],GSE62452L[same,])
  rownames(GSE62452) <- GSE62452[,6]
  outTab1=rbind(geneNames=colnames(GSE62452), GSE62452)
  write.table(outTab1, file = './00.data/PAAD/04.meger/GSE62452_ALL.txt', sep="\t", quote=F, col.names=F)
  GSE62452 <- GSE62452[,-c(1,4:6)]
}
###TCGA##
if(T){
  tcga_cli <- read.table('./00.data/PAAD/01.TCGA/TCGA_PAAD_cli.txt',header = T,sep = '\t',,check.names = F,row.names = 1)
  tcga_cli$futime <- tcga_cli$OS.time
  tcga_cli$fustat <- tcga_cli$OS
  cli<- tcga_cli
  cli <-tcga_cli[,c(10,11,3:9)] 
  cli$futime=cli$futime/365
  colnames(cli) <- c('futime','fustat','age','gender','stage','grade','M','N','T')
  geneExp <- read.table('./00.data/PAAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
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
  write.table(outTab2, file = './00.data/PAAD/04.meger/TCGA_ALL.txt', sep="\t", quote=F, col.names=F)
  TCGA <- TCGA[,-c(3:10)]
}


PAAD_m6A <- rbind(GSE62452,GSE79668,TCGA)
outTab3=rbind(IDs=colnames(PAAD_m6A),PAAD_m6A)
write.table(outTab3, file="./00.data/PAAD/04.meger/PAAD_m6A.txt", sep="\t", quote=F, col.names=F)

if(T){
  #对基因进行循环，找出预后相关的基因
  rt <- read.table('./00.data/PAAD/04.meger/PAAD_m6A.txt', header=T, sep="\t", check.names=F,row.names = 1)
  picDir <- './02.step2/01.survial_m6A/'
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
  write.table(outTab,file="./00.data/PAAD/05.survival/uniCox.txt",sep="\t",row.names=F,quote=F)
}
