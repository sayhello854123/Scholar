library(limma)
library(survival)
library(survminer)
library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")
##GSE39582
if(T){
  geneExp <- read.table('./00.data/COAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  group2 <- grep('GSE39582.*',colnames(geneExp))
  GSE39582 <- t(geneExp[,group2])
  GSE39582L <- data.frame(GSE39582) 
  id <- row.names(GSE39582)
  GSE39582L  <- cbind(id,GSE39582L)
  rownames(GSE39582L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE39582L))
  cli <- read.table('./00.data/COAD/02.GSE39582/clinical.txt',header = T,sep = '\t',)
  rownames(cli) <- cli[,1]
  table(cli$fustat)
  #cli$fustat <- ifelse(cli$fustat=='Alive','0','1')
  clinical <- cli[,c(1,2:9)]
  clinical$futime=clinical$futime/12
  same <- intersect(row.names(clinical),row.names(GSE39582L))
  GSE39582 <- cbind(clinical[same,],GSE39582L[same,])
  rownames(GSE39582) <- GSE39582[,10]
  outTab1=rbind(geneNames=colnames(GSE39582), GSE39582)
  write.table(outTab1, file = './00.data/COAD/04.meger/GSE39582_ALL.txt', sep="\t", quote=F, col.names=F)
  GSE39582 <- GSE39582[,-c(1,4:10)]
}
##GSE17538
if(T){
  geneExp <- read.table('./00.data/COAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  group3 <- grep('GSE17538.*',colnames(geneExp))
  GSE17538 <- t(geneExp[,group3])
  GSE17538L <- data.frame(GSE17538) 
  id <- row.names(GSE17538)
  GSE17538L  <- cbind(id,GSE17538L)
  rownames(GSE17538L)=gsub("(.*?)\\_(.*?)", "\\2", rownames(GSE17538L))
  cli <- read.table('./00.data/COAD/03.GSE17538/clinical.txt',header = T,sep = '\t')
  rownames(cli) <- cli[,1]
  #cli <- cli[,-c(8,9)]
  table(cli$fustat)
  cli$fustat <- ifelse(cli$fustat=='no death','0','1')
  clinical <- cli[,c(1,2:7)]
  clinical$futime=clinical$futime/12
  same <- intersect(row.names(clinical),row.names(GSE17538L))
  GSE17538 <- cbind(clinical[same,],GSE17538L[same,])
  rownames(GSE17538) <- GSE17538[,8]
  outTab1=rbind(geneNames=colnames(GSE17538), GSE17538)
  write.table(outTab1, file = './00.data/COAD/04.meger/GSE17538_ALL.txt', sep="\t", quote=F, col.names=F)
  GSE17538 <- GSE17538[,-c(1,4:8)]
}
#TCGA
if(T){
  tcga_cli <- read.table('./00.data/COAD/01.TCGA/TCGA_cli.txt',header = T,sep = '\t',,check.names = F,row.names = 1)
  tcga_cli$futime <- tcga_cli$OS.time
  tcga_cli$fustat <- tcga_cli$OS
  cli<- tcga_cli
  cli <-tcga_cli[,c(11,12,3:10)] 
  cli$futime=cli$futime/365
  colnames(cli) <- c('futime','fustat','age','gender','stage','grade','T','N','M','Statu')
  geneExp <- read.table('./00.data/COAD/04.meger/m6A_merge.txt',sep = '\t',header = T,check.names = F,row.names = 1)
  colnames(geneExp)
  group3 <- grep('TCGA.*',colnames(geneExp))
  TCGA <- t(geneExp[,group3])
  TCGA1 <- TCGA
  id <- row.names(TCGA)
  TCGA1 <- cbind(id,TCGA1)
  rownames(TCGA1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(TCGA1))
  same <- intersect(row.names(cli),row.names(TCGA1))
  TCGA<- cbind(cli[same,],TCGA1[same,])
  rownames(TCGA) <- TCGA[,11]
  outTab2=rbind(geneNames=colnames(TCGA), TCGA)
  write.table(outTab2, file = './00.data/COAD/04.meger/TCGA_ALL.txt', sep="\t", quote=F, col.names=F)
  TCGA <- TCGA[,-c(3:11)]
}
PAAD_m6A <- rbind(GSE39582,GSE17538,TCGA)
outTab3=rbind(IDs=colnames(PAAD_m6A),PAAD_m6A)
write.table(outTab3, file="./00.data/COAD/04.meger/COAD_m6A.txt", sep="\t", quote=F, col.names=F)

if(T){
  #对基因进行循环，找出预后相关的基因
  rt <- read.table('./00.data/COAD/04.meger/COAD_m6A.txt', header=T, sep="\t", check.names=F,row.names = 1)
  picDir <- './04.step4/01.survial_m6A/'
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
  write.table(outTab,file="./00.data/COAD/05.survival/uniCox.txt",sep="\t",row.names=F,quote=F)
}