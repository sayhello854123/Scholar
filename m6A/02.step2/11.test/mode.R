library(limma)
library(survival)
library(glmnet)
library(survivalROC)
library(survminer)
#读取表达文件，并对输入文件整理

rt=read.table('./02.step2/11.test/data.train.txt', header=T, sep="\t", check.names=F,row.names = 1)
rt1 <- read.table('./02.step2/11.test/data.test.txt', header=T, sep="\t", check.names=F,row.names = 1)
#对基因进行循环，找出预后相关的基因
pFilter=0.05
outTab=data.frame()
sigGenes=c("futime","fustat")
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
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
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
  }
}

write.table(outTab,file="./02.step2/11.test/uniCox.xls",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="./02.step2/11.test/uniSigExp.txt",sep="\t",row.names=F,quote=F)

##lasso##
rt=read.table("./02.step2/11.test/uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件
rt$futime[rt$futime<=0]=1

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("./02.step2/11.test/lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("./02.step2/11.test/cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
TCGA_train=rt[,lassoGene]
meger_test <- rt1[,lassoGene]
lassoGene=c("futime","fustat",lassoGene)
TCGA_lassoSigExp=rt[,lassoGene]
meger_lassoSigExp <- rt1[,lassoGene]
lassoSigExp1=cbind(id=row.names(TCGA_lassoSigExp),TCGA_lassoSigExp)
write.table(lassoSigExp1,file="./02.step2/11.test/TCGA_lassoSigExp.txt",sep="\t",row.names=F,quote=F)

#PCA分析
######Train######
pca1=prcomp(TCGA_train, scale=TRUE)
value1=predict(pca1)
m6Ascore_train=value1[,1]+value1[,2]
m6Ascore_train=as.data.frame(m6Ascore_train)
same1=intersect(row.names(TCGA_lassoSigExp), row.names(m6Ascore_train))
m6AscoreTrain=cbind(TCGA_lassoSigExp[same1,],m6Ascore_train[same1,])
colnames(m6AscoreTrain)[ncol(m6AscoreTrain)] <- 'riskscore'
medianTrainRisk=median(m6AscoreTrain$riskscore)
risk=as.vector(ifelse(m6AscoreTrain$riskscore>medianTrainRisk,"high","low"))
m6AscoreTrain <- cbind(m6AscoreTrain, risk)
m6Ascore_trainOut=cbind(id=row.names(m6AscoreTrain), m6AscoreTrain)
write.table(m6Ascore_trainOut,file="./02.step2/11.test/train.txt",sep="\t",row.names=F,quote=F)

######TEST######
pca2=prcomp(meger_test, scale=TRUE)
value2=predict(pca2)
m6Ascore_test=value2[,1]+value2[,2]
m6Ascore_test=as.data.frame(m6Ascore_test)
same2=intersect(row.names(m6Ascore_test), row.names(meger_lassoSigExp))
m6AscoreTest=cbind(meger_lassoSigExp[same2,],m6Ascore_test[same2,])
colnames(m6AscoreTest)[ncol(m6AscoreTest)] <- 'riskscore'
risk=as.vector(ifelse(m6AscoreTest$riskscore>medianTrainRisk,"high","low"))
m6AscoreTest<- cbind(m6AscoreTest, risk)
m6Ascore_testOut=cbind(id=row.names(m6AscoreTest), m6AscoreTest)
write.table(m6Ascore_testOut,file="./02.step2/11.test/test.txt",sep="\t",row.names=F,quote=F)


###计算生存曲线
surplot <- function(data,plot) {
  data <- read.table(data,header=T, sep="\t", check.names=F, row.names=1)
  #计算高低风险组生存差异
  data$group=factor(data$risk, levels=c("low", "high"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = data)
  length=length(levels(factor(data[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = data)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=data,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=levels(factor(data[,"group"])),
                     #legend = c(0.8, 0.8),
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     #surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)
  
  #保存图片
  pdf(file=plot, onefile = FALSE, width=7, height=5.5)
  print(surPlot)
  dev.off()
}
data <- "./04.step4/12.test/train.txt"
plot <- './04.step4/12.test/survival_train.pdf'
surplot(data,plot)

data <- "./04.step4/12.test/test.txt"
plot <- './04.step4/12.test/survival_test.pdf'
surplot(data,plot)

#ROC曲线下面积
rocPlot=function(inputFile,outPdf){
  risk=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
  #定义颜色
  bioCol=c('#FF0000','#00FF00','#00FFFF')
  roc=timeROC(T=risk$futime, delta=risk$fustat,
              marker=risk$riskScore, cause=1,
              times=c(1,3,5), ROC=TRUE)
  pdf(file=outPdf, width=5, height=5)
  plot(roc,time=1,col=bioCol[1],title=FALSE,lwd=2)
  plot(roc,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
  plot(roc,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",roc$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",roc$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",roc$AUC[3]))),
         col=bioCol[1:3], lwd=2, bty = 'n')
  dev.off()
  
}
inputFile <- "./04.step4/12.test/train.txt"
outPdf <- './04.step4/12.test/rocTrain.pdf'
rocPlot(inputFile,outPdf)

inputFile <- "./04.step4/12.test/test.txt"
outPdf <- './04.step4/12.test/rocTest.pdf'
rocPlot(inputFile,outPdf)
