library(limma)
library(lars)
library(survival)
library(glmnet)
library(survivalROC)
library(survminer)
library(pROC)
library(glmSparseNet)
library(loose.rock)
library(timeROC)

#读取表达文件，并对输入文件整理
rt=read.table('./00.data/READ/09.VENN/interGeneExp.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data1=t(data)
#rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))
cli=read.table("./00.data/READ/04.meger/READ_m6A.txt", header=T, sep="\t", check.names=F, row.names=1)
cli <- cli[,c(1:2)]

#数据合并
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt1=cbind(cli, data1)
group1 <-grep('TCGA.*',row.names(rt1))
rt <- rt1[group1,]
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
write.table(outTab,file="./00.data/READ/11.Model/uniCox.xls",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="./00.data/READ/11.Model/uniSigExp.txt",sep="\t",row.names=F,quote=F)

##lasso##
rt=read.table("./00.data/READ/11.Model/uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件
rt$futime[rt$futime<=0]=1

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("./05.step5/09.model/lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("./05.step5/09.model/cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
TCGA_train=rt[,lassoGene]
#myFun=function(x){crossprod(as.numeric(x),actCoef)}
#trainScore=apply(TCGA_train,1,myFun)
trainScore=predict(cvfit, newx=as.matrix(rt[,c(3:ncol(rt))]), s="lambda.min", type="response")
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore), "high", "low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore))
train=cbind(rt[,outCol], riskScore=as.vector(trainScore), risk)
trainRiskOut=cbind(id=rownames(train), train)
write.table(trainRiskOut, file="./00.data/READ/11.Model/m6Ascore_train.txt", sep="\t", quote=F, col.names=F)














diff=survdiff(Surv(futime, fustat) ~risk, data=train)
pValue=1-pchisq(diff$chisq, df=1)


#ROC曲线下面积
predictTime=3    #预测时间
roc=timeROC(T=train$futime, delta=train$fustat,
            marker=trainScore, cause=1,
            weighting='aalen',
            times=c(predictTime), ROC=TRUE)
roc$AUC[2]







group2 <-grep('GSE17537.*',row.names(rt1))
GSE17537 <- rt1[group2,]
group3 <-grep('GSE17536.*',row.names(rt1))
GSE17536 <- rt1[group3,]
test <- rbind(GSE17537,GSE17536)
write.table(test, file="./00.data/READ/11.Model/test.txt", sep="\t", quote=F)
test <- na.omit(test)
Megre_test <- test[,lassoGene]
test <- read.table('./00.data/READ/11.Model/test.txt',header = T,sep = '\t',check.names = F,row.names = 1)
myFun=function(x){crossprod(as.numeric(x), actCoef)}
testScore=apply(Megre_test, 1, myFun)
risk=as.vector(ifelse(testScore>median(trainScore), "high", "low"))
test=cbind(test[,outCol], riskScore=as.vector(testScore), risk)
diffTest=survdiff(Surv(futime, fustat) ~risk, data=test)
pValueTest=1-pchisq(diffTest$chisq, df=1)
rocTest=timeROC(T=test$futime, delta=test$fustat,
                marker=testScore, cause=1,
                weighting='aalen',
                times=c(1,3,5), ROC=TRUE)	


testScore=apply(Megre_test,1,myFun)
outCol=c("futime","fustat",lassoGene)
outTab1=cbind(TEST[,outCol],riskScore=as.vector(testScore))
scoreOut1=rbind(id=colnames(outTab1), outTab1)
write.table(scoreOut1, file="./00.data/READ/11.Model/m6Ascore_test.txt", sep="\t", quote=F, col.names=F)
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp1=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp1,file="./00.data/PAAD/11.Model/lassoSigExp.txt",sep="\t",row.names=F,quote=F)





if(T){
  #PCA分析
  pca1=prcomp(TCGA_train, scale=TRUE)
  pca2=prcomp(Megre_test, scale=TRUE)
  value1=predict(pca1)
  m6Ascore_train=value1[,1]+value1[,2]
  m6Ascore_train=as.data.frame(m6Ascore_train)
  same1=intersect(row.names(m6Ascore_train), row.names(cli))
  m6Ascore_train=cbind(cli[same1,],m6Ascore_train[same1,])
  colnames(m6Ascore_train) <- c("futime","fustat",'riskScore')
  value2=predict(pca2)
  m6Ascore_test=value2[,1]+value2[,2]
  m6Ascore_test=as.data.frame(m6Ascore_test)
  same2=intersect(row.names(m6Ascore_test), row.names(cli))
  m6Ascore_test=cbind(cli[same2,],m6Ascore_test[same2,])
  colnames(m6Ascore_test) <- c("futime","fustat",'riskScore')
  scoreOut=rbind(id=colnames(m6Ascore_train), m6Ascore_train)
  write.table(scoreOut, file="./00.data/COAD/11.Model/m6Ascore_train.txt", sep="\t", quote=F, col.names=F)
  scoreOut1=rbind(id=colnames(m6Ascore_test), m6Ascore_test)
  write.table(scoreOut1, file="./00.data/COAD/11.Model/m6Ascore_test.txt", sep="\t", quote=F, col.names=F)
  
  
  sur_train=read.table("./00.data/READ/11.Model/m6Ascore_train.txt", header=T, sep="\t", check.names=F, row.names=1)
  res.cut=surv_cutpoint(sur_train, time="futime", event="fustat", variables=c("riskScore"))
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(sur_train[,"riskScore"]<=cutoff, "Low", "High")
  table(Type)
  sur_train$group=Type
  sur_test=read.table("./00.data/READ/11.Model/m6Ascore_test.txt", header=T, sep="\t", check.names=F, row.names=1)
  Type1=ifelse(sur_test[,"riskScore"]<=cutoff, "Low", "High")
  table(Type1)
  sur_test$group=Type1
  outTab=rbind(id=colnames(sur_train), sur_train)
  write.table(outTab, file="./00.data/READ/11.Model/m6Ascore_train_group.txt", sep="\t", quote=F, col.names=F)
  outTab1=rbind(id=colnames(sur_test), sur_test)
  write.table(outTab1, file="./00.data/READ/11.Model/m6Ascore_test_group.txt", sep="\t", quote=F, col.names=F)
  
  
  surplot <- function(data,plot){
    data <- read.table(data,header=T, sep="\t", check.names=F, row.names=1)
    #计算高低风险组生存差异
    data$group=factor(data$group, levels=c("Low", "High"))
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
  data <- "./00.data/READ/11.Model/m6Ascore_train_group.txt"
  plot <- './05.step5/09.model/survival_train.pdf'
  surplot(data,plot)
  
  data <- "./00.data/READ/11.Model/m6Ascore_test_group.txt"
  plot <- './05.step5/09.model/survival_test.pdf'
  surplot(data,plot)
  
  
  rocplot <- function(inputFile,outFile){
    
    
    #读入文件，整理
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    y=colnames(rt)[length(colnames(rt))]
    x= "riskScore"
    
    #绘制
    rocobj1=roc(rt[,y], as.vector(rt[,x]))
    pdf(file=outFile,width=5,height=5)
    plot(rocobj1, print.auc=TRUE, col="red")
    dev.off()
  }
  inputFile <- "./00.data/READ/11.Model/m6Ascore_train_group.txt"
  outPdf <- './05.step5/09.model/rocTrain.pdf'
  rocplot(inputFile,outPdf)
  
  inputFile <- "./00.data/READ/11.Model/m6Ascore_test_group.txt"
  outPdf <- './05.step5/09.model/rocTest.pdf'
  rocplot(inputFile,outPdf)
}

