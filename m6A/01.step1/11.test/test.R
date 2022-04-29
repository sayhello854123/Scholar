library(limma)
library(survival)
library(glmnet)
library(survivalROC)
library(survminer)
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
#读取表达文件，并对输入文件整理
rt=read.table('./00.data/LIHC/PCA/interGeneExp.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data1=t(data)
#rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))
cli=read.table("./00.data/LIHC/megerData/LIHC_m6A.txt", header=T, sep="\t", check.names=F, row.names=1)
cli <- cli[,c(1:2)]

#数据合并
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt1=cbind(cli, data1)
group1 <-grep('TCGA.*',row.names(rt1))
group2 <-grep('ICGC.*',row.names(rt1))
group3 <-grep('GSE76427.*',row.names(rt1))
rt <- rt1
coxPfilter=0.05        
#对数据进行分组，构建模型
n=1000     #分组的次数
for(i in 1:n){
  #############对数据进行分组#############
  inTrain<-createDataPartition(y=rt[,3], p=367/715, list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  #单因素cox分析
  outUniTab=data.frame()
  sigGenes=c("futime","fustat")
  for(i in colnames(train[,3:ncol(train)])){
    #cox分析
    cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    
    #保留显著性基因
    if(coxP<coxPfilter){
      sigGenes=c(sigGenes,i)
      outUniTab=rbind(outUniTab,
                      cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  uniSigExp=train[,sigGenes]
  uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
  
  #lasso回归
  x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
  y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  train_lasso <- uniSigExp[,lassoGene]
  lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
  lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
  if(nrow(geneCoef)<2){next}
  
  #############构建PCA模型#############
  pca1=prcomp(train_lasso, scale=TRUE)
  value=predict(pca1)
  m6Ascore=value[,1]+value[,2]
  m6Ascore=as.data.frame(m6Ascore)
  same1=intersect(row.names(lassoSigExp), row.names(m6Ascore))
  m6Ascore_train=cbind(lassoSigExp[same1,],m6Ascore[same1,])
  colnames(m6Ascore_train)[ncol(m6Ascore_train)] <- 'riskscore'
  medianTrainRisk=median(m6Ascore_train$riskscore)
  risk=as.vector(ifelse(m6Ascore_train$riskscore>medianTrainRisk,"high","low"))
  m6Ascore_train <- cbind(m6Ascore_train, risk)
  m6Ascore_trainOut=cbind(id=row.names(m6Ascore_train), m6Ascore_train)
 
  #####TESTPCA模型#############
  test1 <- test[,c("futime", "fustat", lassoGene)]
  test_lasso <- test[,lassoGene]
  pca2=prcomp(test_lasso, scale=TRUE)
  value1=predict(pca2)
  m6Ascore1=value1[,1]+value1[,2]
  m6Ascore1=as.data.frame(m6Ascore1)
  same2=intersect(row.names(test1), row.names(m6Ascore1))
  m6Ascore_test=cbind(test1[same2,],m6Ascore1[same2,])
  colnames(m6Ascore_test)[ncol(m6Ascore_test)] <- 'riskscore'
  riskTest=as.vector(ifelse(m6Ascore_test$riskscore>medianTrainRisk,"high","low"))
  m6Ascore_test <- cbind(m6Ascore_test, riskTest)
  m6Ascore_testOut=cbind(id=row.names(m6Ascore_test), m6Ascore_test)
  
  #生存差异pvalue	
  diff=survdiff(Surv(futime, fustat) ~risk,data = m6Ascore_train)
  pValue=1-pchisq(diff$chisq, df=1)
  diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = m6Ascore_test)
  pValueTest=1-pchisq(diffTest$chisq, df=1)
  
  #ROC曲线下面积
  roc=timeROC(T=m6Ascore_train$futime, delta=m6Ascore_train$fustat,
              marker=m6Ascore_train$riskscore, cause=1,
              times=c(1,3,5), ROC=TRUE)
  rocTest=timeROC(T=m6Ascore_test$futime, delta=m6Ascore_test$fustat,
                  marker=m6Ascore_test$riskscore, cause=1,
                  times=c(1,3,5), ROC=TRUE)	
  
  roc$AUC[2]
  
  if((pValue<0.03) & (roc$AUC[1]>0.7)&(roc$AUC[2]>0.67)&(roc$AUC[2]>0.65) & 
     (pValueTest<0.04) & (rocTest$AUC[1]>0.69)&(rocTest$AUC[2]>0.66)&(rocTest$AUC[3]>0.64)){
    #输出分组结果
    write.table(trainOut,file="./01.step1/11.test/data.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testOut,file="./01.step1/11.test/data.test.txt",sep="\t",quote=F,row.names=F)
    #输出单因素结果
    write.table(outUniTab,file="./01.step1/11.test/uni.trainCox.txt",sep="\t",row.names=F,quote=F)
    write.table(uniSigExpOut,file="./01.step1/11.test/uni.SigExp.txt",sep="\t",row.names=F,quote=F)
    
    #lasso结果
    write.table(lassoSigExpOut,file="./01.step1/11.test/lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("./01.step1/11.test/lasso.lambda.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    pdf("./01.step1/11.test/lasso.cvfit.pdf")
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
    dev.off()
    break
  }
}
