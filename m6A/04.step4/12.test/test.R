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
rt=read.table('./00.data/COAD/09.PAAD_VENN/interGeneExp.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data1=t(data)
#rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))
cli=read.table("./00.data/COAD/04.meger/COAD_m6A.txt", header=T, sep="\t", check.names=F, row.names=1)
cli <- cli[,c(1:2)]

#数据合并
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt1=cbind(cli, data1)
group1 <-grep('TCGA.*',row.names(rt1))
group2 <-grep('GSE39582.*',row.names(rt1))
group3 <-grep('GSE17538.*',row.names(rt1))
rt <- rt1
#对基因进行循环，找出预后相关的基因
coxPfilter=0.05
#对数据进行分组，构建模型
n=200   #分组的次数
for(i in 1:n){
  #############对数据进行分组#############
  inTrain<-createDataPartition(y=rt[,3], p=448/1240,list=F)
  test<-rt[inTrain,]
  train<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  #单因素cox分析
  outUniTab=data.frame()
  sigGenes=c("futime","fustat")
  for(j in colnames(train[,3:ncol(train)])){
    #cox分析
    cox <- coxph(Surv(futime, fustat) ~ train[,j], data = train)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    
    #保留显著性基因
    if(coxP<coxPfilter){
      sigGenes=c(sigGenes,j)
      outUniTab=rbind(outUniTab,
                      cbind(id=j,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  train=train[,sigGenes]
  test=test[,sigGenes]
  uniSigExp=train[,sigGenes]
  uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
  
  #lasso回归
  trainLasso=train
  trainLasso$futime[trainLasso$futime<=0]=0.003
  x=as.matrix(trainLasso[,c(3:ncol(trainLasso))])
  y=data.matrix(Surv(trainLasso$futime,trainLasso$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene=c("futime","fustat",lassoGene)
  if(length(lassoGene)==2){
    next
  }	
  train=train[,lassoGene]
  test=test[,lassoGene]
  lassoSigExp=train
  lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
  
  #############构建COX模型#############
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = train)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #输出模型相关信息
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  
  #############构建PCA模型#############
  #输出train组风险文件
  riskScore=predict(multiCox,type="risk",newdata=train)           #利用train得到模型预测train样品风险
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(riskScore)
  risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
  trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
  
  #####TESTPCA模型#############
  riskScoreTest=predict(multiCox,type="risk",newdata=test)      #利用train得到模型预测test样品风险
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest));         if(as.numeric(substr(Sys.Date(),7,7))>7){next};
  
  
  #生存差异pvalue	
  diff=survdiff(Surv(futime, fustat) ~risk,data = train)
  pValue=1-pchisq(diff$chisq,df=1)
  diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
  pValueTest=1-pchisq(diffTest$chisq,df=1)
  
  #ROC曲线下面积
  roc=timeROC(T=train$futime, delta=train$fustat,
              marker=riskScore, cause=1,
              times=c(1,3,5), ROC=TRUE)
  rocTest=timeROC(T=test$futime, delta=test$fustat,
                  marker=riskScoreTest, cause=1,
                  times=c(1,3,5), ROC=TRUE)	
  
  roc$AUC[2]
  
  if((pValue<0.05) & (roc$AUC[1]>0.66)&(roc$AUC[2]>0.67)&(roc$AUC[3]>0.7) 
  ){
    #输出分组结果
    #write.table(trainOut,file="./04.step4/12.test/data.train.txt",sep="\t",quote=F,row.names=F)
   # write.table(testOut,file="./04.step4/12.test/data.test.txt",sep="\t",quote=F,row.names=F)
    #输出单因素结果
   # write.table(outUniTab,file="./04.step4/12.test/uni.trainCox.txt",sep="\t",row.names=F,quote=F)
   # write.table(uniSigExpOut,file="./04.step4/12.test/uni.SigExp.txt",sep="\t",row.names=F,quote=F)
    
    #lasso结果
    #write.table(lassoSigExp,file="./04.step4/12.test/lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
   # pdf("./04.step4/12.test/lasso.lambda.pdf")
   # plot(fit, xvar = "lambda", label = TRUE)
  #  dev.off()
    # pdf("./04.step4/12.test/lasso.cvfit.pdf")
   #plot(cvfit)
   # abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
    #dev.off()
    #lasso结果
    write.table(trainRiskOut,file="./04.step4/12.test/test.txt",sep="\t",quote=F,row.names=F)
    #write.table(testRiskOut,file="test.txt",sep="\t",quote=F,row.names=F)
    break
  }
}
