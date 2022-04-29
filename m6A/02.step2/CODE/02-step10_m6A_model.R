library(limma)
library(survival)
library(glmnet)
library(survivalROC)
library(survminer)
#读取表达文件，并对输入文件整理
rt=read.table('./00.data/PAAD/09.PAAD_VENN/interGeneExp.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data1=t(data)
#rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))
cli=read.table("./00.data/PAAD/04.meger/PAAD_m6A.txt", header=T, sep="\t", check.names=F, row.names=1)
cli <- cli[,c(1:2)]

#数据合并
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt1=cbind(cli, data1)
group1 <-grep('TCGA.*',row.names(rt1))
rt <- rt1[group1,]
#对基因进行循环，找出预后相关的基因
pFilter=0.01
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
write.table(outTab,file="./00.data/PAAD/11.Model/uniCox.xls",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="./00.data/PAAD/11.Model/uniSigExp.txt",sep="\t",row.names=F,quote=F)

##lasso##
rt=read.table("./00.data/PAAD/11.Model/uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件
rt$futime[rt$futime<=0]=1

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("./02.step2/09.model/lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("./02.step2/09.model/cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
TCGA_train=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(TCGA_train,1,myFun)
outCol=c("futime","fustat",lassoGene)
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore))
scoreOut=rbind(id=colnames(outTab), outTab)
write.table(scoreOut, file="./00.data/PAAD/11.Model/m6Ascore_train.txt", sep="\t", quote=F, col.names=F)
group2 <-grep('GSE62452.*',row.names(rt1))
GSE62452 <- rt1[group2,]
group3 <-grep('GSE79668.*',row.names(rt1))
GSE79668 <- rt1[group3,]
TEST <- rbind(GSE62452,GSE79668)
Megre_test <- TEST[,lassoGene]
testScore=apply(Megre_test,1,myFun)
outCol=c("futime","fustat",lassoGene)
outTab1=cbind(TEST[,outCol],riskScore=as.vector(testScore))
scoreOut1=rbind(id=colnames(outTab1), outTab1)
write.table(scoreOut1, file="./00.data/PAAD/11.Model/m6Ascore_test.txt", sep="\t", quote=F, col.names=F)
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
colnames(m6Ascore_train) <- c("futime","fustat",'train')
value2=predict(pca2)
m6Ascore_test=value2[,1]+value2[,2]
m6Ascore_test=as.data.frame(m6Ascore_test)
same2=intersect(row.names(m6Ascore_test), row.names(cli))
m6Ascore_test=cbind(cli[same2,],m6Ascore_test[same2,])
colnames(m6Ascore_test) <- c("futime","fustat",'test')
scoreOut=rbind(id=colnames(m6Ascore_train), m6Ascore_train)
write.table(scoreOut, file="./00.data/PAAD/11.Model/m6Ascore_train.txt", sep="\t", quote=F, col.names=F)
scoreOut1=rbind(id=colnames(m6Ascore_test), m6Ascore_test)
write.table(scoreOut1, file="./00.data/PAAD/11.Model/m6Ascore_test.txt", sep="\t", quote=F, col.names=F)


sur_train=read.table("./00.data/PAAD/11.Model/m6Ascore_train.txt", header=T, sep="\t", check.names=F, row.names=1)
res.cut=surv_cutpoint(sur_train, time="futime", event="fustat", variables=c("riskScore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(sur_train[,"riskScore"]<=cutoff, "Low", "High")
table(Type)
sur_train$group=Type
sur_test=read.table("./00.data/PAAD/11.Model/m6Ascore_test.txt", header=T, sep="\t", check.names=F, row.names=1)
Type1=ifelse(sur_test[,"riskScore"]<=cutoff, "Low", "High")
table(Type1)
sur_test$group=Type1
outTab=rbind(id=colnames(sur_train), sur_train)
write.table(outTab, file="./00.data/PAAD/11.Model/m6Ascore_train_group.txt", sep="\t", quote=F, col.names=F)
outTab1=rbind(id=colnames(sur_test), sur_test)
write.table(outTab1, file="./00.data/PAAD/11.Model/m6Ascore_test_group.txt", sep="\t", quote=F, col.names=F)


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
data <- "./00.data/PAAD/11.Model/m6Ascore_train_group.txt"
plot <- './02.step2/09.model/survival_train.pdf'
surplot(data,plot)

data <- "./00.data/PAAD/11.Model/m6Ascore_test_group.txt"
plot <- './02.step2/09.model/survival_test.pdf'
surplot(data,plot)


rocPlot=function(inputFile,outPdf){
  rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
  pdf(file=outPdf,width=5.5,height=5.5)
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], 
                  predict.time =1, method="KM")
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
       xlab="False positive rate", ylab="True positive rate",
       main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  dev.off()
}
inputFile <- "./00.data/PAAD/11.Model/m6Ascore_train_group.txt"
outPdf <- './02.step2/09.model/rocTrain.pdf'
rocPlot(inputFile,outPdf)

inputFile <- "./00.data/PAAD/11.Model/m6Ascore_test_group.txt"
outPdf <- './02.step2/09.model/rocTest.pdf'
rocPlot(inputFile,outPdf)
}

