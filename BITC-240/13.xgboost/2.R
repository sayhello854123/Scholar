library("xgboost")
library(limma)
library(caret)
library(data.table)
library(Matrix)
library(ROCR)
if(T){
rt=read.table('./merge.txt', header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
dim(data)
gene <- read.table('gene.txt',sep = '\t')
same <- intersect(as.vector(gene[,1]),row.names(data))
data <- data[same,]

#读取目录下所有"s1.txt"结尾的文件
sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
  rt=read.table(file, header=F, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  sampleName1=c(sampleName1, uniqGene)
}

#读取目录下所有"s2.txt"结尾的文件
sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
  rt=read.table(file, header=F, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  sampleName2=c(sampleName2, uniqGene)
}

#提取实验组和对照组的数据
conData=data[,sampleName1]
treatData=data[,sampleName2]
conNum=ncol(conData)
treatNum=ncol(treatData)
data=t(cbind(conData,treatData))
write.csv(data,file = 'test.csv')

####输入的测试集
train1 <- data.matrix(data)
train <- Matrix(train1,sparse = T )
lab <- read.csv('./test_lab.csv')
train_lab <- lab$type
traindata <- list(data = train,label=train_lab )
dtrain <- xgb.DMatrix(data = traindata$data,label = traindata$label)

}
if(T){
####构建验证集
rt<-data.table::fread(file ='TCGA-STAD.htseq_fpkm.tsv',data.table = F)
gencode <- data.table::fread('gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
rownames(rt) <- rt[,1]
rt <- rt[,-1]
data <- 2^rt-1
same1 <- intersect(row.names(gencode),row.names(data))
length(same1)
data1 <- cbind(gencode[same1,],data[same1,])
data1 <- data1[,-c(1,3:6)]
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

#FPKM转换为TPM
fpkmToTpm=function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(rt, 2, fpkmToTpm)
test1 <- tpm[same,]
group_list=ifelse(as.numeric(substr(colnames(test1),14,15)) < 10,'tumor','normal')
table(group_list)
Tumor <- test1[,group_list == "tumor"]
Normal <- test1[,group_list == "normal"]
test2 <- t(cbind(Normal,Tumor))
test1_lab1 <- cbind(colnames(Normal),rep(0,length(colnames(Normal))))
test1_lab2 <- cbind(colnames(Tumor),rep(1,length(colnames(Tumor))))
test_lab = data.frame(rbind(test1_lab1,test1_lab2))

####输入的验证集
test2 <- data.matrix(test2)
test <- Matrix(test2,sparse = T )
test_lab <- test_lab$X2
testdata <- list(data = test,label=test_lab)
dtest <- xgb.DMatrix(data = testdata$data,label = testdata$label)
save(dtrain,dtest,file = 'binary_logistic.RData')
}
##构建模型
modle_xgb <- xgboost(data = dtrain,booster = 'gbtree',
                     max.depth = 100, eta = 0.5, nthread = 2, nrounds = 100,
                     objective = "binary:logistic")
importance <- xgb.importance(train@Dimnames[[2]], model = modle_xgb) 
write.csv(importance,file = './importance.csv')
head(importance)
xgb.ggplot.importance(importance)
##模型预测
pre <- predict(modle_xgb,newdata = dtest)
#重要重要性排序 
confusionMatrix(as.factor(pre),as.factor(testdata$label))
xgboost_roc <- roc(test_lab,as.numeric(pre))
plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2),grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,main='ROC curve')
