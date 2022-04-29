library(corrplot)
library(survival)
library(survminer)
if(T){
train=read.table("./00.data/LIHC/Model/m6Ascore_train.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(train) = c("futime", "fustat","m6A")
test=read.table("./00.data/LIHC/Model/m6Ascore_test.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(test) = c("futime", "fustat","m6A")
score <- rbind(train,test)
metabolize=read.table('./00.data/LIHC/cluster/ssGSEA.result.txt', header=T, sep="\t", check.names=F, row.names=1)
metabolize=t(metabolize)

#数据合并
sameSample=intersect(row.names(score), row.names(metabolize))
data=cbind(score[sameSample,,drop=F], metabolize[sameSample,,drop=F])
data <- data[,-c(1,2)]

#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

size=c()
for(i in 1:ncol(M)){
  size=c(size,M[1:i,i])
}
size=as.numeric(size)

col1 <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow",
                           "#FF7F00","red","#7F0000"))
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F")) 
pdata(file="./01.step1/09.m6aLIHC_COR/corpot.pdata",width=7,height=7)
corrplot(M, method = "circle",
         type = "upper",
         number.cex = abs(size),
         col=col1(200),
         tl.col="black",
         addCoef.col="black"
)
dev.off()



train=read.table("./00.data/LIHC/Model/m6Ascore_train.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(train) = c("futime", "fustat","m6A")
test=read.table("./00.data/LIHC/Model/m6Ascore_test.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(test) = c("futime", "fustat","m6A")
score <- rbind(train,test)
metabolize=read.table('./00.data/LIHC/cluster/ssGSEA.result.txt', header=T, sep="\t", check.names=F, row.names=1)
metabolize=t(metabolize)
}
if(T){
#数据合并
sameSample=intersect(row.names(score), row.names(metabolize))
data=cbind(score[sameSample,,drop=F], metabolize[sameSample,,drop=F])

#获取最优cutoff
res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("m6A"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"m6A"]<=cutoff, 0, 1)
table(Type)
data$m6Agroup=Type

#获取最优cutoff
res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("AMINOACID"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"AMINOACID"]<=cutoff, 0, 1)
table(Type)
data$AMINOACIDgroup=Type

data$bestGroup <- paste0(data$m6Agroup,data$AMINOACIDgroup)
#按照两个基因表达量高低的排列组合分成4组
data$bestGroup <- ifelse(data$bestGroup=='10','1',ifelse(data$bestGroup=='00','2',ifelse(data$bestGroup=="01","3","4")))
head(data)

data$group<-data$bestGroup
#定义结局变量，生存时间和结局
y <- Surv(data$futime, data$fustat==1)

#logrank检验两两比较的结果
comp <- pairwise_survdiff(Surv(futime,fustat) ~ group, data = data)
pvalue<-as.vector(unlist(comp$p.value))

#将p value和两两比较的组别生成两列的数据集
name <- array(dim=c(3,3))
for(i in 1:3){
  for(j in 2:4){
    name[i,j-1] = print(paste(i,j,sep = " vs "))
  }
}
pvalue_name <- as.vector(t(name))
logrank <- data.frame(pvalue_name,pvalue)
logrank

#挑选p value小于0.05的记录
logrank_sig <- subset(logrank, pvalue<0.05)
#如果p值太小，就写“<0.0001”，否则保留小数点后4位
logrank_sig$pvalue <- lapply(logrank_sig$pvalue,function(i)
  ifelse (i<0.0001,"<0.0001",round(i, 4)))
logrank_sig

#进行COX回归，导出每组对应的HR值，并将1赋值给对照组
coxph.fit <- coxph(y ~ as.factor(group), data = data) 
hr <- round(coef(summary(coxph.fit))[,2],3)
HR <- c(1,as.vector(unlist(hr)))


#KM曲线的绘制
kmfit <- survfit(y~data$group,)

#写legend
A = c("Low; ", "High;", "Low; ", "High;")
B = c("Low; ", "Low; ", "High;", "High;")
Group = c(1,2,3,4)
text.legend1 <- paste0(colnames(data)[3], " = ", A, colnames(data)[4], " = ", B, " Group = ", Group, ", HR = ", HR)
text.legend2 <- paste0(logrank_sig$pvalue_name," : ",logrank_sig$pvalue)

#自定义足够你用的颜色
mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#画图，并保存到pdf文件
pdf("./01.step1/09.m6aLIHC_COR/nSurv.pdf",width = 7,height = 6)
par(xpd = T, 
    mar = par()$mar + c(0,0,0,5)#,#在右侧给图例留出空间
    #cex.axis=0.8, cex.lab=1, cex.main=1, cex.sub=1 #修改字体大小
)
plot(kmfit, 
     col= mycol, 
     lwd = 1.4,#线的粗细
     #pch = "o", #如果你也想让观测点的形状是o，就运行这行
     xlab="Years from diagnosis", 
     ylab="Survival Distribution Function", 
     mark.time=TRUE)
legend("bottomleft", lty=c("solid","solid","solid","solid"),
       col=mycol, 
       legend=text.legend1, 
       bty="n", 
       lwd = 2,
       cex=0.8)
legend("topright", 
       inset=c(-0.3,0), #图例画到图外面
       legend=c("Pairwise comparison",text.legend2), 
       bty="n", 
       cex=0.8)
dev.off()
}

#数据合并
picDir <- './01.step1/09.m6aLIHC_COR/'
train=read.table("./00.data/LIHC/Model/m6Ascore_train.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(train) = c("futime", "fustat","m6A")
test=read.table("./00.data/LIHC/Model/m6Ascore_test.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(test) = c("futime", "fustat","m6A")
score <- rbind(train,test)
metabolize=read.table('./00.data/LIHC/cluster/ssGSEA.result.txt', header=T, sep="\t", check.names=F, row.names=1)
metabolize=t(metabolize)
sameSample=intersect(row.names(score), row.names(metabolize))
data=cbind(score[sameSample,,drop=F], metabolize[sameSample,,drop=F])

for(i in colnames(data[,4:ncol(data)])){
  rt <- data[,c("futime","fustat","m6A",i)]
#获取最优cutoff
res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=c("m6A"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(rt[,"m6A"]<=cutoff, 0, 1)
table(Type)
rt$m6Agroup=Type

#获取最优cutoff
res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,i]<=cutoff, 0, 1)
table(Type)
rt$metabogroup=Type

rt$bestGroup <- paste0(rt$m6Agroup,rt$metabogroup)
#按照两个基因表达量高低的排列组合分成4组
rt$bestGroup<- ifelse(rt$bestGroup=='10','1',ifelse(rt$bestGroup=='00','2',ifelse(rt$bestGroup=="01","3","4")))


rt$group<-rt$bestGroup
#定义结局变量，生存时间和结局
y <- Surv(rt$futime, rt$fustat==1)

#logrank检验两两比较的结果
comp <- pairwise_survdiff(Surv(futime,fustat) ~ group, data = rt)
pvalue<-as.vector(unlist(comp$p.value))

#将p value和两两比较的组别生成两列的数据集
name <- array(dim=c(3,3))
for(j in 1:3){
  for(k in 2:4){
    name[j,k-1] = print(paste(j,k,sep = " vs "))
  }
}
pvalue_name <- as.vector(t(name))
logrank <- data.frame(pvalue_name,pvalue)
logrank

#挑选p value小于0.05的记录
logrank_sig <- subset(logrank, pvalue<0.05)
#如果p值太小，就写“<0.0001”，否则保留小数点后4位
logrank_sig$pvalue <- lapply(logrank_sig$pvalue,function(i)
  ifelse (i<0.0001,"<0.0001",round(i, 4)))
logrank_sig

#进行COX回归，导出每组对应的HR值，并将1赋值给对照组
coxph.fit <- coxph(y ~ as.factor(group), data = rt) 
hr <- round(coef(summary(coxph.fit))[,2],3)
HR <- c(1,as.vector(unlist(hr)))


#KM曲线的绘制
kmfit <- survfit(y~rt$group,)

#写legend
A = c("Low; ", "High;", "Low; ", "High;")
B = c("Low; ", "Low; ", "High;", "High;")
Group = c(1,2,3,4)
text.legend1 <- paste0(colnames(rt)[3], " = ", A, i, " = ", B, " Group = ", Group, ", HR = ", HR)
text.legend2 <- paste0(logrank_sig$pvalue_name," : ",logrank_sig$pvalue)

#自定义足够你用的颜色
mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#画图，并保存到pdf文件
filename <- paste0(i,'_survial' ,".pdf",sep='')
outfile <- paste(picDir,filename,sep="")
pdf(outfile,width = 7,height = 6)
par(xpd = T, 
    mar = par()$mar + c(0,0,0,5),#在右侧给图例留出空间
    cex.axis=0.8, cex.lab=1, cex.main=1, cex.sub=1 #修改字体大小
)
plot(kmfit, 
     col= mycol, 
     lwd = 1.4,#线的粗细
     #pch = "o", #如果你也想让观测点的形状是o，就运行这行
     xlab="Years from diagnosis", 
     ylab="Survival Distribution Function", 
     mark.time=TRUE)
legend("bottomleft", lty=c("solid","solid","solid","solid"),
       col=mycol, 
       legend=text.legend1, 
       bty="n", 
       lwd = 2,
       cex=0.8)
legend("topright", 
       inset=c(-0.3,0), #图例画到图外面
       legend=c("Pairwise comparison",text.legend2), 
       bty="n", 
       cex=0.8)
dev.off()
}
