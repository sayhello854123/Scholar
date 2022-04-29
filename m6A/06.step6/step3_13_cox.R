train=read.table("./00.data/ESCA/11.Model/m6Ascore_train_group.txt", header=T, sep="\t", check.names=F, row.names=1)
train <- train[,c(1,2,(length(colnames(train))-1))]
colnames(train) = c("futime", "fustat","m6A")
test=read.table("./00.data/ESCA/11.Model/m6Ascore_test_group.txt", header=T, sep="\t", check.names=F, row.names=1)
test <- test[,c(1,2,(length(colnames(test))-1))]
colnames(test) = c("futime", "fustat","m6A")
score <- rbind(train,test)
metabolize=read.table('./00.data/ESCA/07.GSVA/ssGSEA.result.txt', header=T, sep="\t", check.names=F, row.names=1)
metabolize=t(metabolize)
sameSample=intersect(row.names(score), row.names(metabolize))
data=cbind(score[sameSample,,drop=F], metabolize[sameSample,,drop=F])

rt <- data
outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  sigGenes=c(sigGenes,i)
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
rownames(outTab) <- outTab[,1]
outTab <- outTab[,-1]
sigGeneExp=rbind(id=colnames(outTab), outTab)
write.table(sigGeneExp, file="./03.step3/10.cor/ESCA_cox.txt", sep="\t", quote=F, col.names=F)

rt=read.table('./03.step3/10.cor/ESCA_cox.txt',header=T,sep="\t",row.names=1,check.names=F)
gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"HR.95L")
hrHigh=sprintf("%.3f",rt$"HR.95H")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file='./03.step3/10.cor/cox.pdf', width = 8.5, height =5.5)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.9-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.9-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()
