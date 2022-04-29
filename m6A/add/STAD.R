metabolize=read.table('./00.data/STAD/07.GSVA/ssGSEA.result.txt', header=T, sep="\t", check.names=F, row.names=1)
metabolize=t(metabolize)
m6A1 <- read.table('./00.data/STAD/04.meger/m6A_merge.txt',header = T,row.names = 1,sep = '\t',check.names = F)
m6A <- t(m6A1)
dir.create('./add/STAD')
dir <- './add/STAD/'

for (i in colnames(metabolize)) {
  a <- median(metabolize[,i])
  risk <- as.vector(ifelse(metabolize[,i] > a,"high","low"))
  table(risk)
  rt <- cbind(m6A,risk)
  risk="risk"
  dotCol=c("red","blue")
  outTab=data.frame()
  outTab1 <- data.frame()
  for(cell in colnames(rt[,1:(ncol(rt)-1)])){
    rt1=rbind(cell=as.numeric(rt[,cell]),risk=rt[,risk])
    rt1=data.frame(t(rt1))
    rt1$cell <- as.numeric(rt1$cell)
    wilcoxTest=wilcox.test(cell ~ risk, data=rt1)
    pValue=wilcoxTest$p.value
    sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
    outTab=rbind(outTab,cbind(Metabolize=i,gene=cell,pValue=pValue,sig))
    pval=0
    if(pValue<0.001){
      pval=signif(pValue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=sprintf("%0.3f",pValue)
    }
    
    if(pValue<0.05){	  
      b = boxplot(cell ~ risk, data = rt1,outline = FALSE, plot=F) 
      yMin=min(b$stats)
      yMax = max(b$stats/5+b$stats)
      ySeg = max(b$stats/10+b$stats)
      ySeg2 = max(b$stats/12+b$stats)
      n = ncol(b$stats)
      
      pdf(file=paste0(dir,i,'_',cell,".pdf"),width=6,height=5)
      par(mar = c(4,7,3,3))
      boxplot(cell ~ risk, data = rt1, ylab = cell,col=dotCol,xlab="",names=c("High expression ","Low expression"),
              cex.main=1.2, cex.lab=1, cex.axis=1, ylim=c(yMin,yMax), outline = FALSE)
      segments(1,ySeg, n,ySeg);segments(1,ySeg, 1,ySeg2);segments(n,ySeg, n,ySeg2)
      text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1,pos=3)
      title(main = i, sub = "sub title")
      dev.off()
    }
    outTab1 <- rbind(outTab1,outTab)
  }
  
}
