#km分析
data=read.csv(file = './3.csv',header = T,row.names = 1)
colnames(data)=c("futime", "fustat", "gene","group")
#比较高低表达生存差异
diff=survdiff(Surv(futime, fustat) ~group,data =data)
pValue=1-pchisq(diff$chisq, df=1)
#对pvalue<0.05的基因绘制生存曲线
if(pValue<0.05){
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }

  fit <- survfit(Surv(futime, fustat) ~ group, data = data)
  #绘制生存曲线
p1 <-   ggsurvplot(fit,
                     data = data,
                     pval=pValue,
                     pval.size=6,
                     legend.title='NOP2',
                     legend.labs=c("High Expression","Low Expression"),
                     xlab="Time(years)",
                     ylab="Overall survival",
                     palette=c("red", "blue"),
                     break.time.by=1,
                     #conf.int=T,
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
  
pdf(file = './03.survial/NOP2.pdf',width = 8,height = 6)
print(p1)
dev.off()