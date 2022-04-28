library(limma)
library(parallel)

###Obtaining TCGA-LIHC expression data
rt<-data.table::fread(file ='./TCGA-LIHC.htseq_fpkm.tsv',data.table = F)
gencode <- data.table::fread('./gencode.v22.annotation.gene.probeMap',data.table = F)
rownames(gencode) <- gencode[,1]
rownames(rt) <- rt[,1]
rt <- rt[,-1]
data <- 2^rt-1
same <- intersect(row.names(gencode),row.names(data))
length(same)
data1 <- cbind(gencode[same,],data[same,])
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
##筛选肿瘤样本
group=sapply(strsplit(colnames(tpm),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
tumor=tpm[,group==0]
####save tumor
save(tumor,file = './expFile.RData')

## Use Cancer Genome Project (CGP) data
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)##Get built-in datasets
length(unique(drugData2016$Drug.name))
drug2016 <- unique(drugData2016$Drug.name)###Get drug name
###Create Folder
dir.create('./IC50_result')
dir.create('./result_file')
system.time({
  detectCores()##查看有多少线程
  cl <- makeCluster(50)##根据自己的线程进行设定
  myresult <- parLapply(cl,drug2016,
                        function(x){
                          library(pRRophetic)
                          library(ggpubr)
                          library(ggplot2)
                          library(rstatix)
                          set.seed(1248203)##设置随机数
                          dir <- './IC50_result/'
                          pir <- './result_file/'
                          load('./expFile.RData')
                          ##Predicting drug sensitivity
                          senstivity=pRRopheticPredict(testMatrix = tumor, 
                                                       drug=x, 
                                                       selection=1,
                                                       tissueType = "all",
                                                       batchCorrect = 'eb',
                                                       dataset = "cgp2016")##使用最新的数据集分析
                          senstivity=senstivity[senstivity!="NaN"]
                          ##Get TCGA Risk Data
                          risk <- read.table('./tcgaRisk.txt',header = T,sep = '\t',
                                             row.names = 1,check.names = F)
                          ##Combination of risk group and drug sensitivity
                          sameSample=intersect(row.names(risk), names(senstivity))
                          Risk=risk[sameSample, "risk",drop=F]
                          senstivity=senstivity[sameSample]
                          rt=cbind(Risk, senstivity)
                          
                          ##Set comparison groups
                          rt$risk=factor(rt$risk, levels=c("low", "high"))
                          type=levels(factor(rt[,"risk"]))
                          comp=combn(type, 2)
                          my_comparisons=list()##多重分组比较
                          for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
                          
                          #获取高低风险组差异pvalue
                          test=wilcox.test(senstivity~risk, data=rt)#使用wilcox.test
                          
                          #Retention of standard compliant results
                          if(test$p.value<0.05){
                            #绘制箱线图
                            boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
                                              xlab="Risk",
                                              ylab=paste0(x, " senstivity (IC50)"),##y轴名称
                                              legend.title="Risk",
                                              palette=c("#33FF66","#FF3300")##修改颜色
                            )+ 
                              stat_compare_means(comparisons=my_comparisons)
                            ggsave(filename = paste0(dir,"durgSenstivity - ", x, ".pdf") ,
                                   boxplot,width = 5,height = 4.5,dpi = 1000)
                            colnames(rt) <- c('risk',x) 
                            write.csv(rt,file = paste0(pir,x,'.csv'))
                          }

                          return(senstivity)
                        })
   stopCluster(cl)
})