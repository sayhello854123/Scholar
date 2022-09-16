rm(list = ls())  
options(stringsAsFactors = F)#清空环境
library(VennDiagram)               #引用包
outFile="intersectGenes.txt"       #输出交集基因文件
outPic="venn.pdf"                  #输出图片文件
dir.create('10.VennPlot')
setwd('./10.VennPlot/')
files=dir()                        #获取目录下所有文件
files=grep("txt$",files,value=T)   #提取TXT结尾的文件
geneList=list()


#读取所有txt文件中的基因信息，保存到GENELIST
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)        #读取
  geneNames=as.vector(rt[,1])              #提取基因名
  geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
  uniqGene=unique(geneNames)               #基因取unique
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

#绘制vennͼ
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)),
                       col = "transparent" )
pdf(file=outPic, width=8.5, height=7)
grid.draw(venn.plot)
dev.off()

#保存交集基因
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
