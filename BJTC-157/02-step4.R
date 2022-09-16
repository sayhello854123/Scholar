library(VennDiagram)
outFile="intersectGenes.txt" 
files=dir('./02.bulk/venn/')    #获取目录下所有文件
files=grep("txt$",files,value=T)#提取TXT结尾的文件
dir <- './02.bulk/venn/'
files <- paste0(dir,files)
geneList=list()

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

venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file='./02.bulk/venn/venn.pdf', width=5, height=5)
grid.draw(venn.plot)
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file='./02.bulk/venn/intersect.txt',intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
