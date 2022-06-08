data1 <- read.table('./00.data/01.TCGA/TCGA.TPM.txt',sep = '\t',header = T,check.names = F)
rt=as.matrix(data1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data <- cbind(gene_id = rownames(data), description = "na", data)
## 创建gct文件
gct_file = paste0("my.gct")
sink(gct_file)
cat("#1.2\n")
cat(paste0(nrow(data), "\t", (length(colnames(data)) - 2), "\n"))
sink()
write.table(data, gct_file, append = T, quote = F, row.names = F, col.names = T,
            sep = "\t")
cluster <-read.csv('./05.clicor/risk.csv', header=T, row.names=1)
group_list <- cluster[,9]
cls_file = paste0("my.cls")
sink(cls_file)
cat(paste0(length(group_list), "\t", length(unique(group_list)), "\t1\n"))
cat(paste0("#", paste(unique(group_list), collapse = "\t"), "\n"))
cat(paste(group_list, collapse = "\t"), "\n")
sink()


#引用包
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

gseaplot <- function(dir,files,x,outfile){
  files <- paste0(dir,files)
  data=lapply(files, read.delim)                             #读取每个文件
  names(data)=sapply(strsplit(files,'/'),'[',4)
  dataSet=ldply(data, data.frame)
  dataSet$pathway = gsub(".tsv", "", dataSet$.id)
  #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #处理后有73种差异还比较明显的颜色，基本够用
  #gseaCol = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
  
  gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
    geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
    labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
    geom_hline(yintercept = 0) + 
    guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
  pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
    scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
    labs(x = x, y = "", title = "") + 
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)
  gGsea = ggplot_gtable(ggplot_build(pGsea))
  gGene = ggplot_gtable(ggplot_build(pGene))
  maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
  gGsea$widths = as.list(maxWidth)
  gGene$widths = as.list(maxWidth)
  dev.off()
  pdf(file=outfile,    #输出图片的文件
      width=9,                    #设置输出图片高度
      height=5.5)                 #设置输出图片高度
  par(mar=c(5,5,2,5))
  grid.arrange(arrangeGrob(gGsea, gGene, nrow=2, heights=c(.8,.25)))
  dev.off()
}
dir <- './09.ConsensusClusterPlus/result1/'
files=grep(".tsv", dir('./09.ConsensusClusterPlus/result1/'), value=T) #获取目录下的所有tsv文件
x <- "A cluster score<----------->B cluster score"
outfile <- "./09.ConsensusClusterPlus/all1.pdf"
gseaplot(dir=dir,files=files,x=x,outfile=outfile)
