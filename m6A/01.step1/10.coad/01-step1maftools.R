rm(list = ls())
require(maftools)
options(stringsAsFactors = F) 
library(data.table)
picDir='./01.step1/01.maftools'
dir.create(picDir)
geneRT=read.table("./00.data/m6A_gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)
maf_plot <- function(x,y){
tmp <- fread(x)
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
maf = read.maf(maf = tmp,
                     vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))
filename <- paste0(y,'-','oncoplot',".pdf",sep="")
outfile=paste(picDir,filename,sep="/")
pdf(file=outfile, width=6.5, height=6)
oncoplot(maf=maf,genes = gene,draw_titv=T)
dev.off()
}
maf_plot('./00.data/COAD/TCGA-COAD.mutect2_snv.tsv.gz','COAD')
maf_plot('./00.data/ESCA/TCGA-ESCA.mutect2_snv.tsv.gz','ESCA')
maf_plot('./00.data/LIHC/TCGA-LIHC.mutect2_snv.tsv.gz','LIHC')
maf_plot('./00.data/PAAD/TCGA-PAAD.mutect2_snv.tsv.gz','PAAD')
maf_plot('./00.data/STAD/TCGA-STAD.mutect2_snv.tsv.gz','STAD')
maf_plot('./00.data/READ/TCGA-READ.mutect2_snv.tsv.gz','READ')

library(ggimage)
library(grid)
library(ggplot2)
library(patchwork)
library(EBImage)
library(imager)
library("jpeg")
library(ggpubr)

A1 <- readJPEG('./02.maftools/COAD-oncoplot.jpg')
p0<-ggplot()+
  background_image(A1)+
  theme_void()
A2 <- readJPEG('02.maftools/ESCA-oncoplot.jpg')
p1<-ggplot()+
  background_image(A2)+
  theme_void()
A3 <- readJPEG('02.maftools/LIHC-oncoplot.jpg')
p2<-ggplot()+
  background_image(A3)+
  theme_void()
A4 <- readJPEG('02.maftools/PAAD-oncoplot.jpg')
p3<-ggplot()+
  background_image(A4)+
  theme_void()
A5 <- readJPEG('./02.maftools/STAD-oncoplot.jpg')
p4<-ggplot()+
  background_image(A5)+
  theme_void()
A6 <- readJPEG('./02.maftools/all.jpg')
p5<-ggplot()+
  background_image(A6)+
  theme_void()
p<- (p0|p1)/(p2|p3)/(p4|p5) 
p <- p+plot_annotation(tag_levels = 'A')
ggsave(filename="./02.maftools/ALL_maf.jpg",
       p,
       width=8,
       heigh=8,
       dpi = 1000)