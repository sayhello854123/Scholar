rm(list = ls())
options(stringsAsFactors = F) 
gc()
require(maftools)
options(stringsAsFactors = F) 
library(data.table)
picDir <- './01.maftools/'
geneRT=read.table("./00.data/m5c_gene.txt", header=T, sep="\t", check.names=F, row.names=1)
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
maf_plot('./00.data/TCGA-LIHC.mutect2_snv.tsv.gz','LIHC')
