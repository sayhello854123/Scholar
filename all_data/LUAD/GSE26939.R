library(tidyverse)
library(devtools)
library(AnnoProbe)
library(GEOquery)
library(limma)

# 设置连接大小
Sys.setenv("VROOM_CONNECTION_SIZE" = 172157430)

# 创建输出目录
dir.create("step00-data", showWarnings = FALSE)
dir.create("step01-normalize", showWarnings = FALSE)

# ===== 1. 下载数据 =====
cat("=== 下载GSE26939数据 ===\n")
gset <- getGEO('GSE26939', destdir=".",
               AnnotGPL = F,
               getGPL = F) 

a <- gset[[1]]
dat1 <- a@assayData$exprs
cat("原始数据维度:", dim(dat1), "\n")

# ===== 2. 下载并处理平台注释 =====
cat("\n=== 处理平台注释 GPL9053 ===\n")
gpl <- getGEO('GPL9053', destdir=".")

# 查看平台注释列名
cat("平台注释列名:\n")
print(colnames(Table(gpl)))

# 提取探针-基因对应关系
probe2gene <- Table(gpl)[, c("ID", "ORF")]
cat("\n探针-基因对应关系示例:\n")
print(head(probe2gene))

# 创建ID映射
ids <- probe2gene 
colnames(ids) <- c('probe_id', 'symbol')

# 转换为字符型确保匹配
ids$probe_id <- as.character(ids$probe_id)
rownames(dat1) <- as.character(rownames(dat1))

# 去除空symbol
ids <- ids[ids$symbol != '', ]
cat("\n去除空symbol后:", nrow(ids), "个探针\n")

# ===== 3. 匹配探针ID =====
# 检查ID格式
cat("\n=== 检查ID匹配 ===\n")
cat("表达矩阵探针ID示例:\n")
print(head(rownames(dat1)))
cat("\n平台注释探针ID示例:\n")
print(head(ids$probe_id))

# 尝试数字匹配（因为可能存在格式差异）
ids$probe_id_num <- as.numeric(ids$probe_id)
dat1_rownames_num <- as.numeric(rownames(dat1))

# 找到共同探针
common_probes_num <- intersect(ids$probe_id_num, dat1_rownames_num)
cat("\n匹配的探针数:", length(common_probes_num), "\n")

if(length(common_probes_num) < 100) {
  # 如果数字匹配失败，尝试直接匹配
  cat("数字匹配失败，尝试直接匹配...\n")
  ids <- ids[ids$probe_id %in% rownames(dat1), ]
  dat1 <- dat1[ids$probe_id, ]
} else {
  # 使用数字匹配
  cat("使用数字匹配\n")
  ids <- ids[ids$probe_id_num %in% common_probes_num, ]
  dat1 <- dat1[as.character(ids$probe_id_num), ]
}

cat("匹配后探针数:", nrow(dat1), "\n")

# ===== 4. 处理多探针对应同一基因 =====
cat("\n=== 处理重复基因 ===\n")
ids$median <- apply(dat1, 1, median) 
ids <- ids[order(ids$symbol, ids$median, decreasing = TRUE), ]
ids <- ids[!duplicated(ids$symbol), ]

# 提取唯一基因的表达数据
if("probe_id_num" %in% colnames(ids)) {
  dat1 <- dat1[as.character(ids$probe_id_num), ]
} else {
  dat1 <- dat1[ids$probe_id, ]
}

rownames(dat1) <- ids$symbol
rt1 <- dat1

cat("去重后基因数:", nrow(rt1), "\n")

# ===== 5. 提取临床信息 =====
cat("\n=== 提取临床信息 ===\n")
metadata <- pData(a)
clinical = data.frame(gsm=metadata[,2],
          futime= trimws(sapply(as.character(metadata$characteristics_ch1.9),function(x) strsplit(x,":")[[1]][2])),
          futime= trimws(sapply(as.character(metadata$characteristics_ch1.10),function(x) strsplit(x,":")[[1]][2])),
          futime= trimws(sapply(as.character(metadata$characteristics_ch1.11),function(x) strsplit(x,":")[[1]][2])),
          fustat=trimws(sapply(as.character(metadata$characteristics_ch1.12),function(x) strsplit(x,":")[[1]][2])),
          fustat=trimws(sapply(as.character(metadata$characteristics_ch1.13),function(x) strsplit(x,":")[[1]][2])),          
          Age=trimws(sapply(as.character(metadata$characteristics_ch1.1),function(x) strsplit(x,":")[[1]][2])),
          Sex=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2])),
          Stage =trimws(sapply(as.character(metadata$characteristics_ch1.5),function(x) strsplit(x,":")[[1]][2])),
          Grade =trimws(sapply(as.character(metadata$characteristics_ch1.6),function(x) strsplit(x,":")[[1]][2])))

write.csv(clinical,file = './step00-data/GSE26939_cli.csv')


# 保存原始数据
write.csv(rt1, file = './step00-data/GSE26939.csv')


# ===== 7. 数据标准化 =====
cat("\n=== 数据标准化 ===\n")
rt <- read.csv('./step00-data/GSE26939.csv', header = TRUE, row.names = 1)

# 判断是否需要log转换
qx <- as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))

if(LogC) {
  cat("需要log2转换\n")
  rt[rt < 0] <- 0
  rt <- log2(rt + 1)
} else {
  cat("数据已经是log转换\n")
}

# 标准化
data <- normalizeBetweenArrays(rt)
write.csv(data, file = './step00-data/GSE26939_normalize.csv')
