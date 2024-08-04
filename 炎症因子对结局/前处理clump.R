设置工作目录
setwd("E:/代谢物/significant/阳性代谢物对结局/新建文件夹")

# 加载必要的包
library(ieugwasr)
library(data.table)
ldl_gwas <- fread("meta_filt_PDGF_BB_noBMI.txt", sep = "\t", header = TRUE)
head(ldl_gwas)
ldl_gwas$'P-value' <- as.numeric(ldl_gwas$'P-value')

ldl_gwas$effect_allele_frequency <- as.numeric(ldl_gwas$effect_allele_frequency)

ldl_gwas$standard_error <- as.numeric(ldl_gwas$standard_error)

ldl_gwas<-subset(ldl_gwas,ldl_gwas$'P-value'<5e-06) 

ldl_iv <- ldl_gwas[, c("MarkerName", "P-value")] #选取rsid和p值这两列
colnames(ldl_iv) <- c("rsid", "pval" ) #修改列名，列名必须为rsid和pval

a <- ld_clump(
  #dat = X1,
  clump_kb = 10000,
  clump_r2 = 0.001,    #连锁不平衡的标准
  pop = "EUR",
  dplyr::tibble(rsid=ldl_gwas$MarkerName, pval=ldl_gwas$'p-value'),
  #plink.exe的位置
  plink_bin = "C:/Windows/本地plink/plink_win64_20230116/plink.exe",
  #EUR.bed的位置
  bfile = "C:/Windows/本地plink/EUR/EUR"
)


c <- ldl_gwas[which(ldl_gwas$MarkerName%in%a$rsid), ] 

head(c)



b1<-format_data(
  
  c,
  
  type = "exposure",
  
  snps = NULL,
  
  header = TRUE,
  
  snp_col = "MarkerName",
  
  beta_col = "Effect",
  
  se_col = "StdErr",
  
  eaf_col = "Freq1",
  
  effect_allele_col = "Allele1",
  
  other_allele_col = "Allele2",
  
  pval_col = "P-value",
  
  samplesize_col = "Total_N",
  
  log_pval = FALSE)




# 定义计算 PVE 的函数
PVEfx <- function(BETA, SE, N) {
  pve <- (BETA^2) / ((BETA^2) + ((SE^2) * N))
  return(pve)
}

# 应用 PVE 函数到数据框 b1 的相应列，并计算 R²
b1$PVE <- mapply(PVEfx, b1$beta.exposure, b1$se.exposure, N = b1$samplesize.exposure)

# 计算 F 统计量
b1$FSTAT <- ((b1$samplesize.exposure - 2) / 1) * (b1$PVE / (1 - b1$PVE))

# 输出数据框 b1 中的 R² 和 F 值
b1





write.csv(b1, file="PDGF.csv", row.names=F)

