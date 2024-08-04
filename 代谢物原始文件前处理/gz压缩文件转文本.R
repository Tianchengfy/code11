# 确保已加载所需的包
library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(readxl)
library(writexl)

# 设置工作路径
setwd(dir="E:/1400血液代谢物原始数据")  # 更改为你自己的工作路径

# 获取工作路径下所有的.gz文件
FileNames <- list.files(paste0(getwd()), pattern=".gz")
exp_dat_ids <- FileNames
exps <- FileNames

# 读取代谢物信息
dataa <- fread("1400代谢物信息.csv")
Phenotype <- dataa[,3]
Samplesize <- dataa[,4]

# 循环变量初始化
qaq <- 1

# 循环处理每个文件
for (qaq in 1:length(exp_dat_ids)) {
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  # 读取文件
  d1 <- try(fread(paste0(getwd(),"/",FileNames[qaq]), sep = "\t"), silent = T)
  
  # 不进行P值筛选，直接添加表型和样本大小信息
  d1$Phenotype <- Phenotype[qaq]
  d1$samplesize <- Samplesize[qaq]
  colnames(d1) <- c("chr","pos","effect_allele","other_allele","eaf","beta","se","pval","SNP","Phenotype","samplesize")
  
  # 使用 fwrite 写入 txt 文件
  # 使用 fwrite 写入 txt 文件
  fwrite(d1, file = paste0("F:/1400血液代谢物原始数据文本/", exp_dat_ids[qaq], ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
