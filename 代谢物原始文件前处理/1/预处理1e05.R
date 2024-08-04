
setwd(dir="F:/测试") 


library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)


FileNames <-list.files(paste0(getwd()),pattern=".gz")
exp_dat_ids <- FileNames
exps <- FileNames

dataa <- fread("1400代谢物信息.csv")
Phenotype <- dataa[,3]
Samplesize <- dataa[,4]

qaq <- 1

library(data.table) 

for (qaq in 1:length(exp_dat_ids)) {
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  d1 <- try(fread(paste0(getwd(),"/",FileNames[qaq]), sep = "\t"), silent = T)
  d3 <- subset(d1, d1$p_value < 5e-6)
  d3$Phenotype <- Phenotype[qaq]
  d3$samplesize <- Samplesize[qaq]
  colnames(d3) <- c("chr","pos","effect_allele","other_allele","eaf","beta","se","pval","SNP","Phenotype","samplesize")
  
  
  fwrite(d3, file = paste0("deal/", exp_dat_ids[qaq], ".txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
