
library(pacman)
library(TwoSampleMR)
library(data.table)
p_load(MendelianRandomization,purrr,readr,TwoSampleMR)
setwd("E:/代谢物最新/去混杂/clumped")
MRresult_data = read.table("1.txt",header = TRUE,sep = "\t")
immcell_ID = as.vector(MRresult_data$id)
immcell_ID=immcell_ID[1:3]
for (i in immcell_ID){
  setwd("E:/代谢物最新/去混杂/clumped")
  data<-fread(i)
  colnames(data)
  grp_size<-100
  grps<-split(data$SNP,ceiling(seq_along(data$SNP)/grp_size))
  result<-map_dfr(grps,~phenoscanner(snpquery=.x,
                                     catalogue = "GWAS",
                                     pvalue=1e-05,
                                     proxies = "None",
                                     r2=0.8,
                                     build = 38)$results)
  setwd("E:/代谢物最新/去混杂/clumped/out")
  write.table(result,i,sep = "\t",quote = F,row.names = F)
  
}
