setwd("F:/中介分析/代谢物对PDGF/多变量/VTE")

##1.暴露和结局--本地数据处理----------
library(TwoSampleMR)
library(tidyr)
library(data.table)
#结局AF-----

VTE<- fread("VTE整理.txt",header=T)


##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"samplesize"

##https://cvd.hugeamp.org/downloads.html

##29892015

##European only + Readme

colnames(VTE)##查看列名
##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"N"

VTE = VTE[,c(1,2,3,4,5,7,9,10,11,14)]

colnames(VTE) = c("chr","pos","other_allele","effect_allele","SNP","pval","beta","se","eaf","samplesize")

VTE$pval <- as.numeric(VTE$pval)

VTE$beta <- as.numeric(VTE$beta)

VTE$se <- as.numeric(VTE$se)

#暴露1------

PDGF<- fread("meta_filt_PDGF_BB_noBMI.txt",header=T)
##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"samplesize"
head(PDGF)
PDGF = PDGF[,c(1,2,3,4,8,9,10,16)]

colnames(PDGF) = c("SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize")

PDGF$pval <- as.numeric(PDGF$pval)

PDGF$beta <- as.numeric(PDGF$beta)

PDGF$se <- as.numeric(PDGF$se)
#暴露2-HDL-----
rm(ME)
ME<- fread("GCST90200653_buildGRCh38.tsv.gz.txt",header=T)

colnames(ME)##查看列名

##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"samplesize"

ME = ME[,c(1,2,3,4,5,6,7,8,9,10,11)]

colnames(ME) = c("chr","pos","effect_allele","other_allele","eaf","beta","se","pval","SNP","Phenotype","samplesize")

ME$pval <- as.numeric(ME$pval)

ME$beta <- as.numeric(ME$beta)

ME$se <- as.numeric(ME$se)
#save(HDL,file="HDL.rda")

#write.csv(HDL,"HDL.csv")


#读取暴露+结局数据---------------

load("VTE.rda")

load("PDGF.rda")

load("BG.rda")


library(TwoSampleMR)

library(ieugwasr)

#2.多个暴露分别提取强相关SNP+去除LD------------

#①暴露1的工具变量----------

X1<-subset(PDGF,pval<5e-6)

X1<-format_data(X1,type = "exposure")

X111<-clump_data(X1,clump_kb = 10000,
                 
                 clump_r2 = 0.001)#在线clump

###本地clump开始

#安装并加载包

#devtools::install_github("explodecomputer/plinkbinr")

library(plinkbinr)

get_plink_exe() #######此行代码运行之后会显示--[1]"plink.exe的路径"

#devtools::install_github("mrcieu/ieugwasr")

library(ieugwasr)

#运行本地clump代码

X111<- ld_clump(
  
  clump_kb = 10000,
  
  clump_r2 = 0.001,
  
  clump_p = 0.99,
  
  pop ="EUR",
  
  dplyr::tibble(rsid=X1$SNP, pval=X1$pval.exposure, id=X1$id.exposure),
  
  plink_bin = "C:/Windows/本地plink/plink_win64_20230116/plink.exe",
  bfile = "C:/Windows/本地plink/EUR/EUR"
  
)

#索引完整信息

##%in%是啥意思【a1 %in% a2表示的是向量a1中，同时出现在向量a2中的元素】

X111<-subset(X1,SNP %in% X111$rsid)

###本地clump结束

save(X111,file="X111.rda")

#②暴露2的工具变量----------

X2<-subset(ME,pval<1e-5)

X2<-format_data(X2,type = "exposure")

X222<-clump_data(X2,clump_kb =10000,
                 
                 clump_r2 =0.001)#在线clump

###本地clump开始

X222<- ld_clump(
  
  clump_kb = 10000,
  
  clump_r2 = 0.001,
  
  clump_p = 0.99,
  
  pop ="EUR",
  
  dplyr::tibble(rsid=X2$SNP, pval=X2$pval.exposure, id=X2$id.exposure),
  
  plink_bin = "C:/Windows/本地plink/plink_win64_20230116/plink.exe",
  bfile = "C:/Windows/本地plink/EUR/EUR"
  
)

#索引完整信息

##%in%是啥意思【a1 %in% a2表示的是向量a1中，同时出现在向量a2中的元素】

X222<-subset(X2,SNP %in% X222$rsid)

###本地clump结束

save(X222,file="X222.rda")


#3.合并且去LD+去重----------

#SNP合并

library(dplyr)

library(plyr)

colnames(X111)
colnames(X222)
X222 <- X222[,-c(1, 2)]
SNP_list <- rbind.fill(X111, X222)#按行合并，填充缺失的列

SNP_list <- select(SNP_list,-id.exposure)

#去LD

SNP_list_c<- clump_data(SNP_list,clump_kb = 10000,
                        
                        clump_r2 =0.001)

###本地clump开始

SNP_list_c<- ld_clump(
  
  clump_kb = 10000,
  
  clump_r2 = 0.001,
  
  clump_p = 0.99,
  
  pop ="EUR",
  
  dplyr::tibble(rsid=SNP_list$SNP, pval=SNP_list$pval.exposure, id=SNP_list$id.exposure),
  
  plink_bin = "C:/Windows/本地plink/plink_win64_20230116/plink.exe",
  bfile = "C:/Windows/本地plink/EUR/EUR"
  
)

#索引完整信息

##%in%是啥意思【a1 %in% a2表示的是向量a1中，同时出现在向量a2中的元素】

SNP_list_c<-subset(SNP_list,SNP %in% SNP_list_c$rsid)

###本地clump结束

save(SNP_list_c,file="SNP_list_c.rda")

#去除重复的SNP

SNP_u<-unique(SNP_list_c$SNP)

#4.提取SNP在所有暴露、结局中的信息-----------

#确定用哪些SNP

#X1--140

PDGF_X1<-PDGF[PDGF$SNP%in%SNP_u,]

#X2--140

ME_X2<-ME[ME$SNP%in%SNP_u,]


#AF--140

VTE_out<-VTE[VTE$SNP%in%SNP_u,]

#5.取交集【两两之间取交集，确定在所有GWAS数据中都有的SNP】--------

snps <- intersect(PDGF_X1$SNP, ME_X2$SNP)
snps <- intersect(snps, VTE_out$SNP)

#6.等位基因对齐-----------

#①X1暴露----------

d1<-format_data(
  
  PDGF,
  
  type = "exposure",
  
  snps = snps)

#自定义一下exposure+id

d1$exposure<-"PDGF"

d1$id.exposure<-"PDGF"

#②X2暴露----------

d2<-format_data(
  
  ME,
  
  type = "exposure",
  
  snps = snps)

#自定义

d2$exposure<-"ME"

d2$id.exposure<-"ME"


#④X1、X2、X3合并，作为多变量的暴露-----------------------
colnames(d1)
colnames(d2)
d2 <- d2[,-c(1, 2)]
mv_exposures <- rbind(d1, d2)

#④Y作为结局-------------

d4<-format_data(
  
  VTE,
  
  type = "outcome",
  
  snps = snps)

#自定义名字

d4$outcome<-"VTE"

d4$id.outcome<-"VTE"

#⑤等位基因对齐--harmonis暴露【mv_exposures】和结局【dh4】--------------

mv_DAT<-mv_harmonise_data(mv_exposures, d4)

head(mv_DAT)#去除回文SNP等

#7.MVMR分析【TwoSampleMR包】------------

library(TwoSampleMR)

res <-mv_multiple(mv_DAT)

res
write.csv(res, file="多变量结果200653.csv")
res_OR<-generate_odds_ratios(res$result)
write.csv(res_OR, file="多变量结果OR200653.csv")
res_OR

#8.将mv_DAT矩阵转换为environment--data----------

data<-cbind(mv_DAT[["outcome_beta"]],
            
            mv_DAT[["exposure_beta"]][,1],
            
            mv_DAT[["exposure_beta"]][,2],
            
            
            
            mv_DAT[["exposure_se"]][,1],
            
            mv_DAT[["exposure_se"]][,2],
            
            
            
            mv_DAT[["outcome_se"]])

data<-data.frame(data)

colnames(data) = c("betaY","betaX1","betaX2",
                   
                   "seX1","seX2","seY")

#9.MVMR分析【MendelianRandomization包--IVW+Egger+LASSO+WM+robust异质性检验(前两种)------------

#安装MendelianRandomization包

#if (!requireNamespace("MendelianRandomization"))

#install.packages("MendelianRandomization")

library(MendelianRandomization)

#准备此包函数专用输入数据

MRMVInputObject<- mr_mvinput(bx = cbind(data$betaX1,data$betaX2),
                             
                             bxse = cbind(data$seX1,data$seX2,data$seX3),
                             
                             by = data$betaY,
                             
                             byse = data$seY)

MRMVInputObject

#IVW方法

IVW <- mr_mvivw(MRMVInputObject,
                
                model = "default",
                
                correl = FALSE,
                
                distribution = "normal",
                
                alpha = 0.05)

IVW

##缺失OR和OR的95%CI怎么办--使用transform函数

transform<- function(beta,se) {
  
  OR= exp(beta)
  
  LL_OR= exp(beta-se*1.96)
  
  UL_OR=exp(beta+se*1.96)
  
  pval<-2*pnorm(abs(beta/se),lower.tail=FALSE)
  
  OR_P<-cbind(OR,LL_OR,UL_OR,pval)
  
  print(OR_P)
  
}

OR_ivw<- transform(IVW@Estimate,IVW@StdError)

#egger方法

egger<-mr_mvegger(MRMVInputObject,
                  
                  orientate = 1,
                  
                  correl = FALSE,
                  
                  distribution = "normal",
                  
                  alpha = 0.05)

egger

#LASSO

Lasso<-mr_mvlasso(
  
  MRMVInputObject,
  
  orientate = 1,
  
  distribution = "normal",
  
  alpha = 0.05,
  
  lambda = numeric(0)
  
)

Lasso

#median

median<-mr_mvmedian(
  
  MRMVInputObject,
  
  distribution = "normal",
  
  alpha = 0.05,
  
  iterations = 10000,
  
  seed = 314159265
  
)

median

#robust

#无需修改，先定义robust函数[34342032]

library(tidyverse)

library(RColorBrewer)

library(cowplot)

library(gridExtra)

library(robustbase)

mvmr_robust = function(bx, by, seby, k.max = 500, maxit.scale = 500){
  
  robmod = lmrob(by ~ bx - 1, weights = seby^-2, k.max = k.max,
                 
                 maxit.scale = maxit.scale)
  
  coefficients = summary(robmod)$coef[, 1]
  
  se = summary(robmod)$coef[, 2] / min(summary(robmod)$sigma, 1)
  
  return(list("coefficients" = coefficients, "se" = se))
  
}

#获取数据

mvmr_robust <- mvmr_robust(cbind(data$betaX1,data$betaX2,data$betaX3),
                           
                           data$betaY,data$seY,
                           
                           k.max = 1000, maxit.scale = 1000)

mvmr_robust

#函数--根据beta、se计算其余值

transform<- function(beta,se) {
  
  OR= exp(beta)
  
  LL_OR= exp(beta-se*1.96)
  
  UL_OR=exp(beta+se*1.96)
  
  pval<-2*pnorm(abs(beta/se),lower.tail=FALSE)
  
  OR_P<-cbind(OR,LL_OR,UL_OR,pval)
  
  print(OR_P)
  
}

mvmr_robust<- transform(mvmr_robust$coefficients,mvmr_robust$se)

#10.MVMR分析【MVMR包--条件F统计量+异质性检验+IVW】---------------

#install.packages("remotes")

#library(remotes)

#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

library(MVMR)

head(data)

F.data <- format_mvmr(BXGs = data[,c(2,3,4)],
                      
                      BYG = data[,1],
                      
                      seBXGs = data[,c(5,6,7)],
                      
                      seBYG = data[,8],
                      
                      RSID="NULL")

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)