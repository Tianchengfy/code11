
setwd("E:/代谢物最新/中介分析/DVT/多变量")

library(tidyr)
library(data.table)

DVT<- fread("DVT整理.txt",header=T)


colnames(DVT)


DVT = DVT[,c(3,4,5,7,9,10,11,12,13,14)]

colnames(DVT) = c("other_allele","effect_allele","SNP","pval","beta","se","eaf","ncase", "ncontrol","samplesize")

DVT$pval <- as.numeric(DVT$pval)

DVT$beta <- as.numeric(DVT$beta)

DVT$se <- as.numeric(DVT$se)



#暴露1-PDGF-----

PDGF<- fread("meta_filt_PDGF_BB_noBMI.txt",header=T)

colnames(PDGF)##查看列名

##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"samplesize"

PDGF = PDGF[,c(1,2,3,4,8,9,10,16)]

colnames(PDGF) = c("SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize")

PDGF$pval <- as.numeric(PDGF$pval)

PDGF$beta <- as.numeric(PDGF$beta)

PDGF$se <- as.numeric(PDGF$se)



#暴露2-FL-----

FL<- fread("GCST90200155.txt",header=T)

colnames(FL)
FL$samplesize<-6448

##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"samplesize"

FL = FL[,c(3,4,5,6,7,8,9,10)]

colnames(FL) = c("effect_allele","other_allele","eaf","beta","se","pval","SNP","samplesize")

FL$pval <- as.numeric(FL$pval)

FL$beta <- as.numeric(FL$beta)

FL$se <- as.numeric(FL$se)


#暴露3-PH------

PH<- fread("GCST90200771.txt",header=T)

colnames(PH)##查看列名
PH$samplesize<-8209

##标准化列名："SNP"，"beta"，"se"，"eaf"，"effect_allele"，"other_allele"，"pval"，"samplesize"

PH = PH[,c(3,4,5,6,7,8,9,10)]

colnames(PH) = c("effect_allele","other_allele","eaf","beta","se","pval","SNP","samplesize")

PH$pval <- as.numeric(PH$pval)

PH$beta <- as.numeric(PH$beta)

PH$se <- as.numeric(PH$se)


library(TwoSampleMR)

library(ieugwasr)

#2.多个暴露分别提取强相关SNP+去除LD------------

#①暴露1的工具变量----------

X1<-subset(PDGF,pval<5e-6)

X1<-format_data(X1,type = "exposure")

X111<-clump_data(X1,clump_kb = 10000,
                 
                 clump_r2 = 0.001)#在线clump


#②暴露2的工具变量----------

X2<-subset(FL,pval<1e-5)

X2<-format_data(X2,type = "exposure")

X222<-clump_data(X2,clump_kb =10000,
                 
                 clump_r2 =0.001)#在线clump

#③暴露3的工具变量-----------

X3<-subset(PH,pval<1e-5)

X3<-format_data(X3,type = "exposure")

X333<-clump_data(X3,clump_kb =10000,
                 
                 clump_r2 =0.001)#在线clump


library(dplyr)

library(plyr)
colnames(X111)
colnames(X222)
colnames(X333)


SNP_list <- rbind.fill(X111, X222, X333)

SNP_list <- select(SNP_list,-id.exposure)



SNP_list_c<- clump_data(SNP_list)


#去除重复的SNP

SNP_u<-unique(SNP_list_c$SNP)


#X1--46

PDGF_X1<-PDGF[PDGF$SNP%in%SNP_u,]

#X2--52

FL_X2<-FL[FL$SNP%in%SNP_u,]

#X3--52

PH_X3<-PH[PH$SNP%in%SNP_u,]

#DVT--52

DVT_out<-DVT[DVT$SNP%in%SNP_u,]

#5.取交集【两两之间取交集，确定在所有GWAS数据中都有的SNP】--------

snps<-intersect(PDGF_X1$SNP,FL_X2$SNP)

snps<-intersect(snps,PH_X3$SNP)

snps<-intersect(snps,DVT_out$SNP)

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
  
  FL,
  
  type = "exposure",
  
  snps = snps)

#自定义

d2$exposure<-"FL"

d2$id.exposure<-"FL"

#③X3暴露--------

d3<-format_data(
  
  PH,
  
  type = "exposure",
  
  snps = snps)

#自定义

d3$exposure<-"PH"

d3$id.exposure<-"PH"

#④X1、X2、X3合并，作为多变量的暴露-----------------------
colnames(d1)
colnames(d2)
colnames(d3)


mv_exposures <- rbind(d1, d2, d3)

#④Y作为结局-------------

d4<-format_data(
  
  DVT,
  
  type = "outcome",
  
  snps = snps)

#自定义名字

d4$outcome<-"DVT"

d4$id.outcome<-"DVT"

#⑤等位基因对齐--harmonis暴露【mv_exposures】和结局【dh4】--------------

mv_DAT<-mv_harmonise_data(mv_exposures, d4)

head(mv_DAT)#去除回文SNP等
#7.MVMR分析【TwoSampleMR包】------------

library(TwoSampleMR)

res <-mv_multiple(mv_DAT)

res
write.csv(res, file="多变量结果.csv")
res_OR<-generate_odds_ratios(res$result)
write.csv(res_OR, file="多变量结果OR.csv")
res_OR

