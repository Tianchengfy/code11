#注意事项如下
#该代码要修改3个工作目录，其他内容不要改
#finnen_R9.txt为参考结局文件，需要更换为自己的数据，注意修改自己结局数据的列名（第49-64行）
#请看参考文件对数据处理的描述

setwd("E:\\41种炎症因子最新\\血栓\\正向\\DVT")#修改工作目录1

library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(openxlsx)
name="exposure_data"

FileNames <-list.files(paste0(getwd(),"/",name),pattern=".txt")
exp_dat_ids <- FileNames
exps <- FileNames
# 结局
outNames <-list.files("E:\\41种炎症因子最新\\血栓\\正向\\DVT",pattern=".txt")#修改工作目录2

#for (ID in IDs) 
for (outname in outNames) {
  ID=gsub(".txt","",outname)
  out<-fread(paste0("E:\\41种炎症因子最新\\血栓\\正向\\DVT\\",outname),header = T)##修改工作目录1
  
  out$trait <- 'outcome'  
  outcomeid <- out
  rm(out)
  dir.create(path = ID)  
  #######以下为循环代码，不需要进行更改######################
  for (qaq in 1:length(exp_dat_ids)) {{ 
    exp_dat_id <- exp_dat_ids[qaq]
    expname <- exps[qaq]
    exp_data<-read_exposure_data(filename = paste0('./exposure_data/',exps[qaq]),
                                 sep = "\t",
                                 snp_col = "MarkerName",
                                 beta_col = "Effect",
                                 se_col = "StdErr",
                                 effect_allele_col = "Allele1",
                                 other_allele_col = "Allele2",
                                 eaf_col = "Freq1",
                                 pval_col = "P-value",
                                 samplesize_col = "Total_N",
                                 phenotype_col = "PHENO",
    )
    
    
    if(length(exp_data[,1])>2){
      outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "SNP")
      write.csv(outcome_dat,file = "d.csv")
      #注意修改自己结局数据的列名
      out_data <- read_outcome_data(
        snps = exp_data$SNP, 
        filename = "d.csv",
        sep = ",",
        snp_col = "SNP",  
        beta_col = "beta",  
        se_col = "sebeta",  
        effect_allele_col = "alt",  
        other_allele_col = "ref",  
        eaf_col = "eaf.exposure", 
        pval_col = "pval", 
        samplesize_col = "samplesize",
        chr_col="#chrom", pos_col = "pos")
      
      
      dat <- TwoSampleMR::harmonise_data(
        exposure_dat = exp_data,
        outcome_dat = out_data)
      
      ####回文的直接去除
      dat <-subset(dat,palindromic==FALSE)
      
      #计算F值和R2
      get_f<-function(dat,F_value=10){
        log<-is.na(dat$eaf.exposure)
        log<-unique(log)
        if(length(log)==1)
        {if(log==TRUE){
          print("数据不包含eaf，无法计算F统计量")
          return(dat)}
        }
        if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
          return(dat)}
        if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
          return(dat)}
        if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
          return(dat)}
        
        
        if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
          R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
          F<- (dat$samplesize.exposure-2)*R2/(1-R2)
          dat$R2<-R2
          dat$F<-F
          dat<-subset(dat,F>F_value)
          return(dat)
        }
      }
      
      dat <- get_f(dat, F_value = 10)
      
      res=TwoSampleMR::mr(dat)
      
      print(paste0(expname,"_SNP数_",res$nsnp[1]))
      
      results <- TwoSampleMR::generate_odds_ratios(res)
      
      results$estimate <- paste0(
        format(round(results$or, 2), nsmall = 2), " (", 
        format(round(results$or_lci95, 2), nsmall = 2), "-",
        format(round(results$or_uci95, 2), nsmall = 2), ")")
      resdata <- dat
      openxlsx::write.xlsx(dat,file = paste0(ID,"/",expname,"-dat.xlsx"), row.names = FALSE)
      
      names(resdata)
      Assumption13 <- subset(resdata,mr_keep==TRUE,
                             select = c("SNP","pval.exposure",
                                        "pval.outcome", # "F_statistic",
                                        "mr_keep"))
      
      res_hete <- TwoSampleMR::mr_heterogeneity(dat)
      res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
      res_leaveone <- mr_leaveoneout(dat)  # 
      
      #####MR-PRESSO######
      # PRESSO
      mr_Presso<-function(dat,num=10000){
        library(TwoSampleMR)
        library(MRPRESSO)
        library(dplyr)
        set.seed(123)
        try (mr_presso_res<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat,  
                                      SignifThreshold = 0.05,NbDistribution = num))
        return(mr_presso_res)
      }
      mr_presso_pval<-function(mr_presso_res){ 
        try ( mr_presso_main<-mr_presso_res$`Main MR results`)
        try ( mr_presso_main[3:5,]<-NA) 
        return(mr_presso_main)
      }
      mr_presso_snp<-function(mr_presso_res,mr_presso_main,dat,type="list"){
        data_re<-list()
        if(type=="list"){
          for(i in 1:length(mr_presso_res)){
            res<-mr_presso_res[[i]]
            main<-mr_presso_main[[i]]
            data<-dat[[i]]
            try(if(is.na(main[2,6])==FALSE){
              outliers<-which(res$`Outlier Test`$Pvalue<0.05)
              data$mr_keep[outliers]<-FALSE
            })
            data_re[[i]]<-data
            names(data_re)[[i]]<-names(dat)[[i]]
          }
          return(data_re)
        }
        
        if(type=="data"){
          res<-mr_presso_res$`MR-PRESSO results`
          main<-mr_presso_main
          data<-dat
          try(if(is.na(main[2,6])==FALSE){
            outliers<-which(res$`Outlier Test`$Pvalue<0.05)
            data$mr_keep[outliers]<-FALSE
          })
          return(data)
        }
      }
      
      #MR-PRESSO （寻找离群值）
      mr_presso_res <- mr_Presso(dat, num = 1000)
      mr_presso_main <- mr_presso_pval(mr_presso_res)
      dat <- mr_presso_snp(mr_presso_res, mr_presso_main, dat, type = "data")
      
      ###根据beta、p计算se
      resMRPRESSO=mr_presso_res[["Main MR results"]]
      resMRPRESSO
      global_test_p <- mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
      se1=sqrt(((resMRPRESSO[1,3])^2)/qchisq(Pval_raw <- resMRPRESSO[1,6],1,lower.tail=F))
      se2=sqrt(((beta_cor <- resMRPRESSO[2,3])^2)/qchisq(Pval_cor <- resMRPRESSO[2,6],1,lower.tail=F))
      resMRPRESSO <- resMRPRESSO %>%
        dplyr::mutate(se = c(se1,se2))
      
      #生成globe test_P和outliers
      outliers <- dat$SNP[mr_presso_res[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]]
      outliers = as.data.frame(outliers)
      global_test_p = as.data.frame(global_test_p)
      resMRPRESSO
      TTT <- plyr::rbind.fill(resMRPRESSO, global_test_p,outliers)
      TTT <- as.data.frame(TTT)
      TTT
      openxlsx::write.xlsx(TTT,file = paste0(ID,"/",expname,"-MR-PRESSO.xlsx"), row.names = FALSE)
      
      #########
      
      openxlsx::write.xlsx(x = list(
        "main"=results,
        "Assumption13"=Assumption13,
        "pleiotropy"=res_plei,
        "heterogeneity"=res_hete,
        "leaveone"=res_leaveone),
        overwrite = TRUE,
        paste0(ID,"/",expname,"-res.xlsx"))
      
      
      p1 <- mr_scatter_plot(res, dat)
      p1[[1]]
      pdf(paste0(ID,"/",expname,"_scatter.pdf"))
      print(p1[[1]])
      dev.off()
      
      res_single <- mr_singlesnp(dat)
      p2 <- mr_forest_plot(res_single)
      pdf(paste0(ID,"/",expname,"_forest.pdf"))
      print(p2[[1]])
      dev.off()
      
      p3 <- mr_funnel_plot(res_single)
      pdf(paste0(ID,"/",expname,"_funnel.pdf"))
      print(p3[[1]])
      dev.off()
      
      res_loo <- mr_leaveoneout(dat)
      pdf(paste0(ID,"/",expname,"_leave_one_out.pdf"))
      print(mr_leaveoneout_plot(res_loo))
      dev.off()
    }
    
    library(magrittr)
    res2 <- res[1:3,]
    
    judge_1 <- function(mr_res=res2) {
      mr_res$b_direction <- as.numeric(sign(mr_res$b))
      mr_res$b_direction=ifelse(abs(sum(mr_res$b_direction))==3 ,
                                NA,"Inconsistent direction")
      mr_res$p_no <- NA
      mr_res[mr_res$method=="MR Egger","p_no"] <- ifelse(
        mr_res[mr_res$method=="MR Egger","pval"]<0.05," ",
        "MR Egger")
      mr_res[mr_res$method=="Weighted median","p_no"] <- ifelse(
        mr_res[mr_res$method=="Weighted median","pval"]<0.05," ",
        "Weighted median")
      mr_res[mr_res$method=="Inverse variance weighted","p_no"] <- ifelse(
        mr_res[mr_res$method=="Inverse variance weighted","pval"]<0.05,
        " ","Inverse variance weighted")
      mr_res$p_no <- paste(mr_res$p_no,collapse = " ")
      mr_res$p_no=trimws(mr_res$p_no,which = c("both"))
      return(mr_res)
    }
    
    res3 <- judge_1(mr_res = res2)
    
    # 转换成论文格式
    library(magrittr)
    # Main result 
    res4 <- tidyr::pivot_wider(
      res3,names_from ="method",names_vary = "slowest",
      values_from = c("b","se","pval") )
    # Heterogeneity statistics
    res_hete2 <- tidyr::pivot_wider(
      res_hete,names_from ="method",names_vary = "slowest",
      values_from = c("Q","Q_df","Q_pval") ) %>% 
      dplyr::select( -id.exposure,-id.outcome,-outcome,-exposure)
    # Horizontal pleiotropy
    res_plei2 <- dplyr::select(res_plei,
                               egger_intercept,se,pval)
    # Merge 
    res_ALL <- cbind(res4, res_hete2, res_plei2)
    
    write.csv(res_ALL,file = paste0(ID,"/",expname,".csv"), row.names = FALSE)
  }}
  
  #把导出的路径改成自己的工作路径  
  fs=list.files(paste0("./",ID), pattern = "csv",full.names = TRUE) #把导出的路径改成自己的工作路径  
  df = map_dfr(fs, read.csv)
  write.csv(df,paste0("res_",ID,".csv"))
}
