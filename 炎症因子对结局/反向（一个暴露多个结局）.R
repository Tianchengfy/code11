library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)
setwd("E:\\41种炎症因子最新\\血栓\\反向\\PE")
exp_data<-read_exposure_data(
  filename = "PE结局整理.csv",
  sep = ",",
  phenotype_col = "trait",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure")

setwd("./41种炎症因子原始数据")


FileNames <-list.files(paste0(getwd()),pattern=".txt")
outcomeids <- FileNames
out_comes <- FileNames


dir.create(path = "final result")  



for (qaq in 1:length(outcomeids)) { 
  outcomeid <- outcomeids[qaq]
  out_come <- out_comes[qaq]
  
  
  if(length(exp_data[,1])>2){
    
    out<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = T)
    out$PHENO <- FileNames[qaq]
    outcomeid <- out
    head(out)
    outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "MarkerName")
    head(out)
    out_data<-format_data(outcome_dat,
                          type="outcome",
                          snp_col = "SNP",
                          samplesize_col = "Total_N",
                          beta_col = "Effect",
                          se_col = "StdErr",
                          pval_col = "P-value",
                          effect_allele_col = "Allele1",
                          other_allele_col = "Allele2",
                          eaf_col = "Freq1")
    
    
    
    
    
    if(length(out_data[,1])>0){  
      dat <- TwoSampleMR::harmonise_data(
        exposure_dat = exp_data,
        outcome_dat = out_data)
      
     
      dat <-subset(dat,palindromic==FALSE)
      
      
      
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
        if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize，无法计算F统计量")
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
      
      
      res=TwoSampleMR::mr(dat,method_list= c("mr_ivw" ,
                                             "mr_weighted_median" ,
                                             "mr_egger_regression",
                                             "mr_simple_mode",
                                             "mr_weighted_mode",
                                             "mr_wald_ratio"))
      
      print(paste0(out_come,"_SNP数_",res$nsnp[1]))
      
      results <- TwoSampleMR::generate_odds_ratios(res)
      
      results$estimate <- paste0(
        format(round(results$or, 2), nsmall = 2), " (", 
        format(round(results$or_lci95, 2), nsmall = 2), "-",
        format(round(results$or_uci95, 2), nsmall = 2), ")")
      resdata <- dat
      openxlsx::write.xlsx(dat,file = paste0("final result/",out_come,"-dat.xlsx"), row.names = FALSE)
      
      names(resdata)
      Assumption13 <- subset(resdata,mr_keep==TRUE,
                             select = c("SNP","pval.exposure",
                                        "pval.outcome", 
                                        "mr_keep"))
      
      openxlsx::write.xlsx(x = list(
        "main"=results,
        "Assumption13"=Assumption13),
        overwrite = TRUE,
        paste0("final result/",out_come,"-res.xlsx"))
      
    }}
  if(length(dat[,1])>2){
    res_hete <- TwoSampleMR::mr_heterogeneity(dat)
    res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    res_leaveone <- mr_leaveoneout(dat)  # 
    
   
    dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                      dat$se.exposure,
                                      dat$samplesize.exposure)
    dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                     dat$se.outcome,
                                     dat$samplesize.outcome)
    res_steiger <- mr_steiger(
      p_exp = dat$pval.exposure,
      p_out = dat$pval.outcome,
      n_exp = dat$samplesize.exposure,
      n_out = dat$samplesize.outcome,
      r_exp = dat$r.exposure,
      r_out = dat$r.outcome
    )
    res_steiger <- directionality_test(dat)
    
    
    
    
    
    res_presso <- TwoSampleMR::run_mr_presso(dat, NbDistribution = 100)
    
    sink(paste0("final result/",out_come,"_PRESSO.txt"),append=FALSE,split = FALSE) 
    print(res_presso)
    sink()
    print(res_presso)
    
    
    
    
    p1 <- mr_scatter_plot(res, dat)
    p1[[1]]
    pdf(paste0("final result/",out_come,"_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    
    res_single <- mr_singlesnp(dat)
    p2 <- mr_forest_plot(res_single)
    pdf(paste0("final result/",out_come,"_forest.pdf"))
    print(p2[[1]])
    dev.off()
    
    p3 <- mr_funnel_plot(res_single)
    pdf(paste0("final result/",out_come,"_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    
    res_loo <- mr_leaveoneout(dat)
    pdf(paste0("final result/",out_come,"_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
    
    
    library(magrittr)
    res3 <- res[1:3,]
    
    
    
    
    library(magrittr)
    
    res4 <- tidyr::pivot_wider(
      res3,names_from ="method",names_vary = "slowest",
      values_from = c("b","se","pval") )
    
    res_hete2 <- tidyr::pivot_wider(
      res_hete,names_from ="method",names_vary = "slowest",
      values_from = c("Q","Q_df","Q_pval") ) %>% 
      dplyr::select( -id.exposure,-id.outcome,-outcome,-exposure)
    
    res_plei2 <- dplyr::select(res_plei,
                               egger_intercept,se,pval)
    
   
    res_steiger2 <- dplyr::select(res_steiger,
                                  correct_causal_direction,steiger_pval)
    
    
   
    res_ALL <- cbind(res4, res_hete2, res_plei2,res_steiger2)
    
    
    write.csv(res_ALL,file = paste0("final result/",out_come,".csv"), row.names = FALSE)
  }}




folder_path <- "E:\\41种炎症因子最新\\血栓\\反向\\PE"


fs <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)


df <- map_dfr(fs, read.csv)




write.csv(df, "E:\\41种炎症因子最新\\血栓\\反向\\PE\\res.csv", row.names = FALSE)

