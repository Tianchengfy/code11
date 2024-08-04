
library(writexl)
setwd("E:\\41种炎症因子最新\\血栓\\正向\\VTE\\FDR")

res<-read.csv("结果合并.csv")


res_IVW<-res[res$method=="Inverse variance weighted",] 


res_p<-res_IVW[,c("X","pval")]


res_p<-res_p[order(res_p$pval),]


res_p$adjustP<-NA


for (n in 1:nrow(res_p)) {
  
  if(n==nrow(res_p)){
    res_p[n,]$adjustP<-res_p[n,]$pval
    
  }else{  
    res_p[n,]$adjustP<-((res_p[n,]$pval)*(nrow(res_p)/n))
  }
 
  if(res_p[n,]$adjustP>=1)res_p[n,]$adjustP<-1
}


res<-merge(x=res,y=res_p[,c("X","adjustP")],by = "X",all = T)






res_WM <- res[res$method == "Weighted median",]
res_p <- res_WM[,c("X","pval")]
res_p <- res_p[order(res_p$pval),]
res_p$adjustP <- p.adjust(res_p$pval, method = "fdr")
res <- merge(x=res, y=res_p[,c("X","adjustP")], by = "X", all = T)


res_mr_egger <- res[res$method == "MR Egger",]
res_p <- res_mr_egger[,c("X","pval")]
res_p <- res_p[order(res_p$pval),]
res_p$adjustP <- p.adjust(res_p$pval, method = "fdr")
res <- merge(x=res, y=res_p[,c("X","adjustP")], by = "X", all = T)


res_Simple <- res[res$method == "Simple mode",]
res_p_Simple <- res_Simple[,c("X","pval")]
res_p_Simple <- res_p_Simple[order(res_p_Simple$pval),]
res_p_Simple$adjustP <- p.adjust(res_p_Simple$pval, method = "fdr")
res <- merge(x=res, y=res_p_Simple[,c("X","adjustP")], by = "X", all = T)


res_Weighted <- res[res$method == "Weighted mode",]
res_p_Weighted <- res_Weighted[,c("X","pval")]
res_p_Weighted <- res_p_Weighted[order(res_p_Weighted$pval),]
res_p_Weighted$adjustP <- p.adjust(res_p_Weighted$pval, method = "fdr")
res <- merge(x=res, y=res_p_Weighted[,c("X","adjustP")], by = "X", all = T)



write_xlsx(res, "FDR校正后结果.xlsx")

