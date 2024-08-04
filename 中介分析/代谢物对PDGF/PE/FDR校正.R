# The Benjamini-Hochberg method
library(writexl)
setwd("F:/中介分析/代谢物对PDGF/PE")
# 1.读取结果文件，注意，要带X列
res<-read.csv("FDR校正.csv")

# 2.选择IVW方法的结果，其他方法的结果也是一样的操作
res_IVW<-res[res$method=="Inverse variance weighted",] 

# 3.创建一个只有X和pval列的数据框
res_p<-res_IVW[,c("X","pval")]

# 4.按pval从小到大排序
res_p<-res_p[order(res_p$pval),]

# 5.创建一个adjustP值列
res_p$adjustP<-NA

# 6.开始计算FDR
for (n in 1:nrow(res_p)) {
  # 最大的P值就是最大的FDR值，也就是最后一行的P值设置为FDR
  if(n==nrow(res_p)){
    res_p[n,]$adjustP<-res_p[n,]$pval
    
  }else{  # 如果不是最后一行，那么要下列方法计算
    # 计算公式 adjustP = pval * ((p值总个数)/(当前pval的所在行数))
    res_p[n,]$adjustP<-((res_p[n,]$pval)*(nrow(res_p)/n))
  }
  # 如果当前adjustP值大于1，那就设置为1，因为P值范围是0-1，就算纠正之后最大也只能是1
  if(res_p[n,]$adjustP>=1)res_p[n,]$adjustP<-1
}

# 7.将结果返回到res里面
res<-merge(x=res,y=res_p[,c("X","adjustP")],by = "X",all = T)


# 同理，WM方法，mr_egger也是这样子操作



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

# Simple mode FDR校正
res_Simple <- res[res$method == "Simple mode",]
res_p_Simple <- res_Simple[,c("X","pval")]
res_p_Simple <- res_p_Simple[order(res_p_Simple$pval),]
res_p_Simple$adjustP <- p.adjust(res_p_Simple$pval, method = "fdr")
res <- merge(x=res, y=res_p_Simple[,c("X","adjustP")], by = "X", all = T)

# Weighted mode FDR校正
res_Weighted <- res[res$method == "Weighted mode",]
res_p_Weighted <- res_Weighted[,c("X","pval")]
res_p_Weighted <- res_p_Weighted[order(res_p_Weighted$pval),]
res_p_Weighted$adjustP <- p.adjust(res_p_Weighted$pval, method = "fdr")
res <- merge(x=res, y=res_p_Weighted[,c("X","adjustP")], by = "X", all = T)


write_xlsx(res, "FDR校正后结果新.xlsx")

