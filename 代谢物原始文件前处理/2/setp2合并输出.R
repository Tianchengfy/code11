
library(dplyr)
setwd("E:/代谢物最新/去混杂/clumped/out")

txt_files <- list.files(pattern = "\\.txt$", full.names = TRUE)
combined_df <- data.frame()
for (file in txt_files) {
  df <- read.csv(file, sep = "\t") 
  combined_df <- rbind(combined_df, df)
}


write.csv(combined_df, file = "bb.csv")
