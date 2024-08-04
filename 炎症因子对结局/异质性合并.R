library(readr) 
library(writexl) 


setwd("E:\\41种炎症因子最新\\血栓\\反向\\PE\\final result\\异质性合并")


file_list <- list.files(pattern="*.csv")


data <- lapply(file_list, function(file) {
  df <- read_csv(file)
  df$file_name <- file 
  return(df)
})


combined_data <- do.call(rbind, data)


write_xlsx(combined_data, "异质性合并.xlsx")

