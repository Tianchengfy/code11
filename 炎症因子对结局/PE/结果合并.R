library(readxl)
library(writexl)
library(dplyr)

setwd("E:\\41种炎症因子最新\\血栓\\正向\\PE\\PE整理\\初筛合并")

file_list <- list.files(pattern = "\\.xlsx$")

data <- lapply(file_list, function(file) {
  df <- read_excel(file)
  df$file_name <- file
  return(df)
})

combined_data <- bind_rows(data)

write_xlsx(combined_data, "结果合并.xlsx")
