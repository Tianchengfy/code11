library(readr)
library(writexl) 


setwd("E:\\代谢物最新\\正向\\代谢物对PE\\PE整理\\异质性合并")


file_list <- list.files(pattern="*.csv")


data <- lapply(file_list, function(file) {
  df <- read_csv(file)
  df$file_name <- file 
  return(df)
})


combined_data <- do.call(rbind, data)


write_xlsx(combined_data, "异质性合并.xlsx")
