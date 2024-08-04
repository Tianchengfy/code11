library(readxl)
library(writexl)


setwd("E:\\代谢物最新\\正向\\代谢物对VTE\\VTE整理\\多效性合并")


file_list <- list.files(pattern="*.xlsx")


data <- lapply(file_list, function(file) {
  df <- read_excel(file)
  df$file_name <- file 
  return(df)
})


combined_data <- do.call(rbind, data)

write_xlsx(combined_data, "多效性合并.xlsx")
