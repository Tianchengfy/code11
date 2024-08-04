library(readxl)
library(writexl)

# 设置工作目录
setwd("E:\\41种炎症因子最新\\血栓\\正向\\DVT\\DVT整理\\多效性合并")

# 获取文件夹中的所有excel文件
file_list <- list.files(pattern="*.xlsx")

# 读入所有excel文件，并为每个文件添加一列文件名
data <- lapply(file_list, function(file) {
  df <- read_excel(file)
  df$file_name <- file # 添加文件名作为新的列
  return(df)
})

# 合并数据
combined_data <- do.call(rbind, data)
# 将合并后的数据保存为excel文件
write_xlsx(combined_data, "多效性合并.xlsx")
