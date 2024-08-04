library(readxl) # 用于读取Excel文件
library(writexl) # 用于写入Excel文件
library(dplyr) # 提供了一些便捷的数据处理功能

# 设置工作目录到存放Excel文件的文件夹
setwd("E:\\41种炎症因子最新\\血栓\\正向\\DVT\\DVT整理\\初筛合并")

# 获取文件夹中的所有Excel文件
file_list <- list.files(pattern="\\.xlsx$") # 注意：这里的正则表达式确保只匹配.xlsx文件

# 读入所有Excel文件，并为每个文件添加一列文件名
data <- lapply(file_list, function(file) {
  df <- read_excel(file) # 使用read_excel函数读取Excel文件
  df$file_name <- file # 添加文件名作为新的列
  return(df)
})

# 使用dplyr的bind_rows而不是base R的do.call(rbind, ...)来合并数据框
combined_data <- bind_rows(data)

# 将合并后的数据保存为一个新的Excel文件
write_xlsx(combined_data, "结果合并.xlsx")
