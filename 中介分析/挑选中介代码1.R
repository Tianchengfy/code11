# 加载所需的包
library(readr)
library(dplyr)
setwd("E:\\代谢物最新\\中介分析\\VTE")

# 定义要查找的新ID列表
ids_to_find <- c("GCST90200984", "GCST90200976", "GCST90201018", "GCST90199683", "GCST90200416",
                 "GCST90200386", "GCST90200653", "GCST90199622", "GCST90200030", "GCST90200429",
                 "GCST90199969", "GCST90199786", "GCST90200771", "GCST90200155", "GCST90200082",
                 "GCST90200906", "GCST90200325", "GCST90200597", "GCST90200468", "GCST90200915",
                 "GCST90199694", "GCST90199844", "GCST90199944", "GCST90199811", "GCST90200037",
                 "GCST90200307", "GCST90200505", "GCST90200656", "GCST90199760", "GCST90199933",
                 "GCST90200777", "GCST90200598", "GCST90200667", "GCST90200026", "GCST90200475",
                 "GCST90200067", "GCST90200738", "GCST90200951", "GCST90200602", "GCST90200489",
                 "GCST90199689")

# 使用read_csv函数读取CSV文件
df <- read_csv("FDR校正后结果significant.csv")

# 筛选出FileName列中值在ids_to_find列表中的行
filtered_df <- filter(df, FileName %in% ids_to_find)

# 将筛选后的数据写入新的CSV文件
write_csv(filtered_df, "中介.csv")

# 打印消息表示操作完成
print("筛选完成，结果已写入中介.csv")
