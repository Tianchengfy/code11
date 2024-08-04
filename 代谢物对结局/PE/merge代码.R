library(readr)
library(writexl) 

setwd("E:\\41种炎症因子最新\\血栓\\正向\\PE\\结果merge")


file1 <- read_csv("1.csv")
file2 <- read_csv("2.csv")
file3 <- read_csv("3.csv")


merged_data_12 <- merge(file1, file2, by = "file_name")


merged_data <- merge(merged_data_12, file3, by = "file_name")


print(merged_data)


write_csv(merged_data, "所有合并.csv")

