
folder_path <- "E:/代谢物最新/前处理/clumped_10000kb_r20.001"


new_folder_path <- "E:/代谢物最新/前处理/clumped"


file_names <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)


processed_file_names <- c()


for (file_path in file_names) {
 
  data <- read.table(file_path, header = TRUE, sep = "\t", quote = "")
  
 
  file_name <- basename(file_path)
  
 
  extracted_name <- sub("processed_(GCST[0-9]+)_.*", "\\1", file_name)
  
  
  data$processed_name <- extracted_name
  
  
  processed_file_names <- c(processed_file_names, extracted_name)
  
  
  new_file_path <- paste0(new_folder_path, "/", extracted_name, ".txt")
  
  
  write.table(data, new_file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


processed_file_names_with_extension <- paste0(processed_file_names, ".txt")


write(processed_file_names_with_extension, file = "1.txt")
