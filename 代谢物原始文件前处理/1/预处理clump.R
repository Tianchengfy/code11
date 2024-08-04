

setwd(dir="E:/代谢物最新/前处理/1e05全") 


library(ieugwasr)
library(data.table)


gwas_files <- list.files(pattern = "\\.txt$")


for(file_name in gwas_files) {
  ldl_gwas <- fread(file_name, sep = "\t", header = TRUE)
  
  
  ldl_gwas$pval <- as.numeric(ldl_gwas$pval)
  ldl_gwas$eaf <- as.numeric(ldl_gwas$eaf)
  ldl_gwas$se <- as.numeric(ldl_gwas$se)
  ldl_gwas <- subset(ldl_gwas, pval < 1e-05)
  
  ldl_iv <- ldl_gwas[, c("SNP", "pval")]
  colnames(ldl_iv) <- c("rsid", "pval")
  
  
  clumped_data <- ld_clump(
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR",
    dplyr::tibble(rsid=ldl_iv$rsid, pval=ldl_iv$pval),
    plink_bin = "C:/Windows/本地plink/plink_win64_20230116/plink.exe",
    bfile = "C:/Windows/本地plink/EUR/EUR"
  )
  
  filtered_gwas <- ldl_gwas[ldl_gwas$SNP %in% clumped_data$rsid, ]
  
  
  formatted_data <- format_data(
    filtered_gwas,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    samplesize_col = "samplesize",
    phenotype_col = "Phenotype",
    log_pval = FALSE
  )
  
  
  
 
  print(head(formatted_data))
  
  
  fwrite(formatted_data, file = paste0("clumped_10000kb_r20.001/processed_", file_name), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
