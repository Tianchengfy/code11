
library(TwoSampleMR) 
library(data.table)
library(R.utils)
setwd("E:\\VTE投稿返修\\测试")

local_outcome_data <- fread("finngen_R11_I9_VTE.txt", header = TRUE) 


exposure_id <- "ebi-a-GCST004432"
flag <- TRUE


while (flag) {
  tryCatch({
    exposure_dat <- extract_instruments(outcomes = exposure_id, p1 = 5e-06, clump = TRUE, access_token = NULL)
    flag <- FALSE
  }, error = function(e) {
    message("Error: ", e)
  })
}


print(head(exposure_dat))


if (!"rsids" %in% colnames(local_outcome_data)) {
  stop("The local outcome data does not contain 'rsids' column. Please check the column names.")
}


outcome_dat <- merge(exposure_dat, local_outcome_data, by.x = "SNP", by.y = "rsids") 
write.csv(outcome_dat, file = "merged_local_data.csv")


outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "merged_local_data.csv",
  sep = ",",
  pval_col = "pval",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
) 


print(head(outcome_dat))


dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
print(head(dat))


if (length(which(dat$mr_keep == TRUE)) > 0) {
  res <- generate_odds_ratios(mr_res = mr(dat, method_list = c("mr_ivw",
                                                               "mr_egger_regression",
                                                               "mr_weighted_median",
                                                               "mr_weighted_mode",
                                                               "mr_simple_mode",
                                                               "mr_wald_ratio")))
  write.csv(res, file = "res_VTE11.csv", row.names = FALSE)
  write.csv(dat, file = "dat_VTE11.csv", row.names = FALSE)
}
