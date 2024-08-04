
library(TwoSampleMR) 
library(data.table)
library(R.utils)
setwd("E:\\VTE投稿返修\\测试")


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


outcome_id <- "ebi-a-GCST90038615"
flag <- TRUE


while (flag) {
  tryCatch({
    outcome_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = outcome_id,
      access_token = NULL
    )
    if (length(outcome_dat) > 0) {
      flag <- FALSE
    } else {
      flag <- TRUE 
    }
  }, error = function(e) {
    message("Error: ", e)
  })
}


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
  write.csv(res, file = "res_local.csv", row.names = FALSE)
  write.csv(dat, file = "dat_local.csv", row.names = FALSE)
}

write.csv(res, file = "res_DVT.csv", row.names = FALSE)
write.csv(dat, file = "dat_DVT.csv", row.names = FALSE)
