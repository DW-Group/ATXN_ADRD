library(dplyr)
library(data.table)
library(writexl)
library(TwoSampleMR)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/mr.RData")

##-----------------------------------------
## Harmonize data
##-----------------------------------------

harmonized_list <- list()
for (pheno in names(outcome_list)){
  harmonized_list[[pheno]] <- harmonise_data(exposure_dat=instru_pqtl_list[[pheno]], outcome_dat=outcome_list[[pheno]])
}

##-----------------------------------------
## MR analysis
##-----------------------------------------

method_list <- c("mr_wald_ratio","mr_ivw")

mr_res_list <- list()
for (pheno in names(harmonized_list)){
  mr_res_list[[pheno]] <- mr(harmonized_list[[pheno]], method_list=method_list)
}

mr_res_dat <- as.data.frame(do.call("rbind",mr_res_list))

write_xlsx(mr_res_dat, "results/mr/mr_res.xlsx")
