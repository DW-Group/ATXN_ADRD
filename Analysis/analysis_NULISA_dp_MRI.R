library(dplyr)
library(tidyverse)
library(data.table)
library(writexl)

rm(list = ls())

##----------------------------------
## Load data
##----------------------------------

load("data/data_NULISA_dp.RData")

dp_score_dat <- read.csv("data/data_NULISA_dp_score.csv")

npq_all <- npq_all %>% 
  left_join(dp_score_dat, by = "tsno")

protall_list <- unique(c(prot_list, prot_list_add))
protall_list <- protall_list[!grepl("raw", protall_list)]

scoreall_list <- setdiff(colnames(dp_score_dat), "tsno")

mri_raw <- fread("DIRECT-PLUS_Metadata/brain_mri.csv") %>% 
  mutate(
    HippocampusTotalVolume_change =
      100 * (HippocampusTotalVolume_time18 - HippocampusTotalVolume_time0) /
      HippocampusTotalVolume_time0,
    InferiorLateralVentricleTotalVolume_change =
      100 * (InferiorLateralVentricleTotalVolume_time18 - InferiorLateralVentricleTotalVolume_time0) /
      InferiorLateralVentricleTotalVolume_time0,
    HOC_score_change_right =
      100 * (HOC_score_time_18_right - HOC_score_time_0_right) /
      HOC_score_time_0_right,
    HOC_score_change_left =
      100 * (HOC_score_time_18_left - HOC_score_time_0_left) /
      HOC_score_time_0_left,
    HOC_mean_change =
      100 * (HOC_mean_time18 - HOC_mean_time0) /
      HOC_mean_time0,
    group = case_when(
      group == 1 ~ "HDG",
      group == 2 ~ "MED",
      group == 3 ~ "GreenMED"
    )
  )

mri_var_dict <- data.frame(
  time0 = c(
    "HippocampusTotalVolume_time0",
    "InferiorLateralVentricleTotalVolume_time0",
    "HOC_score_time_0_right",
    "HOC_score_time_0_left",
    "HOC_mean_time0"
  ),
  time18 = c(
    "HippocampusTotalVolume_time18",
    "InferiorLateralVentricleTotalVolume_time18",
    "HOC_score_time_18_right",
    "HOC_score_time_18_left",
    "HOC_mean_time18"
  ),
  change = c(
    "HippocampusTotalVolume_change",
    "InferiorLateralVentricleTotalVolume_change",
    "HOC_score_change_right",
    "HOC_score_change_left",
    "HOC_mean_change"
  )
)

outlier_counts <- list()

for (v in c(mri_var_dict[["time0"]], mri_var_dict[["time18"]])) {
  
  if (grepl("HOC_mean", v)) next
  
  z_v <- (mri_raw[[v]] - mean(mri_raw[[v]], na.rm = TRUE)) /
    sd(mri_raw[[v]], na.rm = TRUE)
  
  idx_out <- abs(z_v) > 2.5 & !is.na(z_v)
  outlier_counts[[v]] <- sum(idx_out)
  mri_raw[[v]][idx_out] <- NA
}

mri_outlier_summary <- data.frame(
  var = names(outlier_counts),
  n_outlier = unlist(outlier_counts)
)

mri_raw <- mri_raw %>% 
  mutate(
    HippocampusTotalVolume_change =
      if_else(
        is.na(HippocampusTotalVolume_time0) | is.na(HippocampusTotalVolume_time18),
        NA_real_,
        HippocampusTotalVolume_change
      ),
    InferiorLateralVentricleTotalVolume_change =
      if_else(
        is.na(InferiorLateralVentricleTotalVolume_time0) |
          is.na(InferiorLateralVentricleTotalVolume_time18),
        NA_real_,
        InferiorLateralVentricleTotalVolume_change
      ),
    HOC_score_change_right =
      if_else(
        is.na(HOC_score_time_0_right) | is.na(HOC_score_time_18_right),
        NA_real_,
        HOC_score_change_right
      ),
    HOC_score_change_left =
      if_else(
        is.na(HOC_score_time_0_left) | is.na(HOC_score_time_18_left),
        NA_real_,
        HOC_score_change_left
      ),
    HOC_mean_time0 =
      if_else(
        is.na(HOC_score_time_0_left) | is.na(HOC_score_time_0_right),
        NA_real_,
        HOC_mean_time0
      ),
    HOC_mean_time18 =
      if_else(
        is.na(HOC_score_time_18_left) | is.na(HOC_score_time_18_right),
        NA_real_,
        HOC_mean_time18
      ),
    HOC_mean_change =
      if_else(
        is.na(HOC_mean_time0) | is.na(HOC_mean_time18),
        NA_real_,
        HOC_mean_change
      )
  )

source("scripts/data_genetic_dp.R")

gene_list <- c(
  "APOE4_carrier",
  "PGS002280_scaled",
  "PGS000334_scaled",
  ad_dict$rsid_ref
)

genecore_list <- c(
  "APOE4_carrier",
  "PGS002280_scaled",
  "PGS000334_scaled"
)

source("scripts/data_diet_dp.R")

diet_list <- setdiff(colnames(direct_diet), c("sno", "calor"))
dietscore_list <- c("amed", "MIND", "MIND14")
amedcomp_list <- grep("^amedcomp_", colnames(direct_diet), value = TRUE)
mindcomp_list <- grep("^mindcomp_", colnames(direct_diet), value = TRUE)

covar_base <- c("age", "gender", "HOMAIR_0")
covar_gene <- c("PC1", "PC2", "PC3", "PC4")
covar_diet <- c("calor")

save.image(file = "data/processed_data_dp_MRI.RData")

##----------------------------------
## Baseline NULISA and MRI
##----------------------------------

baseline_data <- npq_all %>% 
  filter(time == "A") %>% 
  left_join(mri_raw, by = c("sno", "group")) %>% 
  left_join(allgeno, by = "sno", suffix = c("", "_genetic")) %>% 
  left_join(direct_diet, by = "sno") %>% 
  filter(
    !is.na(group),
    !is.na(HOMAIR_0)
  ) %>% 
  mutate(
    group = factor(group, levels = c("HDG", "MED", "GreenMED"))
  ) %>%
  mutate(
    dplyr::across(
      all_of(c(protall_list, scoreall_list)),
      ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  )

baseline_data_list <- list(
  "age_ge_50" = list(
    "all" = baseline_data %>% filter(age_baseline >= 50)
  )
)

##----------------------------------
## MRI change analyses
##----------------------------------

mri_analysis_list <- c(
  "HOC_mean_change",
  "InferiorLateralVentricleTotalVolume_change"
)

prot_analysis_list <- unique(c(
  "all_core3_ew_signed",
  core_list[["all"]][["core3"]],
  core_list[["BD"]][["core3"]]
))

base_res_list <- list()
i <- 1

for (prot in prot_analysis_list) {
  
  message(prot)
  
  for (mri in mri_analysis_list) {
    
    for (g in levels(baseline_data$group)) {
      
      tmp_data <- baseline_data_list[["age_ge_50"]][["all"]]
      tmp_data <- tmp_data %>% filter(group == g)
      
      tmp_covar <- c("age_baseline", "gender", "HOMAIR_0")
      
      prot_z <- paste0(prot, "_z")
      
      tmp_mod_scaled <- tryCatch(
        glm(
          as.formula(
            paste(
              mri,
              " ~ ",
              paste(c(prot_z, tmp_covar), collapse = " + ")
            )
          ),
          data = tmp_data
        ),
        error = function(e) e,
        warning = function(w) w
      )
      
      base_res_list[[i]] <- c(
        prot,
        "scaled",
        mri,
        g,
        "age_ge_50",
        "all",
        if (inherits(tmp_mod_scaled, "glm"))
          c(
            nobs(tmp_mod_scaled),
            summary(tmp_mod_scaled)$coefficients[prot_z, c(1, 2, 4)]
          )
        else c(tmp_mod_scaled$message, rep(NA, 3))
      )
      
      i <- i + 1
    }
  }
}

base_res_dat <- as.data.frame(do.call("rbind", base_res_list))

colnames(base_res_dat) <- c(
  "targetName",
  "scale",
  "MRI",
  "group",
  "age_group",
  "gender_group",
  "n",
  "beta",
  "se",
  "p"
)

base_res_dat <- base_res_dat %>% 
  left_join(marker_dict, by = "targetName") %>% 
  mutate(
    n = as.numeric(n),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  )

write_xlsx(base_res_dat, "results/dp_MRI_NULISA_baseline.xlsx")