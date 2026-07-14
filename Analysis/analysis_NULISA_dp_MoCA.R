library(dplyr)
library(data.table)
library(stringr)
library(lubridate)
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

##----------------------------------
## MoCA total z-score
##----------------------------------

fit_cog <- fread("DIRECT-PLUS_brain_microbiome/moca_test.csv") %>% 
  mutate(
    sno = as.numeric(gsub("DIRECT_PLUS_", "", Trial_SNO)),
    MOCA_Date_clean = str_extract(MOCA_Date, "^\\d{2}\\.\\d{2}\\.\\d{2,4}"),
    MOCA_Date_parsed = case_when(
      str_detect(MOCA_Date_clean, "\\.\\d{4}$") ~ dmy(MOCA_Date_clean),
      str_detect(MOCA_Date_clean, "\\.\\d{2}$") ~ dmy(MOCA_Date_clean),
      TRUE ~ as.Date(NA)
    ),
    MOCA_year  = year(MOCA_Date_parsed),
    MOCA_month = month(MOCA_Date_parsed),
    MOCA_months_since_1900 = (MOCA_year - 1900) * 12 + MOCA_month,
    MOCA_months_since_2017 = MOCA_months_since_1900 - (117 * 12 + 5),
    
    moca_visexec = space1_1p + space2_1p + space3_1p + space4_1p + space5_1p,
    moca_naming  = name1_1p  + name2_1p  + name3_1p,
    moca_attention = attention1_1p + attention2_1p + attention3_1p + attention4_3p,
    moca_language  = Lang1_1p + Lang2_1p + Lang3_1p,
    moca_abstraction = Abstract1_1p + Abstract2_1p,
    moca_memory_delayed = Memory1_5p,
    moca_orientation = Orientation1_1p + Orientation2_1p + Orientation3_1p +
      Orientation4_1p + Orientation5_1p + Orientation6_1p,
    
    moca_visexec_z = as.numeric(scale(moca_visexec)),
    moca_naming_z = as.numeric(scale(moca_naming)),
    moca_attention_z = as.numeric(scale(moca_attention)),
    moca_language_z = as.numeric(scale(moca_language)),
    moca_abstraction_z = as.numeric(scale(moca_abstraction)),
    moca_memory_delayed_z = as.numeric(scale(moca_memory_delayed)),
    moca_orientation_z = as.numeric(scale(moca_orientation)),
    
    Total_z = rowMeans(
      cbind(
        moca_visexec_z,
        moca_naming_z,
        moca_attention_z,
        moca_language_z,
        moca_abstraction_z,
        moca_memory_delayed_z,
        moca_orientation_z
      ),
      na.rm = TRUE
    ),
    Total_z = as.numeric(scale(Total_z)),
    edu = `if ≤12_1p`
  ) %>%
  select(
    sno,
    Total_z,
    edu,
    MOCA_months_since_2017
  )

age_group_info <- npq_all %>% 
  select(sno, age_baseline, group, gender) %>% 
  distinct()

##----------------------------------
## FIT NULISA and MoCA
##----------------------------------

cross_data <- npq_all %>% 
  filter(time == "D") %>% 
  select(-group) %>% 
  left_join(age_group_info %>% select(sno, group), by = "sno") %>% 
  left_join(fit_cog, by = "sno") %>% 
  filter(!is.na(group), !is.na(Total_z)) %>% 
  mutate(
    age_FIT = age_baseline + MOCA_months_since_2017 / 12,
    MOCA_months_since_month18 = MOCA_months_since_2017 - 18,
    group = factor(group, levels = c("HDG", "MED", "GreenMED"))
  ) %>%
  mutate(
    across(
      all_of(c(protall_list, scoreall_list)),
      ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  filter(age_FIT >= 50)

##----------------------------------
## DIRECT-PLUS MoCA total z-score analysis
##----------------------------------

main_res_list <- list()
i <- 1

for (prot in c(protall_list, scoreall_list)) {
  
  message(prot)
  
  prot_z <- paste0(prot, "_z")
  
  tmp_data <- cross_data %>%
    filter(
      !is.na(.data[[prot_z]]),
      !is.na(Total_z),
      !is.na(gender),
      !is.na(group),
      !is.na(age_FIT),
      !is.na(MOCA_months_since_month18),
      !is.na(edu)
    )
  
  tmp_mod <- tryCatch(
    lm(
      as.formula(
        paste(
          "Total_z ~",
          paste(
            c(prot_z, "gender", "group", "age_FIT", "MOCA_months_since_month18", "edu"),
            collapse = " + "
          )
        )
      ),
      data = tmp_data
    ),
    error = function(e) e
  )
  
  main_res_list[[i]] <- c(
    targetName = prot,
    scale = "scaled",
    cog_type = "zscore",
    cog = "Total_z",
    group = "all",
    age_group = "age_ge_50",
    age_covar = "linear",
    gender_group = "all",
    n = if (inherits(tmp_mod, "lm")) nobs(tmp_mod) else NA,
    beta = if (inherits(tmp_mod, "lm")) summary(tmp_mod)$coefficients[prot_z, "Estimate"] else NA,
    se = if (inherits(tmp_mod, "lm")) summary(tmp_mod)$coefficients[prot_z, "Std. Error"] else NA,
    p = if (inherits(tmp_mod, "lm")) summary(tmp_mod)$coefficients[prot_z, "Pr(>|t|)"] else NA,
    status = if (inherits(tmp_mod, "lm")) "OK" else paste0("MODEL ERROR: ", tmp_mod$message)
  )
  
  i <- i + 1
}

main_res_dat <- as.data.frame(do.call(rbind, main_res_list), stringsAsFactors = FALSE)

main_res_dat <- main_res_dat %>%
  mutate(
    n = as.numeric(n),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(
  main_res_dat,
  "results/dp_MoCA_NULISA_cross_Total_z.xlsx"
)