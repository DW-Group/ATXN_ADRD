library(dplyr)
library(readxl)
library(writexl)
library(metafor)

rm(list = ls())

load("data/NULISA_core.RData")

##-------------------------------
## Results from NHS substudies
##-------------------------------

diet_marker_list <- unique(c(
  "all_core3_ew_signed",
  core_list[["all"]][["core3"]],
  core_list[["BD"]][["core3"]]
))

diet_var_list <- unique(c(
  dietscore_list,
  amedcomp_list,
  mindcomp_list
))

pros_res <- read_excel("results/nhs_NULISA_diet.xlsx") %>% 
  filter(
    scale == "scaled",
    model == "univariate",
    targetName %in% diet_marker_list,
    diet %in% diet_var_list
  ) %>% 
  dplyr::select(-scale) %>% 
  mutate(
    n = as.numeric(n),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  )

caco_res <- read_excel("results/nhs_caco_NULISA_diet.xlsx") %>% 
  filter(
    analysis_set == "control_only",
    scale == "scaled",
    model == "univariate",
    targetName %in% diet_marker_list,
    diet %in% diet_var_list
  ) %>% 
  dplyr::select(-analysis_set, -scale) %>% 
  mutate(
    n = as.numeric(n),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  )

##-------------------------------
## Fixed-effect meta-analysis
##-------------------------------

id_cols <- c("targetName", "diet", "model")

anno_cols <- c(
  "targetLongName",
  "uniprotID",
  "proteinName",
  "geneName",
  "targetName_orig"
)

meta_dat_long <- bind_rows(
  pros_res %>% mutate(study = "pros"),
  caco_res %>% mutate(study = "caco")
)

meta_res <- meta_dat_long %>%
  group_by(across(all_of(id_cols))) %>%
  group_modify(~{
    
    dat <- .x
    
    if (n_distinct(dat$study) < 2) {
      return(tibble())
    }
    
    fit <- rma.uni(
      yi = dat$beta,
      sei = dat$se,
      method = "FE"
    )
    
    tibble(
      n = sum(dat$n, na.rm = TRUE),
      LOD_flag = any(dat$LOD_flag, na.rm = TRUE),
      
      targetLongName = dplyr::coalesce(dat$targetLongName[1], dat$targetLongName[2]),
      uniprotID      = dplyr::coalesce(dat$uniprotID[1], dat$uniprotID[2]),
      proteinName    = dplyr::coalesce(dat$proteinName[1], dat$proteinName[2]),
      geneName       = dplyr::coalesce(dat$geneName[1], dat$geneName[2]),
      targetName_orig = dplyr::coalesce(dat$targetName_orig[1], dat$targetName_orig[2]),
      
      n_pros    = dat$n[dat$study == "pros"][1],
      beta_pros = dat$beta[dat$study == "pros"][1],
      se_pros   = dat$se[dat$study == "pros"][1],
      p_pros    = dat$p[dat$study == "pros"][1],
      
      n_caco    = dat$n[dat$study == "caco"][1],
      beta_caco = dat$beta[dat$study == "caco"][1],
      se_caco   = dat$se[dat$study == "caco"][1],
      p_caco    = dat$p[dat$study == "caco"][1],
      
      beta_meta  = as.numeric(fit$b),
      se_meta    = fit$se,
      z_meta     = fit$zval,
      p_meta     = fit$pval,
      lower_meta = fit$ci.lb,
      upper_meta = fit$ci.ub,
      QE         = fit$QE,
      QEp        = fit$QEp,
      I2         = fit$I2
    )
  }) %>%
  ungroup() %>%
  dplyr::select(
    all_of(id_cols),
    all_of(anno_cols),
    n,
    beta_meta,
    se_meta,
    z_meta,
    p_meta,
    lower_meta,
    upper_meta,
    QE,
    QEp,
    I2,
    LOD_flag,
    n_pros,
    beta_pros,
    se_pros,
    p_pros,
    n_caco,
    beta_caco,
    se_caco,
    p_caco
  )

write_xlsx(
  meta_res,
  "results/nhs_meta_NULISA_diet.xlsx"
)