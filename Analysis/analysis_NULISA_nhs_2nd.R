library(dplyr)
library(tidyverse)
library(data.table)
library(survival)
library(readxl)
library(writexl)

rm(list = ls())

##-------------------------------
## Load data
##-------------------------------

load("data/processed_data_nhs_2nd_calibrated.RData")

pc <- read_excel("data/nhs_global_change.xlsx")

model_data <- model_data %>% 
  left_join(pc, by = "id")

##-------------------------------
## Change
##-------------------------------

### Paired t-test and Wilcoxon

res_paired <- map_dfr(c(protall_list, scoreall_list), function(p) {
  
  map_dfr(c("all", "0", "1"), function(apoe_grp) {
    
    dat_sub <- if (apoe_grp == "all") {
      model_data
    } else {
      model_data %>% filter(APOE4_carrier == as.numeric(apoe_grp))
    }
    
    x1 <- dat_sub[[paste0(p, "_1")]]
    x2 <- dat_sub[[paste0(p, "_2")]]
    chg <- dat_sub[[paste0(p, "_chg")]]
    chg_yr <- dat_sub[[paste0(p, "_chg_yr")]]
    
    tt <- tryCatch(
      t.test(x2, x1, paired = TRUE),
      error = function(e) NULL
    )
    
    wt <- tryCatch(
      wilcox.test(x2, x1, paired = TRUE, exact = FALSE),
      error = function(e) NULL
    )
    
    tibble(
      marker = p,
      APOE4_group = apoe_grp,
      n = sum(!is.na(x1) & !is.na(x2)),
      mean_1 = mean(x1, na.rm = TRUE),
      mean_2 = mean(x2, na.rm = TRUE),
      mean_change = mean(x2 - x1, na.rm = TRUE),
      median_change = median(chg, na.rm = TRUE),
      median_change_yr = median(chg_yr, na.rm = TRUE),
      t = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
      p_t = if (!is.null(tt)) tt$p.value else NA_real_,
      w = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
      p_w = if (!is.null(wt)) wt$p.value else NA_real_
    )
  })
}) %>%
  group_by(APOE4_group) %>%
  mutate(
    p_t_fdr = p.adjust(p_t, method = "fdr"),
    p_w_fdr = p.adjust(p_w, method = "fdr")
  ) %>%
  ungroup()

write_xlsx(res_paired, "results/nhs_2nd_NULISA_paired_tests.xlsx")

### Global median summaries

chg_summary <- map_dfr(c(protall_list, scoreall_list), function(p) {
  
  tibble(
    marker = p,
    median_chg = median(model_data[[paste0(p, "_chg")]], na.rm = TRUE),
    median_chg_yr = median(model_data[[paste0(p, "_chg_yr")]], na.rm = TRUE),
    IQR_chg = IQR(model_data[[paste0(p, "_chg")]], na.rm = TRUE),
    IQR_chg_yr = IQR(model_data[[paste0(p, "_chg_yr")]], na.rm = TRUE)
  )
})

write_xlsx(chg_summary, "results/nhs_2nd_NULISA_change.xlsx")

### Direction

direction_summary <- map_dfr(c(protall_list, scoreall_list), function(p) {
  
  chg <- model_data[[paste0(p, "_chg")]]
  
  tibble(
    marker = p,
    prop_increase = mean(chg > 0, na.rm = TRUE),
    prop_decrease = mean(chg < 0, na.rm = TRUE)
  )
})

write_xlsx(direction_summary, "results/nhs_2nd_NULISA_direction.xlsx")

##-------------------------------
## Gene and NULISA change
##-------------------------------

gene_res_list <- list()
i <- 1

for (prot in c(protall_list, scoreall_list)) {
  
  message(prot)
  
  for (prot_type in c("_chg", "_chg_yr")) {
    
    for (gene in c("APOE4_carrier", "PGS002280_scaled_pool_rsq0.3", "PGS000334_scaled_pool_rsq0.3")) {
      
      tmp_covar <- c("age_1st", "bld_yr_diff", paste0(prot, "_1"), covar_gene_list)
      
      tmp_mod <- tryCatch(
        glm(
          as.formula(
            paste(
              paste0(prot, prot_type),
              " ~ ",
              paste(c(gene, tmp_covar), collapse = " + ")
            )
          ),
          data = model_data
        ),
        error = function(e) e,
        warning = function(w) w
      )
      
      gene_res_list[[i]] <- c(
        prot,
        gsub("^_", "", prot_type),
        gene,
        if (inherits(tmp_mod, "glm"))
          c(
            nobs(tmp_mod),
            summary(tmp_mod)$coefficients[gene, c(1, 2, 4)]
          )
        else c(tmp_mod$message, rep(NA, 3))
      )
      
      i <- i + 1
    }
  }
}

### Global change

for (gene in c("APOE4_carrier", "PGS002280_scaled_pool_rsq0.3", "PGS000334_scaled_pool_rsq0.3")) {
  
  tmp_covar <- c("age_1st", "bld_yr_diff", covar_gene_list)
  
  tmp_mod <- tryCatch(
    glm(
      as.formula(
        paste(
          "global_change ~",
          paste(c(gene, tmp_covar), collapse = " + ")
        )
      ),
      data = model_data
    ),
    error = function(e) e,
    warning = function(w) w
  )
  
  gene_res_list[[i]] <- c(
    "global_change",
    "chg",
    gene,
    if (inherits(tmp_mod, "glm"))
      c(
        nobs(tmp_mod),
        summary(tmp_mod)$coefficients[gene, c(1, 2, 4)]
      )
    else c(tmp_mod$message, rep(NA, 3))
  )
  
  i <- i + 1
}

gene_res_dat <- as.data.frame(do.call("rbind", gene_res_list))

colnames(gene_res_dat) <- c(
  "targetName",
  "prot_type",
  "gene",
  "n",
  "beta",
  "se",
  "p"
)

gene_res_dat <- gene_res_dat %>%
  mutate(
    n = as.numeric(n),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) %>% 
  left_join(marker_dict, by = "targetName")

write_xlsx(gene_res_dat, "results/nhs_2nd_NULISA_gene.xlsx")

##-------------------------------
## NULISA and dementia risk
##-------------------------------

dementia_res_list <- list()
i <- 1

for (prot in c(protall_list, scoreall_list)) {
  
  message(prot)
  
  for (prot_type in c("_2", "_chg", "_chg_yr")) {
    
    model_data_scaled <- model_data %>% 
      mutate(
        !!paste0(prot, "_1") := as.numeric(scale(.data[[paste0(prot, "_1")]])),
        !!paste0(prot, "_2") := as.numeric(scale(.data[[paste0(prot, "_2")]]))
      )
    
    if (prot_type != "_2") {
      model_data_scaled <- model_data_scaled %>% 
        mutate(
          !!paste0(prot, prot_type) := as.numeric(scale(.data[[paste0(prot, prot_type)]]))
        )
    }
    
    tmp_covar <- c("age_2nd", "fhdem", "dep_antidep")
    
    if (prot_type %in% c("_chg", "_chg_yr")) {
      tmp_covar <- c("bld_yr_diff", paste0(prot, "_1"), tmp_covar)
    }
    
    tmp_mod <- tryCatch(
      coxph(
        as.formula(
          paste(
            "Surv(tad_2nd, adcase) ~",
            paste(c(paste0(prot, prot_type), tmp_covar), collapse = " + ")
          )
        ),
        ties = "breslow",
        data = model_data
      ),
      error = function(e) e,
      warning = function(w) w
    )
    
    tmp_mod_scaled <- tryCatch(
      coxph(
        as.formula(
          paste(
            "Surv(tad_2nd, adcase) ~",
            paste(c(paste0(prot, prot_type), tmp_covar), collapse = " + ")
          )
        ),
        ties = "breslow",
        data = model_data_scaled
      ),
      error = function(e) e,
      warning = function(w) w
    )
    
    dementia_res_list[[i]] <- c(
      prot,
      gsub("^_", "", prot_type),
      "dementia",
      "raw",
      if (inherits(tmp_mod, "coxph")) {
        sm <- summary(tmp_mod)
        c(
          sm$n,
          sm$nevent,
          sm$coefficients[paste0(prot, prot_type), c(1, 3, 5)]
        )
      } else {
        c(tmp_mod$message, rep(NA, 4))
      }
    )
    
    i <- i + 1
    
    dementia_res_list[[i]] <- c(
      prot,
      gsub("^_", "", prot_type),
      "dementia",
      "scaled",
      if (inherits(tmp_mod_scaled, "coxph")) {
        sm <- summary(tmp_mod_scaled)
        c(
          sm$n,
          sm$nevent,
          sm$coefficients[paste0(prot, prot_type), c(1, 3, 5)]
        )
      } else {
        c(tmp_mod_scaled$message, rep(NA, 4))
      }
    )
    
    i <- i + 1
  }
}

dementia_res_dat <- as.data.frame(do.call("rbind", dementia_res_list))

colnames(dementia_res_dat) <- c(
  "targetName",
  "prot_type",
  "outcome",
  "scale",
  "n",
  "n_events",
  "beta",
  "se",
  "p"
)

dementia_res_dat <- dementia_res_dat %>% 
  mutate(
    n = as.numeric(n),
    n_events = as.numeric(n_events),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(dementia_res_dat, "results/nhs_2nd_dementia_NULISA.xlsx")