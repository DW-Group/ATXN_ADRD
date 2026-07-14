library(dplyr)
library(tidyverse)
library(survival)
library(readxl)
library(writexl)
library(Hmisc)

rm(list = ls())

##-------------------------------
## Load data
##-------------------------------

load("data/processed_data_nhs.RData")
dir.create("results", recursive = TRUE, showWarnings = FALSE)

##-------------------------------
## Correlations
##-------------------------------

cor_mat <- model_data %>%
  select(all_of(c(protall_list, scoreall_list))) %>%
  as.matrix()

rc <- rcorr(cor_mat, type = "spearman")

cor_long <- as.data.frame(as.table(rc$r)) %>%
  rename(var1 = Var1, var2 = Var2, r = Freq) %>%
  left_join(
    as.data.frame(as.table(rc$P)) %>%
      rename(var1 = Var1, var2 = Var2, p = Freq),
    by = c("var1", "var2")
  )

write_xlsx(cor_long, "results/nhs_NULISA_spearman_corr.xlsx")

##-------------------------------
## Diet and NULISA
##-------------------------------

load("data/NULISA_core.RData")

diet_marker_list <- unique(c(
  "all_core3_ew_signed",
  core_list[["all"]][["core3"]],
  core_list[["BD"]][["core3"]]
))

diet_marker_list <- diet_marker_list[diet_marker_list %in% colnames(model_data)]

diet_var_list <- unique(c(
  dietscore_list,
  amedcomp_list,
  mindcomp_list
))

diet_var_list <- diet_var_list[diet_var_list %in% colnames(model_data)]

diet_res_list <- list()
i <- 1

for (prot in diet_marker_list) {
  
  message(prot)
  
  model_data_scaled <- model_data %>%
    mutate(!!prot := as.numeric(scale(.data[[prot]])))
  
  for (diet in diet_var_list) {
    
    tmp_mod_scaled <- tryCatch(
      glm(
        as.formula(
          paste(
            prot,
            " ~ ",
            paste(c(diet, "age_final", covar_diet_list), collapse = " + ")
          )
        ),
        data = model_data_scaled
      ),
      error = function(e) e,
      warning = function(w) w
    )
    
    diet_res_list[[i]] <- c(
      prot,
      diet,
      "univariate",
      "scaled",
      if (inherits(tmp_mod_scaled, "glm")) {
        c(
          nobs(tmp_mod_scaled),
          summary(tmp_mod_scaled)$coefficients[diet, c(1, 2, 4)]
        )
      } else {
        c(tmp_mod_scaled$message, rep(NA, 3))
      }
    )
    
    i <- i + 1
  }
}

diet_res_dat <- as.data.frame(do.call("rbind", diet_res_list))

colnames(diet_res_dat) <- c(
  "targetName",
  "diet",
  "model",
  "scale",
  "n",
  "beta",
  "se",
  "p"
)

diet_res_dat <- diet_res_dat %>%
  mutate(
    n = suppressWarnings(as.numeric(n)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    p = suppressWarnings(as.numeric(p))
  ) %>% 
  left_join(marker_dict, by = "targetName")

write_xlsx(diet_res_dat, "results/nhs_NULISA_diet.xlsx")

##-------------------------------
## NULISA and dementia risk
##-------------------------------

dementia_marker_list <- unique(c(protall_list, scoreall_list))
dementia_marker_list <- dementia_marker_list[dementia_marker_list %in% colnames(model_data)]

dementia_res_list <- list()
i <- 1

for (prot in dementia_marker_list) {
  
  message(prot)
  
  model_data_scaled <- model_data %>%
    mutate(!!prot := as.numeric(scale(.data[[prot]])))
  
  tmp_mod_scaled <- tryCatch(
    coxph(
      as.formula(
        paste(
          "Surv(tad_final, adcase) ~",
          paste(c(prot, covar_list), collapse = " + ")
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
    "dementia",
    "scaled",
    if (inherits(tmp_mod_scaled, "coxph")) {
      sm <- summary(tmp_mod_scaled)
      c(
        sm$n,
        sm$nevent,
        sm$coefficients[prot, c(1, 3, 5)]
      )
    } else {
      c(tmp_mod_scaled$message, rep(NA, 4))
    }
  )
  
  i <- i + 1
}

dementia_res_dat <- as.data.frame(do.call("rbind", dementia_res_list))

colnames(dementia_res_dat) <- c(
  "targetName",
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
    n = suppressWarnings(as.numeric(n)),
    n_events = suppressWarnings(as.numeric(n_events)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    p = suppressWarnings(as.numeric(p)),
    HR = exp(beta),
    p_fdr = p.adjust(p, method = "BH")
  ) %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(dementia_res_dat, "results/nhs_dementia_NULISA.xlsx")

##-------------------------------
## NULISA and dementia risk quartiles
##-------------------------------

dementia_quartile_res_list <- list()
i <- 1

for (prot in dementia_marker_list) {
  
  message(prot)
  
  dat <- model_data %>%
    dplyr::select(dplyr::all_of(c("tad_final", "adcase", prot, covar_list))) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.))) %>%
    dplyr::mutate(
      prot_quartile_num = dplyr::ntile(.data[[prot]], 4),
      prot_quartile = factor(
        paste0("Q", prot_quartile_num),
        levels = c("Q1", "Q2", "Q3", "Q4")
      )
    ) %>%
    dplyr::group_by(prot_quartile) %>%
    dplyr::mutate(
      prot_quartile_median = median(.data[[prot]], na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
  
  quartile_counts <- dat %>%
    dplyr::group_by(prot_quartile) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_events = sum(adcase == 1, na.rm = TRUE),
      marker_median = median(.data[[prot]], na.rm = TRUE),
      .groups = "drop"
    )
  
  tmp_mod_cat <- tryCatch(
    coxph(
      as.formula(
        paste(
          "Surv(tad_final, adcase) ~ prot_quartile +",
          paste(covar_list, collapse = " + ")
        )
      ),
      ties = "breslow",
      data = dat
    ),
    error = function(e) e
  )
  
  tmp_mod_trend <- tryCatch(
    coxph(
      as.formula(
        paste(
          "Surv(tad_final, adcase) ~ prot_quartile_median +",
          paste(covar_list, collapse = " + ")
        )
      ),
      ties = "breslow",
      data = dat
    ),
    error = function(e) e
  )
  
  p_trend <- if (inherits(tmp_mod_trend, "coxph")) {
    summary(tmp_mod_trend)$coefficients["prot_quartile_median", "Pr(>|z|)"]
  } else {
    NA_real_
  }
  
  p_nonlinearity <- if (inherits(tmp_mod_cat, "coxph") && inherits(tmp_mod_trend, "coxph")) {
    tryCatch(
      anova(tmp_mod_trend, tmp_mod_cat, test = "LRT")[2, "Pr(>|Chi|)"],
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  
  if (inherits(tmp_mod_cat, "coxph")) {
    
    sm <- summary(tmp_mod_cat)$coefficients
    ci <- summary(tmp_mod_cat)$conf.int
    
    for (quartile in c("Q1", "Q2", "Q3", "Q4")) {
      
      tmp_count <- quartile_counts %>%
        dplyr::filter(prot_quartile == quartile)
      
      if (quartile == "Q1") {
        
        dementia_quartile_res_list[[i]] <- c(
          prot,
          "dementia",
          quartile,
          "Q1",
          tmp_count$n,
          tmp_count$n_events,
          tmp_count$marker_median,
          0,
          NA,
          1,
          NA,
          NA,
          NA,
          p_trend,
          p_nonlinearity
        )
        
      } else {
        
        coef_name <- paste0("prot_quartile", quartile)
        
        dementia_quartile_res_list[[i]] <- c(
          prot,
          "dementia",
          quartile,
          "Q1",
          tmp_count$n,
          tmp_count$n_events,
          tmp_count$marker_median,
          sm[coef_name, "coef"],
          sm[coef_name, "se(coef)"],
          ci[coef_name, "exp(coef)"],
          ci[coef_name, "lower .95"],
          ci[coef_name, "upper .95"],
          sm[coef_name, "Pr(>|z|)"],
          p_trend,
          p_nonlinearity
        )
      }
      
      i <- i + 1
    }
    
  } else {
    
    for (quartile in c("Q1", "Q2", "Q3", "Q4")) {
      
      tmp_count <- quartile_counts %>%
        dplyr::filter(prot_quartile == quartile)
      
      dementia_quartile_res_list[[i]] <- c(
        prot,
        "dementia",
        quartile,
        "Q1",
        tmp_count$n,
        tmp_count$n_events,
        tmp_count$marker_median,
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        p_trend,
        p_nonlinearity
      )
      
      i <- i + 1
    }
  }
}

dementia_quartile_long_dat <- as.data.frame(do.call("rbind", dementia_quartile_res_list))

colnames(dementia_quartile_long_dat) <- c(
  "targetName",
  "outcome",
  "quartile",
  "reference",
  "n_quartile",
  "n_events_quartile",
  "marker_median_quartile",
  "quartile_beta",
  "quartile_se",
  "quartile_HR",
  "quartile_HR_L",
  "quartile_HR_U",
  "quartile_p",
  "p_trend_quartile",
  "p_nonlinearity_quartile"
)

dementia_quartile_long_dat <- dementia_quartile_long_dat %>%
  dplyr::mutate(
    n_quartile = as.numeric(n_quartile),
    n_events_quartile = as.numeric(n_events_quartile),
    marker_median_quartile = as.numeric(marker_median_quartile),
    quartile_beta = as.numeric(quartile_beta),
    quartile_se = as.numeric(quartile_se),
    quartile_HR = as.numeric(quartile_HR),
    quartile_HR_L = as.numeric(quartile_HR_L),
    quartile_HR_U = as.numeric(quartile_HR_U),
    quartile_p = as.numeric(quartile_p),
    p_trend_quartile = as.numeric(p_trend_quartile),
    p_nonlinearity_quartile = as.numeric(p_nonlinearity_quartile)
  )

dementia_quartile_wide_dat <- dementia_quartile_long_dat %>%
  dplyr::filter(quartile %in% c("Q2", "Q3", "Q4")) %>%
  dplyr::select(
    targetName,
    quartile,
    n_quartile,
    n_events_quartile,
    marker_median_quartile,
    quartile_beta,
    quartile_se,
    quartile_HR,
    quartile_HR_L,
    quartile_HR_U,
    quartile_p
  ) %>%
  tidyr::pivot_wider(
    id_cols = targetName,
    names_from = quartile,
    values_from = c(
      n_quartile,
      n_events_quartile,
      marker_median_quartile,
      quartile_beta,
      quartile_se,
      quartile_HR,
      quartile_HR_L,
      quartile_HR_U,
      quartile_p
    ),
    names_glue = "{quartile}_vs_Q1_{.value}"
  ) %>%
  dplyr::left_join(
    dementia_quartile_long_dat %>%
      dplyr::distinct(targetName, p_trend_quartile, p_nonlinearity_quartile),
    by = "targetName"
  ) %>%
  dplyr::left_join(marker_dict, by = "targetName")

write_xlsx(
  dementia_quartile_wide_dat,
  "results/nhs_dementia_NULISA_quartile.xlsx"
)

##-------------------------------
## NULISA and TICS
##-------------------------------

tics_marker_list <- unique(c(protall_list, scoreall_list))
tics_marker_list <- tics_marker_list[tics_marker_list %in% colnames(model_data)]

tics_res_list <- list()
i <- 1

for (prot in tics_marker_list) {
  
  message(prot)
  
  model_data_scaled <- model_data %>%
    mutate(!!prot := as.numeric(scale(.data[[prot]])))
  
  for (tics in tics_list) {
    
    tmp_mod_scaled <- tryCatch(
      glm(
        reformulate(c(prot, covar_list), response = tics),
        data = model_data_scaled
      ),
      error = function(e) e,
      warning = function(w) w
    )
    
    tics_res_list[[i]] <- c(
      prot,
      tics,
      "scaled",
      if (inherits(tmp_mod_scaled, "glm")) {
        c(
          sum(!is.na(model_data_scaled[[tics]])),
          nobs(tmp_mod_scaled),
          summary(tmp_mod_scaled)$coefficients[prot, c(1, 2, 4)]
        )
      } else {
        c(tmp_mod_scaled$message, rep(NA, 4))
      }
    )
    
    i <- i + 1
  }
}

tics_res_dat <- as.data.frame(do.call("rbind", tics_res_list))

colnames(tics_res_dat) <- c(
  "targetName",
  "outcome",
  "scale",
  "n",
  "n_model",
  "beta",
  "se",
  "p"
)

tics_res_dat <- tics_res_dat %>% 
  mutate(
    n = suppressWarnings(as.numeric(n)),
    n_model = suppressWarnings(as.numeric(n_model)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    p = suppressWarnings(as.numeric(p))
  ) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(tics_res_dat, "results/nhs_tics_NULISA.xlsx")

##-------------------------------
## R2
##-------------------------------

load("data/NULISA_core.RData")

core_vec <- unique(c(
  core_list[["all"]][["core3"]],
  core_list[["BD"]][["core3"]]
))

joint_name <- "core_joint"

tics_res_list <- list()
i <- 1

for (prot in c(core_vec, joint_name)) {
  
  message(prot)
  
  for (tics in tics_list) {
    
    prot_vars <- if (prot == joint_name) core_list[["all"]][["core3"]] else prot
    
    dat <- model_data %>%
      dplyr::select(dplyr::all_of(c(tics, prot_vars, covar_list))) %>%
      dplyr::filter(if_all(everything(), ~ !is.na(.)))
    
    n_total <- sum(!is.na(model_data[[tics]]))
    n_used  <- nrow(dat)
    
    base_form <- reformulate(covar_list, response = tics)
    full_form <- reformulate(c(prot_vars, covar_list), response = tics)
    
    base_mod <- tryCatch(lm(base_form, data = dat), error = function(e) e)
    full_mod <- tryCatch(lm(full_form, data = dat), error = function(e) e)
    
    if (inherits(base_mod, "lm") && inherits(full_mod, "lm")) {
      
      sm_base <- summary(base_mod)
      sm_full <- summary(full_mod)
      
      r2_base  <- sm_base$r.squared
      r2_full  <- sm_full$r.squared
      delta_r2 <- r2_full - r2_base
      
      adj_r2_base  <- sm_base$adj.r.squared
      adj_r2_full  <- sm_full$adj.r.squared
      delta_adj_r2 <- adj_r2_full - adj_r2_base
      
    } else {
      
      r2_base <- NA_real_
      r2_full <- NA_real_
      delta_r2 <- NA_real_
      adj_r2_base <- NA_real_
      adj_r2_full <- NA_real_
      delta_adj_r2 <- NA_real_
    }
    
    tics_res_list[[i]] <- data.frame(
      targetName = prot,
      outcome = tics,
      n = n_total,
      n_used = n_used,
      r2_base = r2_base,
      r2_full = r2_full,
      delta_r2 = delta_r2,
      adj_r2_base = adj_r2_base,
      adj_r2_full = adj_r2_full,
      delta_adj_r2 = delta_adj_r2,
      note = if (prot == joint_name) paste(core_list[["all"]][["core3"]], collapse = " + ") else NA_character_
    )
    
    i <- i + 1
  }
}

tics_res_dat <- dplyr::bind_rows(tics_res_list)

write_xlsx(
  tics_res_dat,
  "results/nhs_tics_NULISA_R2.xlsx"
)