library(dplyr)
library(writexl)

##-----------------------------
## Age associations
##-----------------------------

rm(list = ls())

dir.create("results", recursive = TRUE, showWarnings = FALSE)

## bootstrap settings
nboot <- 2000
set.seed(20260321)

run_age_analysis <- function(model_data2, age_var, prot_list, marker_dict) {
  
  protall_list <- prot_list[
    !grepl("raw", prot_list) &
      !(prot_list %in% (marker_dict %>%
                          filter(LOD_flag) %>%
                          pull(targetName)))
  ]
  
  model_data2 <- model_data2 %>%
    dplyr::mutate(age10 = .data[[age_var]] / 10)
  
  age_res_list <- list()
  i <- 1
  
  for (prot in protall_list) {
    
    message(prot)
    
    tmp_dat <- model_data2 %>%
      dplyr::select(all_of(prot), age10) %>%
      dplyr::filter(
        !is.na(.data[[prot]]),
        !is.na(age10)
      )
    
    n_cc <- nrow(tmp_dat)
    
    if (n_cc < 40) {
      age_res_list[[i]] <- tibble(
        targetName = prot,
        n = n_cc,
        beta_age10 = NA_real_,
        se_age10 = NA_real_,
        p_age10 = NA_real_,
        beta_lcl = NA_real_,
        beta_ucl = NA_real_,
        r2 = NA_real_,
        r2_lcl = NA_real_,
        r2_ucl = NA_real_,
        status = "SKIP"
      )
      i <- i + 1
      next
    }
    
    fit <- tryCatch(
      lm(reformulate("age10", response = prot), data = tmp_dat),
      error = function(e) e
    )
    
    if (!inherits(fit, "lm")) {
      age_res_list[[i]] <- tibble(
        targetName = prot,
        n = n_cc,
        beta_age10 = NA_real_,
        se_age10 = NA_real_,
        p_age10 = NA_real_,
        beta_lcl = NA_real_,
        beta_ucl = NA_real_,
        r2 = NA_real_,
        r2_lcl = NA_real_,
        r2_ucl = NA_real_,
        status = paste0("MODEL ERROR: ", fit$message)
      )
      i <- i + 1
      next
    }
    
    sm <- summary(fit)
    
    ## bootstrap beta and R2
    boot_beta <- rep(NA_real_, nboot)
    boot_r2   <- rep(NA_real_, nboot)
    
    for (b in seq_len(nboot)) {
      idx <- sample(seq_len(n_cc), size = n_cc, replace = TRUE)
      boot_dat <- tmp_dat[idx, , drop = FALSE]
      
      boot_fit <- tryCatch(
        lm(reformulate("age10", response = prot), data = boot_dat),
        error = function(e) NULL
      )
      
      if (!is.null(boot_fit)) {
        boot_sm <- summary(boot_fit)
        if ("age10" %in% rownames(boot_sm$coefficients)) {
          boot_beta[b] <- boot_sm$coefficients["age10", "Estimate"]
          boot_r2[b]   <- boot_sm$r.squared
        }
      }
    }
    
    beta_ci <- quantile(boot_beta, probs = c(0.025, 0.975), na.rm = TRUE)
    r2_ci   <- quantile(boot_r2, probs = c(0.025, 0.975), na.rm = TRUE)
    
    age_res_list[[i]] <- tibble(
      targetName = prot,
      n = n_cc,
      beta_age10 = sm$coefficients["age10", "Estimate"],
      se_age10   = sm$coefficients["age10", "Std. Error"],
      p_age10    = sm$coefficients["age10", "Pr(>|t|)"],
      beta_lcl   = beta_ci[1],
      beta_ucl   = beta_ci[2],
      r2         = sm$r.squared,
      r2_lcl     = r2_ci[1],
      r2_ucl     = r2_ci[2],
      status     = "OK"
    )
    
    i <- i + 1
  }
  
  age_res_dat <- bind_rows(age_res_list) %>%
    dplyr::left_join(marker_dict, by = "targetName") %>%
    dplyr::arrange(dplyr::desc(r2))
  
  age_res_dat
}

##-----------------------------
## NHS prospective cohort
##-----------------------------

load("data/processed_data_nhs.RData")

age_res_dat <- run_age_analysis(
  model_data2 = model_data,
  age_var = "age_final",
  prot_list = prot_list,
  marker_dict = marker_dict
)

write_xlsx(age_res_dat, "results/nhs_age_NULISA_boot.xlsx")

rm(model_data, age_res_dat)

##-----------------------------
## NHS nested case-control
##-----------------------------

load("data/processed_data_nhs_caco.RData")

age_res_dat <- run_age_analysis(
  model_data2 = model_data,
  age_var = "age_final",
  prot_list = prot_list,
  marker_dict = marker_dict
)

write_xlsx(age_res_dat, "results/nhs_caco_age_NULISA_boot.xlsx")

age_res_dat <- run_age_analysis(
  model_data2 = model_data %>% filter(caco == 0),
  age_var = "age_final",
  prot_list = prot_list,
  marker_dict = marker_dict
)

write_xlsx(age_res_dat, "results/nhs_caco_age_NULISA_boot_ctrl.xlsx")

rm(model_data, age_res_dat)

##-----------------------------
## DIRECT-PLUS baseline
##-----------------------------

load("data/processed_data_dp.RData")

age_res_dat <- run_age_analysis(
  model_data2 = merged_data %>% filter(time == "A"),
  age_var = "age_baseline",
  prot_list = prot_list,
  marker_dict = marker_dict
)

write_xlsx(age_res_dat, "results/dp_age_NULISA_boot.xlsx")
