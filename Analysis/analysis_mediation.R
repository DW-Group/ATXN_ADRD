library(dplyr)
library(mediation)
library(writexl)

##---------------------------------------------------------
## NHS 
##---------------------------------------------------------

rm(list = ls())

load("data/processed_data_nhs.RData")
load("data/NULISA_core.RData")

## diet exposure 
diet_var_list <- c("MIND90_cumavg", mindcomp_list)

diet_var_list <- diet_var_list[diet_var_list %in% colnames(model_data)]

## outcomes and mediators  
tics_list <- "av_ntotal"

protall_list <- unique(c(
  "all_core3_ew_signed",
  core_list[["all"]][["core3"]],
  core_list[["BD"]][["core3"]]
))

protall_list <- protall_list[protall_list %in% colnames(model_data)]

## covariates
med_covar_list <- unique(c(covar_list, covar_diet_list))

## bootstrap reps
nsim_boot <- 2000

set.seed(20260210)

med_res_list <- list()
i <- 1

for (prot in protall_list) {
  
  message("Mediator: ", prot)
  
  for (tics in tics_list) {
    
    for (diet_var in diet_var_list) {
      
      vars_needed <- unique(c(
        tics,
        diet_var,
        prot,
        "age_final",
        covar_diet_list,
        med_covar_list
      ))
      
      tmp_dat <- model_data %>%
        dplyr::select(all_of(vars_needed)) %>%
        dplyr::filter(
          !is.na(.data[[tics]]),
          !is.na(.data[[diet_var]]),
          !is.na(.data[[prot]])
        )
      
      if (length(med_covar_list) > 0) {
        for (cv in med_covar_list) {
          tmp_dat <- tmp_dat %>% 
            dplyr::filter(!is.na(.data[[cv]]))
        }
      }
      
      n_cc <- nrow(tmp_dat)
      
      if (n_cc < 40) {
        med_res_list[[i]] <- c(
          targetName = prot,
          outcome    = tics,
          diet       = diet_var,
          n          = n_cc,
          status     = "SKIP: too few complete cases",
          a_beta = NA, a_se = NA, a_p = NA,
          b_beta = NA, b_se = NA, b_p = NA,
          cprime_beta = NA, cprime_se = NA, cprime_p = NA,
          acme       = NA, acme_lcl = NA, acme_ucl = NA, acme_p = NA,
          ade        = NA, ade_lcl  = NA, ade_ucl  = NA, ade_p  = NA,
          te         = NA, te_lcl   = NA, te_ucl   = NA, te_p   = NA,
          pm         = NA, pm_lcl   = NA, pm_ucl   = NA, pm_p   = NA
        )
        i <- i + 1
        next
      }
      
      tmp_dat <- tmp_dat %>%
        mutate(
          diet_z = as.numeric(scale(.data[[diet_var]])),
          prot_z = as.numeric(scale(.data[[prot]]))
        )
      
      med_formula <- reformulate(
        c("diet_z", "age_final", covar_diet_list),
        response = "prot_z"
      )
      
      med_fit <- tryCatch(
        suppressWarnings(lm(med_formula, data = tmp_dat)),
        error = function(e) e
      )
      
      out_formula <- reformulate(
        c("diet_z", "prot_z", med_covar_list),
        response = tics
      )
      
      out_fit <- tryCatch(
        suppressWarnings(lm(out_formula, data = tmp_dat)),
        error = function(e) e
      )
      
      if (!inherits(med_fit, "lm") || !inherits(out_fit, "lm")) {
        
        msg <- if (!inherits(med_fit, "lm")) med_fit$message else out_fit$message
        
        med_res_list[[i]] <- c(
          targetName = prot,
          outcome    = tics,
          diet       = diet_var,
          n          = n_cc,
          status     = paste0("MODEL ERROR: ", msg),
          a_beta = NA, a_se = NA, a_p = NA,
          b_beta = NA, b_se = NA, b_p = NA,
          cprime_beta = NA, cprime_se = NA, cprime_p = NA,
          acme       = NA, acme_lcl = NA, acme_ucl = NA, acme_p = NA,
          ade        = NA, ade_lcl  = NA, ade_ucl  = NA, ade_p  = NA,
          te         = NA, te_lcl   = NA, te_ucl   = NA, te_p   = NA,
          pm         = NA, pm_lcl   = NA, pm_ucl   = NA, pm_p   = NA
        )
        i <- i + 1
        next
      }
      
      med_obj <- tryCatch(
        suppressWarnings(
          mediate(
            model.m = med_fit,
            model.y = out_fit,
            treat   = "diet_z",
            mediator = "prot_z",
            control.value = -0.5,
            treat.value   =  0.5,
            boot = TRUE,
            sims = nsim_boot
          )
        ),
        error = function(e) e
      )
      
      if (!inherits(med_obj, "mediate")) {
        
        med_res_list[[i]] <- c(
          targetName = prot,
          outcome    = tics,
          diet       = diet_var,
          n          = n_cc,
          status     = paste0("MEDIATE ERROR: ", med_obj$message),
          a_beta = NA, a_se = NA, a_p = NA,
          b_beta = NA, b_se = NA, b_p = NA,
          cprime_beta = NA, cprime_se = NA, cprime_p = NA,
          acme       = NA, acme_lcl = NA, acme_ucl = NA, acme_p = NA,
          ade        = NA, ade_lcl  = NA, ade_ucl  = NA, ade_p  = NA,
          te         = NA, te_lcl   = NA, te_ucl   = NA, te_p   = NA,
          pm         = NA, pm_lcl   = NA, pm_ucl   = NA, pm_p   = NA
        )
        i <- i + 1
        next
      }
      
      s <- summary(med_obj)
      med_coef <- summary(med_fit)$coefficients
      out_coef <- summary(out_fit)$coefficients
      
      med_res_list[[i]] <- c(
        targetName = prot,
        outcome    = tics,
        diet       = diet_var,
        n          = n_cc,
        status     = "OK",
        
        a_beta = unname(med_coef["diet_z", "Estimate"]),
        a_se   = unname(med_coef["diet_z", "Std. Error"]),
        a_p    = unname(med_coef["diet_z", "Pr(>|t|)"]),
        
        b_beta = unname(out_coef["prot_z", "Estimate"]),
        b_se   = unname(out_coef["prot_z", "Std. Error"]),
        b_p    = unname(out_coef["prot_z", "Pr(>|t|)"]),
        
        cprime_beta = unname(out_coef["diet_z", "Estimate"]),
        cprime_se   = unname(out_coef["diet_z", "Std. Error"]),
        cprime_p    = unname(out_coef["diet_z", "Pr(>|t|)"]),
        
        acme       = unname(s$d.avg),
        acme_lcl   = unname(s$d.avg.ci[1]),
        acme_ucl   = unname(s$d.avg.ci[2]),
        acme_p     = unname(s$d.avg.p),
        
        ade        = unname(s$z.avg),
        ade_lcl    = unname(s$z.avg.ci[1]),
        ade_ucl    = unname(s$z.avg.ci[2]),
        ade_p      = unname(s$z.avg.p),
        
        te         = unname(s$tau.coef),
        te_lcl     = unname(s$tau.ci[1]),
        te_ucl     = unname(s$tau.ci[2]),
        te_p       = unname(s$tau.p),
        
        pm         = unname(s$n.avg),
        pm_lcl     = unname(s$n.avg.ci[1]),
        pm_ucl     = unname(s$n.avg.ci[2]),
        pm_p       = unname(s$n.avg.p)
      )
      
      i <- i + 1
    }
  }
}

med_res_dat <- as.data.frame(do.call("rbind", med_res_list), stringsAsFactors = FALSE)

num_cols <- c(
  "n",
  "a_beta", "a_se", "a_p",
  "b_beta", "b_se", "b_p",
  "cprime_beta", "cprime_se", "cprime_p",
  "acme", "acme_lcl", "acme_ucl", "acme_p",
  "ade", "ade_lcl", "ade_ucl", "ade_p",
  "te", "te_lcl", "te_ucl", "te_p",
  "pm", "pm_lcl", "pm_ucl", "pm_p"
)

for (cc in intersect(num_cols, colnames(med_res_dat))) {
  med_res_dat[[cc]] <- suppressWarnings(as.numeric(med_res_dat[[cc]]))
}

med_res_dat <- med_res_dat %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(
  med_res_dat,
  "results/nhs_tics_mediation.xlsx"
)