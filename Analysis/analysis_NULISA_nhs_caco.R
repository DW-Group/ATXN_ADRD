library(dplyr)
library(tidyverse)
library(data.table)
library(survival)
library(writexl)
library(pROC)
library(fastshap)

rm(list = ls())

##-------------------------------
## Load data
##-------------------------------

load("data/processed_data_nhs_caco.RData")

dir.create("results", recursive = TRUE, showWarnings = FALSE)

##-------------------------------
## Diet and NULISA
##-------------------------------

diet_res_list <- list()
i <- 1

analysis_data_list <- list(
  all = model_data,
  control_only = model_data %>% filter(adcase == 0)
)

for (analysis_set in names(analysis_data_list)) {
  
  message("Analysis set: ", analysis_set)
  
  tmp_model_data_raw <- analysis_data_list[[analysis_set]]
  
  for (prot in c(protall_list, scoreall_list)) {
    
    message(prot)
    
    tmp_model_data_scaled <- tmp_model_data_raw %>%
      mutate(!!prot := as.numeric(scale(.data[[prot]])))
    
    for (scale_type in c("raw", "scaled")) {
      
      tmp_model_data <- if (scale_type == "raw") tmp_model_data_raw else tmp_model_data_scaled
      
      for (diet in c(dietscore_list, amedcomp_list, mindcomp_list)) {
        
        tmp_mod <- tryCatch(
          glm(
            as.formula(
              paste(
                prot,
                " ~ ",
                paste(c(diet, "age_final", covar_diet_list), collapse = " + ")
              )
            ),
            data = tmp_model_data
          ),
          error = function(e) e,
          warning = function(w) w
        )
        
        diet_res_list[[i]] <- c(
          analysis_set,
          prot,
          diet,
          "univariate",
          scale_type,
          if (inherits(tmp_mod, "glm")) {
            c(
              nobs(tmp_mod),
              summary(tmp_mod)$coefficients[diet, c(1, 2, 4)]
            )
          } else {
            c(tmp_mod$message, rep(NA, 3))
          }
        )
        
        i <- i + 1
      }
      
      ## Mutual aMED components only
      
      tmp_mod_mut_amed <- tryCatch(
        glm(
          as.formula(
            paste(
              prot,
              " ~ ",
              paste(c(amedcomp_list, "age_final", covar_diet_list), collapse = " + ")
            )
          ),
          data = tmp_model_data
        ),
        error = function(e) e,
        warning = function(w) w
      )
      
      if (inherits(tmp_mod_mut_amed, "glm")) {
        sm <- summary(tmp_mod_mut_amed)$coefficients
        
        for (diet in amedcomp_list) {
          diet_res_list[[i]] <- c(
            analysis_set,
            prot,
            diet,
            "mutual_amed",
            scale_type,
            c(
              nobs(tmp_mod_mut_amed),
              if (diet %in% rownames(sm)) sm[diet, c(1, 2, 4)] else c(NA, NA, NA)
            )
          )
          
          i <- i + 1
        }
      } else {
        for (diet in amedcomp_list) {
          diet_res_list[[i]] <- c(
            analysis_set,
            prot,
            diet,
            "mutual_amed",
            scale_type,
            c(tmp_mod_mut_amed$message, rep(NA, 3))
          )
          
          i <- i + 1
        }
      }
      
      ## Mutual MIND components only
      
      tmp_mod_mut_mind <- tryCatch(
        glm(
          as.formula(
            paste(
              prot,
              " ~ ",
              paste(c(mindcomp_list, "age_final", covar_diet_list), collapse = " + ")
            )
          ),
          data = tmp_model_data
        ),
        error = function(e) e,
        warning = function(w) w
      )
      
      if (inherits(tmp_mod_mut_mind, "glm")) {
        sm <- summary(tmp_mod_mut_mind)$coefficients
        
        for (diet in mindcomp_list) {
          diet_res_list[[i]] <- c(
            analysis_set,
            prot,
            diet,
            "mutual_mind",
            scale_type,
            c(
              nobs(tmp_mod_mut_mind),
              if (diet %in% rownames(sm)) sm[diet, c(1, 2, 4)] else c(NA, NA, NA)
            )
          )
          
          i <- i + 1
        }
      } else {
        for (diet in mindcomp_list) {
          diet_res_list[[i]] <- c(
            analysis_set,
            prot,
            diet,
            "mutual_mind",
            scale_type,
            c(tmp_mod_mut_mind$message, rep(NA, 3))
          )
          
          i <- i + 1
        }
      }
    }
  }
}

diet_res_dat <- as.data.frame(do.call("rbind", diet_res_list))

colnames(diet_res_dat) <- c(
  "analysis_set",
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
    n = as.numeric(n),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(
  diet_res_dat,
  "results/nhs_caco_NULISA_diet.xlsx"
)

##-------------------------------
## NULISA and dementia risk
##-------------------------------

dementia_res_list <- list()
i <- 1

for (prot in c(protall_list, scoreall_list)) {
  
  message(prot)
  
  model_data_scaled <- model_data %>%
    mutate(!!prot := as.numeric(scale(.data[[prot]])))
  
  tmp_mod <- tryCatch(
    clogit(
      as.formula(
        paste(
          "caco ~",
          paste(c(prot, covar_list, "strata(caseid)"), collapse = " + ")
        )
      ),
      data = model_data
    ),
    error = function(e) e,
    warning = function(w) w
  )
  
  tmp_mod_scaled <- tryCatch(
    clogit(
      as.formula(
        paste(
          "caco ~",
          paste(c(prot, covar_list, "strata(caseid)"), collapse = " + ")
        )
      ),
      data = model_data_scaled
    ),
    error = function(e) e,
    warning = function(w) w
  )
  
  dementia_res_list[[i]] <- c(
    prot,
    "dementia",
    "raw",
    if (inherits(tmp_mod, "coxph")) {
      c(
        tmp_mod$n,
        tmp_mod$nevent,
        summary(tmp_mod)$coefficients[prot, c(1, 3, 5)]
      )
    } else {
      c(tmp_mod$message, rep(NA, 4))
    }
  )
  
  i <- i + 1
  
  dementia_res_list[[i]] <- c(
    prot,
    "dementia",
    "scaled",
    if (inherits(tmp_mod_scaled, "coxph")) {
      c(
        tmp_mod_scaled$n,
        tmp_mod_scaled$nevent,
        summary(tmp_mod_scaled)$coefficients[prot, c(1, 3, 5)]
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
    n = as.numeric(n),
    n_events = as.numeric(n_events),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(
  dementia_res_dat,
  "results/nhs_caco_dementia_NULISA.xlsx"
)

##-------------------------------
## Prediction
##-------------------------------

load("data/NULISA_core.RData")

predictor_list <- list(
  "m1" = c(
    "age_final"
  ),
  "m2" = c(
    "age_final",
    "fhdem",
    "dep_antidep"
  ),
  "m3" = c(
    "age_final",
    "fhdem",
    "dep_antidep",
    "APOE4_carrier",
    "PGS002280_scaled_pool_rsq0.3"
  ),
  "m13" = c(
    "age_final",
    "fhdem",
    "dep_antidep",
    "APOE4_carrier",
    "PGS002280_scaled_pool_rsq0.3",
    core_list[["all"]][["core3"]]
  ),
  "m5" = c(
    "age_final",
    "fhdem",
    "dep_antidep",
    "APOE4_carrier",
    "PGS002280_scaled_pool_rsq0.3",
    "pTau_217"
  ),
  "m6" = c(
    "age_final",
    "fhdem",
    "dep_antidep",
    "APOE4_carrier",
    "PGS002280_scaled_pool_rsq0.3",
    "BD_pTau_217"
  ),
  "m7" = c(
    "age_final",
    "fhdem",
    "dep_antidep",
    "APOE4_carrier",
    "PGS002280_scaled_pool_rsq0.3",
    core_list[["BD"]][["core3"]]
  ),
  "m8" = c(
    "age_final",
    "fhdem",
    "dep_antidep",
    core_list[["all"]][["core3"]]
  )
)

compare_pairs <- list(
  c("m1", "m2"),
  c("m2", "m3"),
  c("m3", "m13"),
  c("m3", "m5"),
  c("m3", "m6"),
  c("m3", "m7"),
  c("m2", "m8")
)

needed_models <- unique(c(unlist(compare_pairs), "m13"))

prediction_vars <- unique(c(
  "caco",
  unlist(predictor_list[needed_models])
))

prediction_data <- model_data %>%
  dplyr::select(all_of(prediction_vars)) %>%
  tidyr::drop_na()

fit_glm_binomial <- function(model_name, data, predictor_list) {
  
  preds <- predictor_list[[model_name]]
  
  glm(
    as.formula(
      paste(
        "caco ~",
        paste(preds, collapse = " + ")
      )
    ),
    family = binomial(),
    data = data
  )
}

calc_idi_cnri <- function(y, p0, p1) {
  
  stopifnot(length(y) == length(p0), length(y) == length(p1))
  stopifnot(all(y %in% c(0, 1)))
  
  case_idx <- y == 1
  ctrl_idx <- y == 0
  
  dp <- p1 - p0
  
  idi <- mean(dp[case_idx], na.rm = TRUE) -
    mean(dp[ctrl_idx], na.rm = TRUE)
  
  up_case   <- mean(dp[case_idx] > 0, na.rm = TRUE)
  down_case <- mean(dp[case_idx] < 0, na.rm = TRUE)
  down_ctrl <- mean(dp[ctrl_idx] < 0, na.rm = TRUE)
  up_ctrl   <- mean(dp[ctrl_idx] > 0, na.rm = TRUE)
  
  cnri <- (up_case - down_case) + (down_ctrl - up_ctrl)
  
  data.frame(
    idi = idi,
    cnri = cnri,
    mean_dp_case = mean(dp[case_idx], na.rm = TRUE),
    mean_dp_ctrl = mean(dp[ctrl_idx], na.rm = TRUE),
    up_case = up_case,
    down_case = down_case,
    down_ctrl = down_ctrl,
    up_ctrl = up_ctrl
  )
}

boot_pvalue <- function(x) {
  
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  
  2 * min(mean(x <= 0), mean(x >= 0))
}

boot_ci <- function(x, probs = c(0.025, 0.975)) {
  
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA_real_, NA_real_))
  
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

bootstrap_compare_models <- function(
    data,
    base_model,
    new_model,
    predictor_list,
    B = 2000,
    seed = 123
) {
  
  set.seed(seed)
  
  n <- nrow(data)
  out <- matrix(NA_real_, nrow = B, ncol = 3)
  colnames(out) <- c("delta_auc", "idi", "cnri")
  
  b_ok <- 0L
  
  for (b in seq_len(B)) {
    
    idx <- sample.int(n, size = n, replace = TRUE)
    d <- data[idx, , drop = FALSE]
    
    if (length(unique(d$caco)) < 2) next
    
    fit0 <- tryCatch(
      fit_glm_binomial(base_model, d, predictor_list),
      error = function(e) NULL
    )
    
    fit1 <- tryCatch(
      fit_glm_binomial(new_model, d, predictor_list),
      error = function(e) NULL
    )
    
    if (is.null(fit0) || is.null(fit1)) next
    
    p0 <- tryCatch(
      predict(fit0, type = "response"),
      error = function(e) NULL
    )
    
    p1 <- tryCatch(
      predict(fit1, type = "response"),
      error = function(e) NULL
    )
    
    if (is.null(p0) || is.null(p1)) next
    
    auc0 <- tryCatch(
      as.numeric(pROC::auc(d$caco, p0, quiet = TRUE)),
      error = function(e) NA_real_
    )
    
    auc1 <- tryCatch(
      as.numeric(pROC::auc(d$caco, p1, quiet = TRUE)),
      error = function(e) NA_real_
    )
    
    met <- tryCatch(
      calc_idi_cnri(d$caco, p0, p1),
      error = function(e) NULL
    )
    
    if (!is.finite(auc0) || !is.finite(auc1) || is.null(met)) next
    
    b_ok <- b_ok + 1L
    
    out[b_ok, "delta_auc"] <- auc1 - auc0
    out[b_ok, "idi"] <- met$idi
    out[b_ok, "cnri"] <- met$cnri
  }
  
  out <- as.data.frame(out[seq_len(b_ok), , drop = FALSE])
  
  list(
    boot = out,
    n_success = b_ok,
    delta_auc_ci = boot_ci(out$delta_auc),
    idi_ci = boot_ci(out$idi),
    cnri_ci = boot_ci(out$cnri),
    p_delta_auc = boot_pvalue(out$delta_auc),
    p_idi = boot_pvalue(out$idi),
    p_cnri = boot_pvalue(out$cnri)
  )
}

logistic_mod_list <- setNames(vector("list", length(needed_models)), needed_models)
pred_prob_list    <- setNames(vector("list", length(needed_models)), needed_models)
roc_list          <- setNames(vector("list", length(needed_models)), needed_models)

for (m in needed_models) {
  
  mod <- fit_glm_binomial(m, prediction_data, predictor_list)
  phat <- predict(mod, type = "response")
  
  logistic_mod_list[[m]] <- mod
  pred_prob_list[[m]] <- phat
  
  roc_list[[m]] <- roc(
    response = prediction_data$caco,
    predictor = phat,
    quiet = TRUE,
    direction = "auto"
  )
}

compare_two_models <- function(
    base_model,
    new_model,
    logistic_mod_list,
    pred_prob_list,
    roc_list,
    data,
    predictor_list,
    B = 2000,
    seed = 123
) {
  
  y <- data$caco
  
  mod0 <- logistic_mod_list[[base_model]]
  mod1 <- logistic_mod_list[[new_model]]
  
  p0 <- pred_prob_list[[base_model]]
  p1 <- pred_prob_list[[new_model]]
  
  roc0 <- roc_list[[base_model]]
  roc1 <- roc_list[[new_model]]
  
  auc0 <- as.numeric(auc(roc0))
  auc1 <- as.numeric(auc(roc1))
  
  ci0 <- ci.auc(roc0, method = "delong")
  ci1 <- ci.auc(roc1, method = "delong")
  
  delong_obj <- roc.test(
    roc0,
    roc1,
    paired = TRUE,
    method = "delong"
  )
  
  lrt_obj <- anova(mod0, mod1, test = "LRT")
  lrt_chisq <- lrt_obj$Deviance[2]
  lrt_df <- lrt_obj$Df[2]
  lrt_p <- lrt_obj$`Pr(>Chi)`[2]
  
  addval <- calc_idi_cnri(y, p0, p1)
  
  boot_obj <- bootstrap_compare_models(
    data = data,
    base_model = base_model,
    new_model = new_model,
    predictor_list = predictor_list,
    B = B,
    seed = seed
  )
  
  data.frame(
    base_model = base_model,
    new_model = new_model,
    
    auc_base = auc0,
    auc_base_lcl = as.numeric(ci0[1]),
    auc_base_ucl = as.numeric(ci0[3]),
    
    auc_new = auc1,
    auc_new_lcl = as.numeric(ci1[1]),
    auc_new_ucl = as.numeric(ci1[3]),
    
    delta_auc = auc1 - auc0,
    delta_auc_boot_lcl = boot_obj$delta_auc_ci[1],
    delta_auc_boot_ucl = boot_obj$delta_auc_ci[2],
    p_delong = delong_obj$p.value,
    p_delta_auc_boot = boot_obj$p_delta_auc,
    
    lrt_chisq = lrt_chisq,
    lrt_df = lrt_df,
    p_lrt = lrt_p,
    
    idi = addval$idi,
    idi_boot_lcl = boot_obj$idi_ci[1],
    idi_boot_ucl = boot_obj$idi_ci[2],
    p_idi_boot = boot_obj$p_idi,
    
    cnri = addval$cnri,
    cnri_boot_lcl = boot_obj$cnri_ci[1],
    cnri_boot_ucl = boot_obj$cnri_ci[2],
    p_cnri_boot = boot_obj$p_cnri,
    
    boot_success = boot_obj$n_success
  )
}

compare_res <- bind_rows(
  lapply(compare_pairs, function(x) {
    compare_two_models(
      base_model = x[1],
      new_model = x[2],
      logistic_mod_list = logistic_mod_list,
      pred_prob_list = pred_prob_list,
      roc_list = roc_list,
      data = prediction_data,
      predictor_list = predictor_list,
      B = 2000,
      seed = 123
    )
  })
)

write_xlsx(
  compare_res,
  "results/nhs_caco_dementia_NULISA_prediction.xlsx"
)

roc_df <- do.call(
  rbind,
  lapply(names(roc_list), function(m) {
    
    df <- data.frame(
      tpr = roc_list[[m]]$sensitivities,
      fpr = 1 - roc_list[[m]]$specificities,
      model = m
    )
    
    df[order(df$fpr, df$tpr), ]
  })
)

auc_labels <- sapply(names(roc_list), function(m) {
  
  auc_val <- as.numeric(pROC::auc(roc_list[[m]]))
  ci_val  <- as.numeric(pROC::ci.auc(roc_list[[m]]))
  
  paste0(
    round(auc_val, 2),
    " (",
    round(ci_val[1], 2),
    ", ",
    round(ci_val[3], 2),
    ")"
  )
})

auc_only <- sapply(names(roc_list), function(m) {
  round(as.numeric(pROC::auc(roc_list[[m]])), 2)
})

names(auc_labels) <- names(roc_list)
names(auc_only) <- names(roc_list)

roc_df$model_label <- auc_labels[roc_df$model]
roc_df$auc_only <- auc_only[roc_df$model]

write_xlsx(
  roc_df,
  "results/nhs_caco_dementia_NULISA_auc_for_fig.xlsx"
)

## SHAP

predict_function <- function(object, newdata) {
  predict(object, newdata = newdata, type = "response")
}

set.seed(123)

shap_x <- as.data.frame(prediction_data)[
  ,
  names(coef(logistic_mod_list[["m13"]]))[-1],
  drop = FALSE
]

shap_values_all <- fastshap::explain(
  logistic_mod_list[["m13"]],
  X = shap_x,
  pred_wrapper = predict_function,
  nsim = 500
)

shap_values <- as.data.frame(shap_values_all)

shap_values$atxn_shap <- rowSums(
  dplyr::select(
    shap_values,
    all_of(core_list[["all"]][["core3"]])
  ),
  na.rm = TRUE
)

shap_importance_all <- as.data.frame(colMeans(abs(shap_values))) %>%
  tibble::rownames_to_column(var = "predictor")

colnames(shap_importance_all) <- c("predictor", "overall")

write_xlsx(
  shap_values,
  "results/nhs_caco_dementia_NULISA_shap_for_fig.xlsx"
)

write_xlsx(
  shap_importance_all,
  "results/nhs_caco_dementia_NULISA_importance_for_fig.xlsx"
)

##-------------------------------
## NULISA and TICS
##-------------------------------

tics_res_list <- list()
i <- 1

for (prot in c(protall_list, scoreall_list)) {
  
  message(prot)
  
  for (subset in c("all", "control")) {
    
    if (subset == "all") {
      tmp_data <- model_data
    } else if (subset == "control") {
      tmp_data <- model_data %>% filter(caco == 0)
    }
    
    tmp_data_scaled <- tmp_data %>%
      mutate(!!prot := as.numeric(scale(.data[[prot]])))
    
    for (tics in tics_list) {
      
      tmp_mod <- tryCatch(
        glm(
          reformulate(c(prot, covar_list), response = tics),
          data = tmp_data
        ),
        error = function(e) e,
        warning = function(w) w
      )
      
      tmp_mod_scaled <- tryCatch(
        glm(
          reformulate(c(prot, covar_list), response = tics),
          data = tmp_data_scaled
        ),
        error = function(e) e,
        warning = function(w) w
      )
      
      tics_res_list[[i]] <- c(
        prot,
        tics,
        subset,
        "raw",
        if (inherits(tmp_mod, "glm")) {
          c(
            sum(!is.na(tmp_data[[tics]])),
            nobs(tmp_mod),
            summary(tmp_mod)$coefficients[prot, c(1, 2, 4)]
          )
        } else {
          c(tmp_mod$message, rep(NA, 4))
        }
      )
      
      i <- i + 1
      
      tics_res_list[[i]] <- c(
        prot,
        tics,
        subset,
        "scaled",
        if (inherits(tmp_mod_scaled, "glm")) {
          c(
            sum(!is.na(tmp_data_scaled[[tics]])),
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
}

tics_res_dat <- as.data.frame(do.call("rbind", tics_res_list))

colnames(tics_res_dat) <- c(
  "targetName",
  "outcome",
  "subset",
  "scale",
  "n",
  "n_model",
  "beta",
  "se",
  "p"
)

tics_res_dat <- tics_res_dat %>% 
  mutate(
    n = as.numeric(n),
    n_model = as.numeric(n_model),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) %>%
  left_join(marker_dict, by = "targetName")

write_xlsx(
  tics_res_dat,
  "results/nhs_caco_tics_NULISA.xlsx"
)