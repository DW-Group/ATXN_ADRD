library(dplyr)
library(tidyverse)
library(readxl)

rm(list = ls())

##-------------------------------
## Load results
##-------------------------------

nulisa_dementia_res <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(scale == "scaled")

##-------------------------------
## Select markers
##-------------------------------

ab_list <- c("Ab42_div_Ab40_log2")

tau_list <- c("pTau_217", "pTau_231", "pTau_181")
tauBD_list <- c("BD_pTau_217", "BD_pTau_231", "BD_pTau_181")

n_list <- c("NEFL", "GFAP")
x_list <- c("ICAM1", "GDF15")

core_list <- list(
  "all" = list(
    "core3" = c(tau_list, n_list, x_list, ab_list)
  ),
  "BD" = list(
    "core3" = c(tauBD_list, n_list, x_list, ab_list)
  )
)

save(core_list, file = "data/NULISA_core.RData")

##-------------------------------
## Derive weights
##-------------------------------

beta_vec <- nulisa_dementia_res %>%
  select(targetName, beta) %>%
  mutate(beta = as.numeric(beta)) %>%
  deframe()

sign_vec <- nulisa_dementia_res %>%
  transmute(
    targetName,
    sign = ifelse(as.numeric(beta) >= 0, 1, -1)
  ) %>%
  deframe()

##-------------
## NHS
##-------------

load("data/data_NULISA_nhs.RData")
nhs_npq_all <- npq_all

nhs_score_list <- list()
i <- 1

for (tau_set in names(core_list)) {
  for (set in names(core_list[[tau_set]])) {
    
    markers <- core_list[[tau_set]][[set]]
    betas   <- beta_vec[markers]
    signs   <- sign_vec[markers]
    
    tmp <- nhs_npq_all %>%
      select(id, all_of(markers)) %>%
      mutate(across(all_of(markers), ~ as.numeric(scale(.x)))) %>%
      mutate(
        !!paste0(tau_set, "_", set, "_ew") :=
          rowSums(across(all_of(markers)), na.rm = TRUE),
        !!paste0(tau_set, "_", set, "_ew_signed") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% signs),
        !!paste0(tau_set, "_", set, "_beta") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% betas)
      )
    
    nhs_score_list[[i]] <- tmp %>%
      select(id, starts_with(paste0(tau_set, "_", set)))
    
    i <- i + 1
  }
}

nhs_score_dat <- Reduce(function(x, y) full_join(x, y, by = "id"), nhs_score_list)

write.csv(
  nhs_score_dat,
  "data/data_NULISA_nhs_score.csv",
  row.names = FALSE
)

##-------------
## NHS 2nd
##-------------

load("data/data_NULISA_nhs_2nd_calibrated.RData")
nhs_npq_all <- npq_all

nhs_score_list <- list()
i <- 1

for (tau_set in names(core_list)) {
  for (set in names(core_list[[tau_set]])) {
    
    markers <- core_list[[tau_set]][[set]]
    betas   <- beta_vec[markers]
    signs   <- sign_vec[markers]
    
    tmp <- nhs_npq_all %>%
      select(id, all_of(markers)) %>%
      mutate(across(all_of(markers), ~ as.numeric(scale(.x)))) %>%
      mutate(
        !!paste0(tau_set, "_", set, "_ew") :=
          rowSums(across(all_of(markers)), na.rm = TRUE),
        !!paste0(tau_set, "_", set, "_ew_signed") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% signs),
        !!paste0(tau_set, "_", set, "_beta") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% betas)
      )
    
    nhs_score_list[[i]] <- tmp %>%
      select(id, starts_with(paste0(tau_set, "_", set)))
    
    i <- i + 1
  }
}

nhs_score_dat <- Reduce(function(x, y) full_join(x, y, by = "id"), nhs_score_list)

write.csv(
  nhs_score_dat,
  "data/data_NULISA_nhs_2nd_score_calibrated.csv",
  row.names = FALSE
)

##-------------
## NHS caco
##-------------

load("data/data_NULISA_nhs_caco.RData")
nhs_npq_all <- npq_all

nhs_score_list <- list()
i <- 1

for (tau_set in names(core_list)) {
  for (set in names(core_list[[tau_set]])) {
    
    markers <- core_list[[tau_set]][[set]]
    betas   <- beta_vec[markers]
    signs   <- sign_vec[markers]
    
    tmp <- nhs_npq_all %>%
      select(id, all_of(markers)) %>%
      mutate(across(all_of(markers), ~ as.numeric(scale(.x)))) %>%
      mutate(
        !!paste0(tau_set, "_", set, "_ew") :=
          rowSums(across(all_of(markers)), na.rm = TRUE),
        !!paste0(tau_set, "_", set, "_ew_signed") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% signs),
        !!paste0(tau_set, "_", set, "_beta") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% betas)
      )
    
    nhs_score_list[[i]] <- tmp %>%
      select(id, starts_with(paste0(tau_set, "_", set)))
    
    i <- i + 1
  }
}

nhs_score_dat <- Reduce(function(x, y) full_join(x, y, by = "id"), nhs_score_list)

write.csv(
  nhs_score_dat,
  "data/data_NULISA_nhs_caco_score.csv",
  row.names = FALSE
)

##-------------
## DIRECT-PLUS
##-------------

load("data/data_NULISA_dp.RData")
dp_npq_all <- npq_all

dp_score_list <- list()
i <- 1

for (tau_set in names(core_list)) {
  for (set in names(core_list[[tau_set]])) {
    
    markers <- core_list[[tau_set]][[set]]
    betas   <- beta_vec[markers]
    signs   <- sign_vec[markers]
    
    tmp <- dp_npq_all %>%
      select(tsno, all_of(markers)) %>%
      mutate(across(all_of(markers), ~ as.numeric(scale(.x)))) %>%
      mutate(
        !!paste0(tau_set, "_", set, "_ew") :=
          rowSums(across(all_of(markers)), na.rm = TRUE),
        !!paste0(tau_set, "_", set, "_ew_signed") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% signs),
        !!paste0(tau_set, "_", set, "_beta") :=
          as.numeric(as.matrix(across(all_of(markers))) %*% betas)
      )
    
    dp_score_list[[i]] <- tmp %>%
      select(tsno, starts_with(paste0(tau_set, "_", set)))
    
    i <- i + 1
  }
}

dp_score_dat <- Reduce(function(x, y) full_join(x, y, by = "tsno"), dp_score_list)

write.csv(
  dp_score_dat,
  "data/data_NULISA_dp_score.csv",
  row.names = FALSE
)