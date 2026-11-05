library(dplyr)
library(tidyr)
library(haven)
library(readxl)
library(writexl)

rm(list=ls())

##---------------------------------------------------------
## Load data
##---------------------------------------------------------

baselab <- read_sas("lab/duabasemindlabres.sas7bdat")
mrihippov <- read_sas("mri/duamrihippov.sas7bdat")
sex <- read_sas("comorbidities/duabasechar.sas7bdat") |> 
  select(mid, psv_gender)

##---------------------------------------------------------
## Process MRI
##---------------------------------------------------------

mri_long <- mrihippov %>%
  mutate(hippo_raw = if_else(qc == 0,
                              (lhippovolume + rhippovolume) / 2,
                              NA_real_)) %>%
  select(mid, visit, hippo_raw) %>%
  mutate(hippo_norm = hippo_raw / icv * 1000)

mri_wide <- mri_long %>%
  mutate(time = case_when(
    visit == "baseline" ~ "base",
    visit == "month36" ~ "m36",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(time)) %>%
  select(mid, assignment, time, mri_age, scanner, hippo_norm) %>%
  pivot_wider(
    names_from = time,
    values_from = c(mri_age, scanner, hippo_norm),
    names_sep = "_"
  ) %>%
  mutate(hippo_change = hippo_norm_m36 - hippo_norm_base)

##---------------------------------------------------------
## NHS-derived directions for AT(N) score
##---------------------------------------------------------

nulisa_dementia_res <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(scale == "scaled")

sign_vec <- nulisa_dementia_res %>%
  transmute(targetName, sign = ifelse(as.numeric(beta) >= 0, 1, -1)) %>%
  tibble::deframe()

atn_signs <- c(
  ab42_40 = sign_vec[["Ab42_div_Ab40_log2"]],
  ptau181 = sign_vec[["pTau_181"]],
  gfap    = sign_vec[["GFAP"]],
  nfl     = sign_vec[["NEFL"]]
)

##---------------------------------------------------------
## Process biomarkers
##---------------------------------------------------------

markers <- c("Abeta 40", "Abeta 42", "pTau-181", "GFAP", "NF-L")

bio0 <- baselab %>%
  filter(LABTEST %in% markers) %>%
  select(mid, LABTEST, RESULT) %>%
  pivot_wider(names_from = LABTEST, values_from = RESULT) %>%
  rename(
    ab40 = `Abeta 40`, ab42 = `Abeta 42`,
    ptau181 = `pTau-181`, gfap = GFAP, nfl = `NF-L`
  ) %>%
  mutate(
    ab42_40 = ab42 / ab40,
    z_ab42_40_raw = as.numeric(scale(ab42_40)),
    z_ptau181_log = as.numeric(scale(log(ptau181))),
    z_gfap_log = as.numeric(scale(log(gfap))),
    z_nfl_log = as.numeric(scale(log(nfl))),
    z_ab42_40_log = as.numeric(scale(log(ab42_40))),
    atn_ew_signed =
      z_ab42_40_log * atn_signs["ab42_40"] +
      z_ptau181_log * atn_signs["ptau181"] +
      z_gfap_log * atn_signs["gfap"] +
      z_nfl_log * atn_signs["nfl"],
    z_atn_ew_signed = as.numeric(scale(atn_ew_signed))
  )

##---------------------------------------------------------
## Model
##---------------------------------------------------------

model_data <- mri_wide %>%
  inner_join(bio0, by = "mid") %>%
  left_join(sex, by = "mid") %>%
  mutate(
    assignment = factor(assignment,
                        levels = c("Control Diet Group", "MIND Diet Group")),
    scanner_base = factor(scanner_base),
    sex = factor(psv_gender)
  )

marker_list <- c(
  "z_atn_ew_signed", "z_ab42_40_raw", "z_ptau181_log",
  "z_gfap_log", "z_nfl_log"
)

mind_res_list <- list()
i <- 1

for (marker in marker_list) {
  form_int <- as.formula(paste0(
    "hippo_change ~ ", marker,
    " * assignment + mri_age_base + sex + scanner_base"
  ))
  fit_int <- lm(form_int, data = model_data)
  coef_int <- summary(fit_int)$coefficients
  int_row <- grep(paste0("^", marker, ":assignment"), rownames(coef_int))

  group_levels <- levels(model_data$assignment)
  strat_res <- list()
  for (g in group_levels) {
    tmp <- model_data %>% filter(assignment == g)
    form_strat <- as.formula(paste0(
      "hippo_change ~ ", marker, " + mri_age_base + sex + scanner_base"
    ))
    fit_strat <- lm(form_strat, data = tmp)
    coef_strat <- summary(fit_strat)$coefficients
    strat_res[[g]] <- c(
      n = nobs(fit_strat),
      beta = coef_strat[marker, "Estimate"],
      se = coef_strat[marker, "Std. Error"],
      p = coef_strat[marker, "Pr(>|t|)"]
    )
  }

  mind_res_list[[i]] <- data.frame(
    marker = marker,
    outcome = "hippo_change",
    n_int = nobs(fit_int),
    beta_int = if (length(int_row) == 1) coef_int[int_row, "Estimate"] else NA_real_,
    se_int = if (length(int_row) == 1) coef_int[int_row, "Std. Error"] else NA_real_,
    p_int = if (length(int_row) == 1) coef_int[int_row, "Pr(>|t|)"] else NA_real_,
    group1 = group_levels[1],
    n_group1 = strat_res[[group_levels[1]]]["n"],
    beta_group1 = strat_res[[group_levels[1]]]["beta"],
    se_group1 = strat_res[[group_levels[1]]]["se"],
    p_group1 = strat_res[[group_levels[1]]]["p"],
    group2 = group_levels[2],
    n_group2 = strat_res[[group_levels[2]]]["n"],
    beta_group2 = strat_res[[group_levels[2]]]["beta"],
    se_group2 = strat_res[[group_levels[2]]]["se"],
    p_group2 = strat_res[[group_levels[2]]]["p"]
  )
  i <- i + 1
}

mind_res <- bind_rows(mind_res_list)
write_xlsx(mind_res, "results/mind_mri_res.xlsx")