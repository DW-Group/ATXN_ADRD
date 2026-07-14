library(dplyr)
library(writexl)

rm(list = ls())

##-------------------------------
## Load data
##-------------------------------

load("data/processed_data_dp.RData")

dir.create("results", recursive = TRUE, showWarnings = FALSE)

##-------------------------------
## Baseline NULISA
##-------------------------------

base_data <- merged_data %>% 
  filter(time == "A")

base_data_list <- list(
  "age_ge_50" = list(
    "all" = base_data %>% filter(age_baseline >= 50)
  )
)

##-------------------------------
## DIRECT-PLUS baseline diet analyses
##-------------------------------

diet_analysis_list <- c("MIND14", "amed")

prot_analysis_list <- unique(c(
  "all_core3_ew_signed",
  core_list[["all"]][["core3"]],
  core_list[["BD"]][["core3"]]
))

base_diet_res_list <- list()
i <- 1

for (prot in prot_analysis_list) {
  
  message(prot)
  
  for (age_group in names(base_data_list)) {
    
    for (gender_group in names(base_data_list[[age_group]])) {
      
      tmp_dat <- base_data_list[[age_group]][[gender_group]]
      
      tmp_covar <- c(covar_base, covar_diet)
      
      for (diet in diet_analysis_list) {
        
        tmp_mod <- tryCatch(
          glm(
            as.formula(
              paste(
                prot,
                " ~ ",
                paste(unique(c(diet, tmp_covar)), collapse = " + ")
              )
            ),
            data = tmp_dat
          ),
          error = function(e) e,
          warning = function(w) w
        )
        
        base_diet_res_list[[i]] <- c(
          prot,
          diet,
          age_group,
          gender_group,
          "univariate",
          if (inherits(tmp_mod, "glm"))
            c(
              nobs(tmp_mod),
              summary(tmp_mod)$coefficients[diet, c(1, 2, 4)]
            )
          else c(tmp_mod$message, rep(NA, 3))
        )
        
        i <- i + 1
      }
    }
  }
}

base_diet_res_dat <- as.data.frame(do.call("rbind", base_diet_res_list))

colnames(base_diet_res_dat) <- c(
  "targetName",
  "diet",
  "age_group",
  "gender_group",
  "model",
  "n",
  "beta",
  "se",
  "p"
)

base_diet_res_dat <- base_diet_res_dat %>% 
  mutate(
    n = suppressWarnings(as.numeric(n)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    p = suppressWarnings(as.numeric(p))
  ) %>% 
  left_join(marker_dict, by = "targetName")

write_xlsx(
  base_diet_res_dat,
  "results/dp_NULISA_baseline_diet.xlsx"
)