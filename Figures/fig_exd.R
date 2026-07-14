library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(readxl)
library(cowplot)
library(patchwork)
library(pheatmap)
library(vegan)
library(writexl)

##---------------------------------------------------------
## Exd Data Fig 1
##---------------------------------------------------------

##------------------
## 1a
##------------------

rm(list=ls())

source("scripts/source_rm_markers.R")

outcome_dict <- c(
  "av_nverbl" = "Verbal memory",
  "av_zscre"  = "TICS score"
)

source_data_ext1a_list <- list()

plot_tics_outcome <- function(outcome_use, title_use) {
  
  fig_dat <- read_excel("results/nhs_tics_NULISA.xlsx") %>%
    filter(
      scale == "scaled",
      LOD_flag != TRUE | is.na(LOD_flag),
      !(targetName %in% ipc_bcg_fails),
      outcome == outcome_use
    ) %>%
    mutate(
      targetName_anno = case_when(
        !is.na(targetName_orig) ~ targetName_orig,
        targetName == "Ab42_div_Ab40_log2" ~ "Aβ42/Aβ40",
        targetName == "AT_comp" ~ "Amyloid-Tau score",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(targetName_anno)) %>%
    mutate(
      p         = as.numeric(p),
      beta      = as.numeric(beta),
      fdr       = p.adjust(p, method = "BH"),
      neglog10p = -log10(p),
      sig_cat   = case_when(
        fdr < 0.05 ~ "FDR < 0.05",
        p < 0.05   ~ "p < 0.05",
        TRUE       ~ "NS"
      )
    )
  
  lab_dat <- fig_dat %>% filter(p < 0.05)
  
  source_data_ext1a_list[[outcome_use]] <<- fig_dat %>%
    mutate(
      outcome_label = title_use,
      label_shown = ifelse(p < 0.05, targetName_anno, "")
    )
  
  yline <- -log10(0.05)
  
  ggplot(fig_dat, aes(x = beta, y = neglog10p)) +
    geom_hline(
      yintercept = yline,
      linetype = 2,
      color = "gray31",
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = max(fig_dat$beta, na.rm = TRUE),
      y = yline,
      label = "italic(P) == 0.05",
      parse = TRUE,
      color = "gray31",
      vjust = -0.5,
      hjust = 0.8,
      size = 3.2
    ) +
    geom_point(
      aes(
        shape = sig_cat,
        alpha = sig_cat,
        fill  = sig_cat,
        color = sig_cat
      ),
      size = 2.5,
      stroke = 0.4
    ) +
    scale_shape_manual(values = c(
      "FDR < 0.05" = 23,
      "p < 0.05"   = 22,
      "NS"         = 16
    )) +
    scale_color_manual(values = c(
      "FDR < 0.05" = "black",
      "p < 0.05"   = "black",
      "NS"         = "#BDBDBD"
    )) +
    scale_alpha_manual(
      values = c(
        "FDR < 0.05" = 0.95,
        "p < 0.05"   = 0.95,
        "NS"         = 0.75
      ),
      guide = "none"
    ) +
    scale_fill_manual(values = c(
      "FDR < 0.05" = "#6E597B",
      "p < 0.05"   = "#b0cba9",
      "NS"         = "#BDBDBD"
    )) +
    geom_text_repel(
      data = lab_dat,
      aes(label = targetName_anno),
      segment.size = 0.2,
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      seed = 1
    ) +
    labs(
      x = expression(beta ~ "coefficient"),
      y = expression(-log[10](italic(P))),
      title = title_use
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.2),
      axis.ticks.y = element_line(color = "black", linewidth = 0.2),
      axis.title = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 12)
    )
}

p_verbal <- plot_tics_outcome("av_nverbl", "Cognitive function (verbal memory)")
p_tics   <- plot_tics_outcome("av_zscre",  "Cognitive function (TICS score)")

p_ext1a <- p_verbal + p_tics + plot_layout(ncol = 2)

ggsave(
  "results/figures/exd_fig1a_secondary_cognitive_outcomes.pdf",
  p_ext1a,
  width = 7,
  height = 3.5,
  dpi = 600,
  bg = "white"
)

### Source data

source_data_ext1a_points <- bind_rows(source_data_ext1a_list) %>%
  transmute(
    outcome,
    protein_label = targetName_anno,
    beta,
    p_value = p
  ) %>%
  arrange(outcome, protein_label)

write_xlsx(ext1a_points = source_data_ext1a_points, "results/source_data/source_data_extended_data_fig1a.xlsx")

##------------------
## 1b-d
##------------------

rm(list = ls())

load("data/processed_data_nhs.RData")
source("scripts/source_rm_markers.R")

table(model_data$hiedu2_3, useNA = "ifany")
table(model_data$ckd_flag, useNA = "ifany")
summary(model_data$bmi88_imp)

top_dat <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(
    scale == "scaled",
    LOD_flag != TRUE | is.na(LOD_flag),
    !(targetName %in% ipc_bcg_fails),
    p < 0.05
  ) %>%
  mutate(
    p = as.numeric(p),
    beta_primary = as.numeric(beta),
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "AT_comp" ~ "Amyloid-Tau score",
      TRUE ~ targetName
    )
  ) %>%
  arrange(beta_primary)

top_markers <- top_dat$targetName
marker_order <- top_dat$targetName_anno

top_dat %>%
  select(targetName, targetName_anno, beta, p)

### Function to run sensitivity analysis

run_sensitivity <- function(data_use, covars_use, analysis_label) {
  
  res_list <- list()
  
  for (prot in top_markers) {
    
    message(paste0(analysis_label, ": ", prot))
    
    tmp_dat <- data_use %>%
      mutate(prot_z = as.numeric(scale(.data[[prot]])))
    
    tmp_formula <- as.formula(
      paste(
        "Surv(tad_final, adcase) ~ prot_z +",
        paste(covars_use, collapse = " + ")
      )
    )
    
    tmp_mod <- tryCatch(
      coxph(
        tmp_formula,
        ties = "breslow",
        data = tmp_dat
      ),
      error = function(e) e
    )
    
    if (inherits(tmp_mod, "coxph")) {
      
      tmp_coef <- summary(tmp_mod)$coefficients["prot_z", ]
      
      res_list[[prot]] <- tibble(
        targetName = prot,
        analysis = analysis_label,
        n = tmp_mod$n,
        n_events = tmp_mod$nevent,
        beta = as.numeric(tmp_coef["coef"]),
        se = as.numeric(tmp_coef["se(coef)"]),
        p = as.numeric(tmp_coef["Pr(>|z|)"])
      )
      
    } else {
      
      res_list[[prot]] <- tibble(
        targetName = prot,
        analysis = analysis_label,
        n = NA_integer_,
        n_events = NA_integer_,
        beta = NA_real_,
        se = NA_real_,
        p = NA_real_
      )
      
    }
  }
  
  bind_rows(res_list) %>%
    left_join(
      top_dat %>%
        select(
          targetName,
          targetName_anno,
          beta_primary
        ),
      by = "targetName"
    )
}

### Run sensitivity analyses

## Additional adjustment for BMI
res_bmi <- run_sensitivity(
  data_use = model_data,
  covars_use = c(covar_list, "bmi88_imp"),
  analysis_label = "Additional adjustment for BMI"
)

## Additional adjustment for education
res_edu <- run_sensitivity(
  data_use = model_data,
  covars_use = c(covar_list, "hiedu2_3"),
  analysis_label = "Additional adjustment for education"
)

## Exclude the participant with CKD
res_ckd <- run_sensitivity(
  data_use = model_data %>%
    filter(is.na(ckd_flag) | !ckd_flag),
  covars_use = covar_list,
  analysis_label = "Excluding participant with CKD"
)

## Combine

sens_res <- bind_rows(
  res_bmi,
  res_edu,
  res_ckd
) %>%
  mutate(
    targetName_anno = factor(
      targetName_anno,
      levels = marker_order
    )
  )

sens_res

common_lim <- range(
  c(sens_res$beta_primary, sens_res$beta),
  na.rm = TRUE
)

common_pad <- 0.05 * diff(common_lim)

common_lim <- c(
  common_lim[1] - common_pad,
  common_lim[2] + common_pad
)

## Scatter plots

plot_sensitivity <- function(dat, title_use, panel_tag) {
  
  r_val <- cor(
    dat$beta_primary,
    dat$beta,
    use = "complete.obs",
    method = "pearson"
  )
  
  ggplot(
    dat,
    aes(
      x = beta_primary,
      y = beta
    )
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = "gray40",
      linewidth = 0.4
    ) +
    geom_point(
      shape = 21,
      size = 2.5,
      stroke = 0.4,
      fill = "#B24C4C",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      aes(label = targetName_anno),
      segment.size = 0.2,
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      seed = 1
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = paste0("italic(r) == ", sprintf("%.2f", r_val)),
      parse = TRUE,
      hjust = -0.1,
      vjust = 1.3,
      size = 3.5
    ) +
    labs(
      x = expression(beta ~ "coefficient in primary analysis"),
      y = expression(beta ~ "coefficient in sensitivity analysis"),
      title = title_use,
      tag = panel_tag
    ) +
    coord_equal(
      xlim = common_lim,
      ylim = common_lim
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 11, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      ),
      plot.title = element_text(
        size = 11,
        hjust = 0.5
      ),
      plot.tag = element_text(
        size = 12,
        face = "bold"
      )
    )
}

p_bmi <- plot_sensitivity(
  dat = res_bmi,
  title_use = "Additional adjustment for BMI",
  panel_tag = "b"
)

p_edu <- plot_sensitivity(
  dat = res_edu,
  title_use = "Additional adjustment for education",
  panel_tag = "c"
)

p_ckd <- plot_sensitivity(
  dat = res_ckd,
  title_use = "Excluding participant with CKD",
  panel_tag = "d"
)

p_ext1_bcd <- p_bmi + p_edu + p_ckd +
  plot_layout(ncol = 3,
              widths = c(1, 1, 1))

ggsave(
  filename = "results/figures/exd_fig1b_d_dementia_sensitivity.pdf",
  plot = p_ext1_bcd,
  width = 10.5,
  height = 3.5,
  dpi = 600,
  bg = "white"
)

### Source data

source_data_ext1b_d_points <- sens_res %>%
  mutate(
    neg_log10_p = -log10(p),
    label_shown = as.character(targetName_anno)
  ) %>%
  transmute(
    protein_label = as.character(targetName_anno),
    beta_primary,
    beta_sensitivity = beta
  )

write_xlsx(
  list(
    ext1bd_points = source_data_ext1b_d_points
  ),
  "results/source_data/source_data_extended_data_fig1b_d.xlsx"
)

##---------------------------------------------------------
## Exd Data Fig 2
##---------------------------------------------------------

##------------------
## 2a
##------------------

rm(list = ls())

load("data/processed_data_dp.RData")
load("data/NULISA_core.RData")
source("scripts/source_col.R")

marker_dict %>% filter(LOD_flag)

protall_list <- prot_list[
  !grepl("raw", prot_list) &
    !(prot_list %in% (marker_dict %>% filter(LOD_flag) %>% pull(targetName)))
]

##-------------------------------
## Keep baseline + Month 18 complete cases only
##-------------------------------

fig_dat <- npq_all %>%
  dplyr::select(sno, time, group, all_of(protall_list)) %>%
  filter(time %in% c("A", "C")) %>%   ## A = baseline, C = Month 18
  pivot_wider(
    names_from  = time,
    values_from = all_of(protall_list),
    names_sep   = "_"
  ) %>%
  filter(
    if_all(ends_with("_A"), ~ !is.na(.)) &
      if_all(ends_with("_C"), ~ !is.na(.))
  ) %>%
  mutate(
    group_anno = factor(group, levels = c("HDG", "MED", "GreenMED"))
  )

##-------------------------------
## PCA using baseline + Month 18 only
##-------------------------------

X_for_pca <- bind_rows(
  fig_dat %>% dplyr::select(ends_with("_A")) %>% setNames(protall_list),
  fig_dat %>% dplyr::select(ends_with("_C")) %>% setNames(protall_list)
)

stopifnot(!anyNA(X_for_pca))

pca_fit <- prcomp(X_for_pca, center = TRUE, scale. = TRUE)

pc_A <- as.data.frame(
  predict(
    pca_fit,
    newdata = fig_dat %>% dplyr::select(ends_with("_A")) %>% setNames(protall_list)
  )
) %>%
  mutate(
    sno       = fig_dat$sno,
    group     = fig_dat$group_anno,
    time_anno = "Baseline"
  )

pc_C <- as.data.frame(
  predict(
    pca_fit,
    newdata = fig_dat %>% dplyr::select(ends_with("_C")) %>% setNames(protall_list)
  )
) %>%
  mutate(
    sno       = fig_dat$sno,
    group     = fig_dat$group_anno,
    time_anno = "Month 18"
  )

pc_long <- bind_rows(pc_A, pc_C) %>%
  mutate(
    time_anno = factor(time_anno, levels = c("Baseline", "Month 18")),
    group = factor(group, levels = c("HDG", "MED", "GreenMED"))
  )

pc_var <- data.frame(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  var_explained = (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
)

##-------------------------------
## Group-level PC trajectory summary
##-------------------------------

pc_summary <- pc_long %>%
  group_by(group, time_anno) %>%
  dplyr::summarise(
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC2_mean = mean(PC2, na.rm = TRUE),
    PC1_sd   = sd(PC1, na.rm = TRUE),
    PC2_sd   = sd(PC2, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    PC1_se = PC1_sd / sqrt(n),
    PC2_se = PC2_sd / sqrt(n),
    PC1_l  = PC1_mean - 1.96 * PC1_se,
    PC1_u  = PC1_mean + 1.96 * PC1_se,
    PC2_l  = PC2_mean - 1.96 * PC2_se,
    PC2_u  = PC2_mean + 1.96 * PC2_se
  )

##-------------------------------
## PERMANOVA on baseline-to-Month 18 change by intervention group
##-------------------------------

X_A <- fig_dat %>%
  dplyr::select(ends_with("_A")) %>%
  setNames(protall_list) %>%
  as.matrix()

X_C <- fig_dat %>%
  dplyr::select(ends_with("_C")) %>%
  setNames(protall_list) %>%
  as.matrix()

delta_mat <- X_C - X_A
delta_mat_scaled <- scale(delta_mat)

dist_delta <- dist(delta_mat_scaled, method = "euclidean")

adonis_res <- adonis2(
  dist_delta ~ group_anno,
  data = fig_dat %>% dplyr::select(group_anno),
  permutations = 999
)

adonis_tab <- as.data.frame(adonis_res)
adonis_tab$term <- rownames(adonis_tab)

group_row <- adonis_tab %>%
  filter(term == "Model")

permanova_label <- paste0(
  "PERMANOVA R² = ",
  sprintf("%.1f", 100 * group_row$R2),
  "%, P = ",
  format.pval(group_row$`Pr(>F)`, digits = 2, eps = 0.001)
)

print(adonis_res)

##-------------------------------
## Plot
##-------------------------------

p_pc_traj_mean_AC <- ggplot(
  pc_summary,
  aes(x = PC1_mean, y = PC2_mean, color = group)
) +
  
  geom_path(
    aes(group = group),
    linewidth = 0.7,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  
  geom_errorbar(
    aes(ymin = PC2_l, ymax = PC2_u),
    width = 0,
    linewidth = 0.5,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  
  geom_errorbarh(
    aes(xmin = PC1_l, xmax = PC1_u),
    height = 0,
    linewidth = 0.5,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  
  geom_point(
    aes(fill = group, shape = time_anno),
    size = 3.6,
    stroke = 0.4,
    color = "black",
    alpha = 0.9
  ) +
  
  scale_fill_manual(values = group_col, name = "Intervention") +
  scale_color_manual(values = group_col, guide = "none") +
  scale_shape_manual(
    values = c(Baseline = 21, `Month 18` = 22),
    name = "Time"
  ) +
  
  xlab(paste0("PC1 (", round(pc_var$var_explained[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(pc_var$var_explained[2] * 100, 1), "%)")) +
  ggtitle("Global proteomic changes in DIRECT-PLUS") +
  
  annotate(
    "text",
    x = max(pc_summary$PC1_u, na.rm = TRUE),
    y = max(pc_summary$PC2_u, na.rm = TRUE) -
      0.03 * diff(range(pc_summary$PC2_u, pc_summary$PC2_l, na.rm = TRUE)),
    label = permanova_label,
    hjust = 1,
    size = 4.2,
    color = "black"
  ) +
  
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 4, alpha = 1)
    ),
    shape = guide_legend(
      override.aes = list(fill = "grey70")
    )
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 14, colour = "black"),
    axis.text.y        = element_text(size = 14, colour = "black"),
    axis.title.x       = element_text(size = 14, colour = "black"),
    axis.title.y       = element_text(size = 14, colour = "black"),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title       = element_text(size = 12),
    legend.text        = element_text(size = 11),
    plot.title         = element_text(size = 15, color = "black", hjust = 0.5),
    legend.position    = c(0.16, 0.3),
    legend.background  = element_blank(),
    legend.box         = "vertical"
  )

ggsave(
  "results/figures/exd_fig2a_pca_dp_baseline_month18_by_group.pdf",
  p_pc_traj_mean_AC,
  width = 5,
  height = 5
)

### Source data

source_data_ext2a_trajectory <- pc_summary %>%
  transmute(
    intervention_group = as.character(group),
    time = as.character(time_anno),
    PC1_mean,
    PC2_mean,
    PC1_lower_95CI = PC1_l,
    PC1_upper_95CI = PC1_u,
    PC2_lower_95CI = PC2_l,
    PC2_upper_95CI = PC2_u
  ) %>%
  arrange(intervention_group, time)

write_xlsx(
  list(
    ext2a_trajectory = source_data_ext2a_trajectory
  ),
  "results/source_data/source_data_extended_data_fig2a.xlsx"
)

##------------------
## 2b
##------------------

rm(list = ls())

load("data/processed_data_dp.RData")
load("data/NULISA_core.RData")

marker_dict %>% filter(LOD_flag)

protall_list <- prot_list[
  !grepl("raw", prot_list) &
    !(prot_list %in% (marker_dict %>% filter(LOD_flag | ipc_bcg_flag) %>% pull(targetName)))
]

fig_dat <- npq_all %>%
  filter(group == "HDG") %>%   ## control only
  dplyr::select(sno, time, all_of(protall_list)) %>%
  pivot_wider(
    names_from  = time,
    values_from = all_of(protall_list),
    names_sep   = "_"
  ) %>%
  filter(
    if_all(ends_with("_A"), ~ !is.na(.)) &
      if_all(ends_with("_D"), ~ !is.na(.))
  )

X_A <- fig_dat %>%
  dplyr::select(ends_with("_A")) %>%
  setNames(protall_list)

X_D <- fig_dat %>%
  dplyr::select(ends_with("_D")) %>%
  setNames(protall_list)

X_pca <- bind_rows(X_A, X_D)
stopifnot(!anyNA(X_pca))

pca_fit <- prcomp(X_pca, center = TRUE, scale. = TRUE)

pc_A <- as.data.frame(predict(pca_fit, newdata = X_A)) %>%
  mutate(
    sno = fig_dat$sno,
    time_anno = "Baseline"
  )

pc_D <- as.data.frame(predict(pca_fit, newdata = X_D)) %>%
  mutate(
    sno = fig_dat$sno,
    time_anno = "FIT"
  )

pc_long <- bind_rows(pc_A, pc_D) %>%
  mutate(
    time_anno = factor(time_anno, levels = c("Baseline", "FIT"))
  )

pc_var <- data.frame(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  var_explained = (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
)

### PERMANOVA for overall time effect, control only

dist_mat <- dist(X_pca, method = "euclidean")

meta <- data.frame(
  sno  = rep(fig_dat$sno, 2),
  time = rep(c("Baseline", "FIT"), each = nrow(fig_dat))
)

adonis_res <- adonis2(
  dist_mat ~ time,
  data = meta,
  permutations = 9999,
  strata = meta$sno
)

permanova_label <- paste0(
  "PERMANOVA ",
  "R² = ", sprintf("%.1f", adonis_res$R2[1] * 100), "%\n",
  "P = ", sprintf("%.5f", adonis_res$`Pr(>F)`[1])
)

### Plot

p_pc_dp_AD_ctrl <- ggplot(pc_long, aes(x = PC1, y = PC2)) +
  
  geom_line(
    aes(group = sno),
    color = "gray31",
    alpha = 0.35,
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  
  stat_ellipse(
    aes(fill = time_anno),
    geom = "polygon",
    type = "norm",
    alpha = 0.10,
    color = NA,
    show.legend = FALSE
  ) +
  
  geom_point(
    aes(fill = time_anno),
    shape = 21,
    size = 3.4,
    stroke = 0.4,
    color = "gray31",
    alpha = 0.85
  ) +
  
  scale_fill_manual(
    values = c(
      "Baseline" = "#8DB7A8",  
      "FIT"      = "#D7A07A"  
    ),
    name = ""
  ) +
  
  xlab(paste0("Proteomic PC1 (", round(pc_var$var_explained[1] * 100, 1), "%)")) +
  ylab(paste0("Proteomic PC2 (", round(pc_var$var_explained[2] * 100, 1), "%)")) +
  ggtitle("Global proteomic changes in DIRECT-PLUS control arm") +
  
  annotate(
    "text",
    x = max(pc_long$PC1, na.rm = TRUE),
    y = max(pc_long$PC2, na.rm = TRUE) -
      0.05 * diff(range(pc_long$PC2, na.rm = TRUE)),
    label = permanova_label,
    hjust = 1,
    size = 4.2,
    color = "black"
  ) +
  
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 4, alpha = 1)
    )
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 14, colour = "black"),
    axis.text.y        = element_text(size = 14, colour = "black"),
    axis.title.x       = element_text(size = 14, colour = "black"),
    axis.title.y       = element_text(size = 14, colour = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    legend.title       = element_text(size = 11),
    legend.text        = element_text(size = 12),
    plot.title         = element_text(size = 15, color = "black", hjust = 0.5),
    legend.position    = c(0.85, 0.12),
    legend.background  = element_blank(),
    legend.box         = "vertical"
  )

ggsave(
  "results/figures/exd_fig2b_pca_dp_baseline_fit_ctrl_only.pdf",
  p_pc_dp_AD_ctrl,
  width = 5,
  height = 5
)

### Source data

source_data_ext2b_points <- pc_long %>%
  transmute(
    pair_id = as.integer(factor(sno)),
    time = as.character(time_anno),
    PC1,
    PC2
  ) %>%
  arrange(pair_id, time)

write_xlsx(
  list(
    ext2b_points = source_data_ext2b_points
  ),
  "results/source_data/source_data_extended_data_fig2b.xlsx"
)

##---------------------------------------------------------
## Exd Data Fig 3
##---------------------------------------------------------

##------------------
## 3a 
##------------------

rm(list = ls())

source("scripts/source_col.R")

### Load prospective

load("data/processed_data_nhs.RData")

pros_dat <- model_data %>%
  mutate(
    cohort = "NHS prospective",
    all_core3_ew_signed_z = as.numeric(scale(all_core3_ew_signed)),
    all_core3_ew_signed_resid = residuals(
      lm(
        all_core3_ew_signed_z ~ age_final + calor90n_cumavg,
        data = model_data,
        na.action = na.exclude
      )
    )
  ) %>%
  select(
    cohort,
    amed_90_cumavg,
    all_core3_ew_signed_resid,
    age_final,
    calor90n_cumavg
  )

### Load case-control, controls only

load("data/processed_data_nhs_caco.RData")

caco_control_dat <- model_data %>%
  filter(adcase == 0) %>%
  mutate(
    cohort = "NHS case-control controls",
    all_core3_ew_signed_z = as.numeric(scale(all_core3_ew_signed)),
    all_core3_ew_signed_resid = residuals(
      lm(
        all_core3_ew_signed_z ~ age_final + calor90n_cumavg,
        data = .,
        na.action = na.exclude
      )
    )
  ) %>%
  select(
    cohort,
    amed_90_cumavg,
    all_core3_ew_signed_resid,
    age_final,
    calor90n_cumavg
  )

### Pool after cohort-specific residualization

plot_dat <- bind_rows(
  pros_dat,
  caco_control_dat
)

plot_dat$cohort <- factor(
  plot_dat$cohort,
  levels = c("NHS prospective", "NHS case-control controls")
)

### Pooled descriptive linear fit

lm_fit_pool <- lm(
  all_core3_ew_signed_resid ~ amed_90_cumavg,
  data = plot_dat
)

sm_pool <- summary(lm_fit_pool)$coefficients

beta_pool <- sm_pool["amed_90_cumavg", "Estimate"]
p_pool <- sm_pool["amed_90_cumavg", "Pr(>|t|)"]

pred_df <- data.frame(
  amed_90_cumavg = seq(
    min(plot_dat$amed_90_cumavg, na.rm = TRUE),
    max(plot_dat$amed_90_cumavg, na.rm = TRUE),
    length.out = 200
  )
)

pred_mat <- predict(lm_fit_pool, newdata = pred_df, interval = "confidence")
pred_df$fit <- pred_mat[, "fit"]
pred_df$lwr <- pred_mat[, "lwr"]
pred_df$upr <- pred_mat[, "upr"]

### Plot

p_scatter_pool <- ggplot(
  plot_dat,
  aes(x = amed_90_cumavg, y = all_core3_ew_signed_resid)
) +
  geom_point(
    shape = 21,
    fill  = "grey66",
    color = "black",
    alpha = 0.3,
    size  = 2,
    stroke = 0.15
  ) +
  geom_ribbon(
    data = pred_df,
    aes(x = amed_90_cumavg, ymin = lwr, ymax = upr),
    fill = "#6E8FA8",
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = amed_90_cumavg, y = fit),
    color = "#6E8FA8",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = amed_90_cumavg, y = lwr),
    color = "#6E8FA8",
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = amed_90_cumavg, y = upr),
    color = "#6E8FA8",
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_rug(
    sides = "bl",
    alpha = 0.10,
    linewidth = 0.20
  ) +
  labs(
    x = "aMED score",
    y = "ATX(N) score (residual)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black"),
    axis.line  = element_line(linewidth = 0.25, color = "black"),
    plot.margin = margin(6, 8, 6, 6),
    panel.grid = element_blank()
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1.05, vjust = 1.15,
    size = 3,
    color = "grey20"
  )

ggsave(
  "results/figures/exd_fig3a_amed_scatter_pooled.pdf",
  p_scatter_pool,
  width = 2.6,
  height = 2.8
)

### Source data

source_data_ext3a_points <- plot_dat %>%
  transmute(
    cohort,
    aMED_score = amed_90_cumavg,
    ATXN_score_residual = all_core3_ew_signed_resid
  ) %>%
  arrange(cohort, aMED_score)

source_data_ext3a_fit <- pred_df %>%
  transmute(
    aMED_score = amed_90_cumavg,
    fitted_ATXN_score_residual = fit,
    lower_95CI = lwr,
    upper_95CI = upr
  )

write_xlsx(
  list(
    ext3a_points = source_data_ext3a_points,
    ext3a_fit = source_data_ext3a_fit
  ),
  "results/source_data/source_data_extended_data_fig3a.xlsx"
)

##------------------
## 3b 
##------------------

rm(list=ls())

### Data

load("data/NULISA_core.RData")
source("scripts/source_dietcomp.R")

marker_order <- c(
  "ATX(N) score",
  "pTau-217", "pTau-181", "pTau-231",
  "BD-pTau-217", "BD-pTau-181", "BD-pTau-231",
  "Abeta42/Abeta40",
  "ICAM1", "GDF15", "NEFL", "GFAP"
)

marker_order_bot <- marker_order[marker_order != "ATX(N) score"]

fig_dat <- read_excel("results/nhs_meta_NULISA_diet.xlsx") %>%
  left_join(nhs_dietcomp_dict, by="diet") %>%
  filter(
    targetName %in% c("all_core3_ew_signed", core_list[["all"]][["core3"]], core_list[["BD"]][["core3"]]),
    LOD_flag != TRUE | is.na(LOD_flag),
    !(grepl("comp", diet) & grepl("mutual", model)),
    grepl("aMED", score)
  ) %>%
  mutate(
    beta = as.numeric(beta_meta),
    se   = as.numeric(se_meta),
    p    = as.numeric(p_meta),
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      grepl("core", targetName) ~ "ATX(N) score",
      TRUE ~ NA_character_
    ),
    targetName_anno = factor(
      targetName_anno,
      levels = marker_order
    )
  ) %>% 
  group_by(diet) %>% 
  mutate(
    fdr = p.adjust(p, method = "BH"),
    sig = case_when(
      p < 0.01 ~ "#",
      p < 0.05 ~ "+",
      TRUE ~ ""
    )
  ) %>%
  arrange(targetName_anno)

wins_top_p <- 0
wins_bot_p <- 0.1

beta_top_q <- quantile(
  fig_dat %>% filter(targetName_anno == "ATX(N) score") %>% pull(beta),
  probs = c(wins_top_p, 1 - wins_top_p),
  na.rm = TRUE
)

beta_top_low  <- as.numeric(beta_top_q[1])
beta_top_high <- as.numeric(beta_top_q[2])

beta_bot_q <- quantile(
  fig_dat %>% filter(targetName_anno != "ATX(N) score") %>% pull(beta),
  probs = c(wins_bot_p, 1 - wins_bot_p),
  na.rm = TRUE
)

beta_bot_low  <- as.numeric(beta_bot_q[1])
beta_bot_high <- as.numeric(beta_bot_q[2])

### aMED heatmap

amed_dat <- fig_dat %>% 
  mutate(
    diet_anno = factor(
      diet_anno,
      levels = unique(nhs_dietcomp_dict %>% filter(grepl("aMED", score)) %>% pull(diet_anno))
    )
  )

nhs_amed_dat <- amed_dat %>% 
  ungroup() %>% 
  filter(diet == "amed_90_cumavg") %>% 
  dplyr::select(targetName_anno, beta_nhs = beta, se_nhs = se, p_nhs = p, sig_nhs = sig)

amed_top <- amed_dat %>% filter(targetName_anno == "ATX(N) score")
amed_bot <- amed_dat %>% filter(targetName_anno != "ATX(N) score")

p_amed_top <- ggplot(amed_top, aes(x = diet_anno, y = targetName_anno)) +
  geom_tile(aes(fill = beta), color = "gray31", linewidth = 0.2) +
  geom_text(aes(label = sig), color = "black", size = 3) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low  = "#7B9FC4",
    mid  = "white",
    high = "#C97A6D",
    midpoint = 0, limits = c(beta_top_low, beta_top_high)
  ) +
  labs(fill = expression(beta["ATX(N) score"])) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_text(color = "black", size = 7),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    panel.border = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  coord_fixed(ratio = 1)

p_amed_bot <- ggplot(amed_bot, aes(x = diet_anno, y = targetName_anno)) +
  geom_tile(aes(fill = beta), color = "gray42", linewidth = 0.2) +
  geom_text(aes(label = sig), color = "black", size = 3) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = rev(marker_order_bot), expand = c(0, 0)) +
  scale_fill_gradient2(
    low  = "#8EC0B8",
    mid  = "white",
    high = "#D8A35D",
    midpoint = 0, limits = c(beta_bot_low, beta_bot_high),
    oob = scales::squish
  ) +
  labs(fill = expression(beta["components"])) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x.bottom = element_text(color = "black", angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(color = "black", size = 7),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    panel.border = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  coord_fixed(ratio = 1)

### Save

final_plot <- p_amed_top / plot_spacer() / p_amed_bot +
  plot_layout(heights = c(1, 0.05, length(unique(amed_bot$targetName_anno))))

ggsave("results/figures/exd_fig3b_meta.pdf", final_plot, height = 4.5)

p_atxn_bar <- p_amed_top +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0)
  )
atxn_bar <- cowplot::get_legend(p_atxn_bar)
ggsave(
  "results/figures/exd_fig3b_meta_colorbar_ATXN.pdf",
  cowplot::ggdraw(atxn_bar),
  width = 2,
  height = 2
)

p_comp_bar <- p_amed_bot +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0)
  )
comp_bar <- cowplot::get_legend(p_comp_bar)
ggsave(
  "results/figures/exd_fig3b_meta_colorbar_components.pdf",
  cowplot::ggdraw(comp_bar),
  width = 2,
  height = 2
)

sig_leg <- data.frame(
  sig   = factor(c("#", "+"), levels = c("#", "+")),
  label = c("P < 0.01", "P < 0.05")
)
p_sig_leg <- ggplot(sig_leg, aes(x = 1, y = sig)) +
  geom_text(aes(label = sig), size = 5.5, vjust = 0.5) +
  geom_text(
    aes(x = 1.3, label = label),
    hjust = 0,
    size  = 4
  ) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(0.9, 2.4), clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(2, 4, 2, 2)
  )
ggsave(
  "results/figures/exd_fig3b_meta_sig_legend.pdf",
  p_sig_leg,
  width  = 1.5,
  height = 0.5
)

### Source data

source_data_ext3b_heatmap <- amed_dat %>%
  ungroup() %>%
  transmute(
    diet_label = as.character(diet_anno),
    protein_label = as.character(targetName_anno),
    beta,
    p_value = p
  ) %>%
  arrange(diet_label, protein_label)

write_xlsx(
  list(
    ext3b_heatmap = source_data_ext3b_heatmap
  ),
  "results/source_data/source_data_extended_data_fig3b.xlsx"
)

##------------------
## 3c 
##------------------

rm(list = ls())

### Data

load("data/NULISA_core.RData")
source("scripts/source_dietcomp.R")

fig_dat <- read_excel("results/nhs_meta_NULISA_diet.xlsx") %>%
  left_join(nhs_dietcomp_dict, by = "diet") %>%
  filter(
    targetName %in% c(
      "all_core3_ew_signed",
      core_list[["all"]][["core3"]],
      core_list[["BD"]][["core3"]]
    ),
    LOD_flag != TRUE | is.na(LOD_flag),
    !(grepl("comp", diet) & grepl("mutual", model)),
    grepl("aMED", score)
  ) %>%
  mutate(
    beta = as.numeric(beta_meta),
    se   = as.numeric(se_meta),
    p    = as.numeric(p_meta),
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      grepl("core", targetName) ~ "ATX(N) score",
      TRUE ~ NA_character_
    ),
    targetName_anno = factor(
      targetName_anno,
      level = c(
        "ATX(N) score", "pTau-217", "pTau-181", "pTau-231",
        "BD-pTau-217", "BD-pTau-181", "BD-pTau-231",
        "Abeta42/Abeta40", "ICAM1", "GDF15", "NEFL", "GFAP"
      )
    )
  ) %>% 
  group_by(diet) %>% 
  mutate(
    fdr = p.adjust(p, method = "BH"),
    sig = case_when(
      p < 0.01 ~ "#",
      p < 0.05 ~ "+",
      TRUE ~ ""
    )
  ) %>%
  arrange(targetName_anno)

amed_dat <- fig_dat %>% 
  mutate(
    diet_anno = factor(
      diet_anno,
      level = unique(
        nhs_dietcomp_dict %>% 
          filter(grepl("aMED", score)) %>% 
          pull(diet_anno)
      )
    )
  )

nhs_amed_dat <- amed_dat %>% 
  ungroup() %>% 
  filter(diet == "amed_90_cumavg") %>% 
  dplyr::select(
    targetName_anno,
    beta_nhs = beta,
    se_nhs   = se,
    p_nhs    = p,
    sig_nhs  = sig
  )

dp_amed_dat <- read_excel("results/dp_NULISA_baseline_diet.xlsx") %>%
  filter(
    targetName %in% c(
      "all_core3_ew_signed",
      core_list[["all"]][["core3"]],
      core_list[["BD"]][["core3"]]
    ),
    LOD_flag != TRUE | is.na(LOD_flag),
    diet == "amed",
    age_group == "age_ge_50",
    gender_group == "all"
  ) %>%
  mutate(
    beta = as.numeric(beta),
    se   = as.numeric(se),
    p    = as.numeric(p),
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      grepl("core", targetName) ~ "ATX(N) score",
      TRUE ~ NA_character_
    ),
    targetName_anno = factor(
      targetName_anno,
      levels = c(
        "ATX(N) score", "pTau-217", "pTau-181", "pTau-231",
        "BD-pTau-217", "BD-pTau-181", "BD-pTau-231",
        "Abeta42/Abeta40", "ICAM1", "GDF15", "NEFL", "GFAP"
      )
    )
  ) %>%
  filter(!is.na(targetName_anno)) %>%
  dplyr::select(
    targetName_anno,
    beta_dp = beta,
    se_dp   = se,
    p_dp    = p
  ) %>%
  arrange(targetName_anno)

merge_amed_dat <- merge(
  nhs_amed_dat,
  dp_amed_dat,
  by = "targetName_anno"
) %>%
  mutate(
    neglog10_pdp = ifelse(
      !is.na(p_dp) & p_dp > 0,
      -log10(p_dp),
      NA_real_
    )
  )

### Correlation

amed_cor <- cor.test(
  merge_amed_dat$beta_nhs,
  merge_amed_dat$beta_dp,
  method = "pearson"
)

x_lim <- range(merge_amed_dat$beta_nhs, na.rm = TRUE)
y_lim <- range(merge_amed_dat$beta_dp, na.rm = TRUE)
lim_all <- range(c(x_lim, y_lim), na.rm = TRUE)

## Fitted line

fit <- lm(beta_dp ~ beta_nhs, data = merge_amed_dat)

beta_fit <- coef(fit)[2]
p_fit    <- summary(fit)$coefficients[2, 4]
r_fit    <- cor(merge_amed_dat$beta_nhs, merge_amed_dat$beta_dp)

pred_df <- data.frame(
  beta_nhs = seq(
    -1,
    1,
    length.out = 200
  )
)

pred_mat <- predict(
  fit,
  newdata = pred_df,
  interval = "confidence"
)

pred_df$fit <- pred_mat[, "fit"]
pred_df$lwr <- pred_mat[, "lwr"]
pred_df$upr <- pred_mat[, "upr"]

lim_all <- range(
  c(
    merge_amed_dat$beta_nhs,
    merge_amed_dat$beta_dp
  ),
  na.rm = TRUE
)

p_scatter_amed <- ggplot(
  merge_amed_dat,
  aes(x = beta_nhs, y = beta_dp)
) +
  
  geom_ribbon(
    data = pred_df,
    aes(
      x = beta_nhs,
      ymin = lwr,
      ymax = upr
    ),
    fill = "#D9D9D9",
    alpha = 0.35,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = pred_df,
    aes(
      x = beta_nhs,
      y = fit
    ),
    color = "gray55",
    linewidth = 0.5,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    aes(
      fill = neglog10_pdp
    ),
    size = 2.6,
    shape = 21,
    stroke = 0.35,
    color = "black",
    alpha = 0.9
  ) +
  
  scale_fill_gradient(
    low  = "#7f9c6c",
    high = "#E3ECD5",
    name = expression(-log[10](italic(P)[DIRECT-PLUS]))
  ) +
  
  geom_text_repel(
    aes(label = targetName_anno),
    size = 3.2,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 1
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = -Inf,
    hjust = -0.2,
    vjust = -0.5,
    size = 3.5,
    color = "black",
    label = paste0(
      "italic(r) == ",
      round(r_fit, 2)
    ),
    parse = TRUE
  ) +
  
  labs(
    x = expression(beta["aMED"] ~ "in NHS"),
    y = expression(beta["aMED"] ~ "in DIRECT-PLUS")
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      size = 10,
      color = "black"
    ),
    legend.text = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.25,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.25,
      color = "black"
    ),
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    panel.grid = element_blank()
  ) +
  coord_equal(
    xlim = c(lim_all[1] - 0.02, lim_all[2] + 0.02),
    ylim = c(lim_all[1] - 0.02, lim_all[2] + 0.02)
  )

ggsave(
  "results/figures/exd_fig3c_amed_nhs_dp.pdf",
  p_scatter_amed,
  width = 5.5,
  height = 3.2
)

### Source data

source_data_ext3c_points <- merge_amed_dat %>%
  transmute(
    protein_label = as.character(targetName_anno),
    beta_aMED_NHS = beta_nhs,
    p_value_aMED_NHS = p_nhs,
    beta_aMED_DIRECT_PLUS = beta_dp,
    p_value_aMED_DIRECT_PLUS = p_dp,
    neg_log10_p_DIRECT_PLUS = neglog10_pdp
  ) %>%
  arrange(protein_label)

source_data_ext3c_fit <- pred_df %>%
  transmute(
    beta_aMED_NHS = beta_nhs,
    fitted_beta_aMED_DIRECT_PLUS = fit,
    lower_95CI = lwr,
    upper_95CI = upr
  )

write_xlsx(
  list(
    ext3c_points = source_data_ext3c_points,
    ext3c_fit = source_data_ext3c_fit
  ),
  "results/source_data/source_data_extended_data_fig3c.xlsx"
)

##------------------
## 3d 
##------------------

rm(list = ls())

load("data/NULISA_core.RData")
source("scripts/source_col.R")
source("scripts/source_dietcomp.R")

marker_order <- c(
  "ATX(N) score",
  "BD-pTau-217", "BD-pTau-181", "BD-pTau-231",
  "pTau-217", "pTau-181", "pTau-231",
  "Abeta42/Abeta40",
  "GFAP", "NEFL", "GDF15", "ICAM1"
)

med_marker_res <- read_excel("results/nhs_tics_mediation.xlsx") %>% 
  filter(
    diet == "MIND90_cumavg",
    targetName %in% c(
      "all_core3_ew_signed",
      core_list[["all"]][["core3"]],
      core_list[["BD"]][["core3"]]
    )
  ) %>%
  mutate(
    acme     = as.numeric(acme),
    acme_p   = as.numeric(acme_p),
    acme_lcl = as.numeric(acme_lcl),
    acme_ucl = as.numeric(acme_ucl),
    sig_signed = sign(acme) * -log10(acme_p),
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      grepl("core", targetName) ~ "ATX(N) score",
      TRUE ~ NA_character_
    ),
    targetName_anno = factor(targetName_anno, levels = marker_order)
  ) %>%
  filter(
    !is.na(targetName_anno),
    !is.na(acme),
    !is.na(acme_p)
  ) %>%
  arrange(targetName_anno)

p_mind_marker_med <- ggplot(
  med_marker_res,
  aes(x = acme, y = targetName_anno)
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray31",
    linewidth = 0.4
  ) +
  
  geom_errorbarh(
    aes(xmin = acme_lcl, xmax = acme_ucl),
    height = 0,
    linewidth = 0.45,
    color = "gray45",
    alpha = 0.85
  ) +
  
  geom_point(
    aes(fill = sig_signed),
    color = "gray31",
    shape = 21,
    size = 3.2,
    stroke = 0.25
  ) +
  
  scale_y_discrete(limits = rev(marker_order)) +
  
  scale_fill_gradient2(
    name = expression(sign(ACME) %*% -log[10](italic(P))),
    low = "#3B6FB6",
    mid = "white",
    high = "#B23A48",
    midpoint = 0
  ) +
  
  labs(
    x = "ACME",
    y = NULL
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.2, color = "black"),
    axis.ticks.y = element_blank(),
    axis.line.x  = element_line(linewidth = 0.2, color = "black"),
    axis.line.y  = element_blank(),
    legend.title = element_text(size = 9, color = "black"),
    legend.text  = element_text(size = 8, color = "black"),
    legend.position = "bottom"
  ) +
  
  guides(
    fill = guide_colorbar(
      barheight = unit(0.35, "cm"),
      barwidth  = unit(3, "cm")
    )
  )

ggsave(
  "results/figures/exd_fig3d_mind_marker_mediation.pdf",
  p_mind_marker_med,
  width = 3,
  height = 3.8
)

### Source data

source_data_ext3d_mediation <- med_marker_res %>%
  transmute(
    protein_label = as.character(targetName_anno),
    ACME = acme,
    ACME_lower_95CI = acme_lcl,
    ACME_upper_95CI = acme_ucl,
    ACME_p_value = acme_p,
    signed_neg_log10_p = sig_signed
  ) %>%
  arrange(protein_label)

write_xlsx(
  list(
    ext3d_mediation = source_data_ext3d_mediation
  ),
  "results/source_data/source_data_extended_data_fig3d.xlsx"
)

##---------------------------------------------------------
## Exd Data Fig 4
##---------------------------------------------------------

##------------------
## 4a
##------------------

rm(list = ls())

load("data/NULISA_core.RData")

mri_var <- "InferiorLateralVentricleTotalVolume"

marker_order <- c(
  "pTau-217", "pTau-181", "pTau-231",
  "BD-pTau-217", "BD-pTau-181", "BD-pTau-231",
  "Abeta42/Abeta40",
  "ICAM1", "GDF15", "NEFL", "GFAP"
)

forest_dat <- read_excel("results/dp_MRI_NULISA_baseline.xlsx") %>%
  filter(
    scale == "scaled",
    grepl("all_core3_ew_signed", targetName),
    grepl(paste0(mri_var, "_change"), MRI),
    !(group == "all" & MRI != paste0(mri_var, "_time0")),
    gender_group == "all",
    age_group == "age_ge_50"
  ) %>% 
  mutate(
    p     = as.numeric(p),
    beta  = as.numeric(beta),
    se    = as.numeric(se),
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    group = factor(group, levels = c("all", "HDG", "MED", "GreenMED"))
  )

forest_lim_min <- min(forest_dat$lower)
forest_lim_max <- max(forest_dat$upper)

heatmap_dat <- read_excel("results/dp_MRI_NULISA_baseline.xlsx") %>%
  filter(
    scale == "scaled",
    targetName %in% c(core_list[["all"]][["core3"]], core_list[["BD"]][["core3"]]),
    grepl(paste0(mri_var, "_change"), MRI),
    !(group == "all" & MRI != paste0(mri_var, "_time0")),
    gender_group == "all",
    age_group == "age_ge_50"
  ) %>%
  mutate(
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      TRUE ~ NA_character_
    ),
    targetName_anno = factor(
      targetName_anno,
      levels = marker_order
    ),
    p    = as.numeric(p),
    beta = as.numeric(beta),
    sig = case_when(
      p < 0.05 ~ "+",
      TRUE ~ ""
    ),
    group = factor(group, levels = c("all", "HDG", "MED", "GreenMED"))
  ) %>%
  arrange(targetName_anno)

heatmap_lim_min <- min(heatmap_dat$beta, na.rm = TRUE)
heatmap_lim_max <- max(heatmap_dat$beta, na.rm = TRUE)

p_forest1 <- ggplot(
  data = forest_dat %>% filter(MRI == paste0(mri_var, "_change")),
  aes(y = group, x = beta, xmin = lower, xmax = upper)
) +
  geom_vline(
    xintercept = 0,
    color = "gray31",
    linetype = "dashed",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_errorbarh(
    size = 0.3,
    height = 0.5,
    color = "black"
  ) +
  geom_point(
    aes(fill = -log10(p)),
    shape = 21,
    size = 3,
    stroke = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  scale_fill_gradient(
    low  = "#334f24",
    high = "#E3ECD5",
    name = expression(-log[10](italic(P)))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9, hjust = 1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.5),
    axis.text.x = element_text(color = "black", size = 9),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(4, "mm"),
    title = element_text(colour = "black", size = 9)
  ) +
  scale_y_discrete(limits = rev)

p_heatmap1 <- ggplot(
  heatmap_dat %>% filter(MRI == paste0(mri_var, "_change")),
  aes(x = targetName_anno, y = group)
) +
  geom_tile(aes(fill = beta), color = "gray31", linewidth = 0.3) +
  geom_text(aes(label = sig), color = "black", size = 5) +
  scale_x_discrete(limits = marker_order, position = "top") +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#7EA6C7",
    mid = "white",
    high = "#E09A5C",
    midpoint = 0,
    limits = c(heatmap_lim_min, heatmap_lim_max)
  ) +
  labs(fill = expression(beta["components"])) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, color = "black", size = 9),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(4, "mm"),
    panel.border = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  coord_fixed(ratio = 1)

p_final <- (p_forest1 | p_heatmap1) +
  plot_layout(widths = c(0.8, 2))

ggsave(
  "results/figures/exd_fig4a_LVV.pdf",
  p_final,
  width = 7,
  height = 2.2
)

p_forest_bar <- p_forest1 +
  guides(
    fill = guide_colorbar(
      barheight = unit(0.4, "cm"),
      barwidth  = unit(2.5, "cm")
    )
  ) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0)
  )

forest_bar <- cowplot::get_legend(p_forest_bar)

ggsave(
  "results/figures/exd_fig4a_LVV_forest_legend.pdf",
  cowplot::ggdraw(forest_bar),
  width = 2,
  height = 2
)

p_heatmap_bar <- p_heatmap1 +
  guides(
    fill = guide_colorbar(
      barheight = unit(0.4, "cm"),
      barwidth  = unit(2.5, "cm")
    )
  ) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0)
  )

heatmap_bar <- cowplot::get_legend(p_heatmap_bar)

ggsave(
  "results/figures/exd_fig4a_LVV_heatmap_legend.pdf",
  cowplot::ggdraw(heatmap_bar),
  width = 2,
  height = 2
)

### Source data

source_data_ext4a_forest <- forest_dat %>%
  transmute(
    intervention_group = as.character(group),
    beta,
    lower_95CI = lower,
    upper_95CI = upper,
    p_value = p,
    neg_log10_p = -log10(p)
  ) %>%
  arrange(intervention_group)

source_data_ext4a_heatmap <- heatmap_dat %>%
  transmute(
    intervention_group = as.character(group),
    protein_label = as.character(targetName_anno),
    beta,
    p_value = p
  ) %>%
  arrange(intervention_group, protein_label)

write_xlsx(
  list(
    ext4a_forest = source_data_ext4a_forest,
    ext4a_heatmap = source_data_ext4a_heatmap
  ),
  "results/source_data/source_data_extended_data_fig4a.xlsx"
)

##------------------
## 4b 
##------------------

load("data/processed_data_dp_MRI.RData")

baseline_data <- npq_all %>% 
  filter(time == "A") %>% 
  left_join(mri_raw, by=c("sno","group")) %>% 
  left_join(allgeno, by="sno", suffix = c("","_genetic")) %>% 
  left_join(direct_diet, by="sno") %>% 
  filter(!is.na(group),
         !is.na(HOMAIR_0)) %>% 
  mutate(group = factor(group, levels=c("HDG","MED","GreenMED"))) %>%
  mutate(
    dplyr::across(
      all_of(c(protall_list, scoreall_list)),
      ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  )

source("scripts/source_col.R")

mri_var <- "InferiorLateralVentricleTotalVolume"
mri_anno <- "LVV"

fig_data <- baseline_data %>% 
  filter(age >= 50) %>% 
  mutate(all_core3_ew_signed_t = ntile(all_core3_ew_signed, 3)) %>% 
  select(
    sno, group, age_visit, gender, all_core3_ew_signed_t,
    time0  = .data[[paste0(mri_var, "_time0")]],
    time18 = .data[[paste0(mri_var, "_time18")]]
  )

fig_data_t0 <- fig_data %>% 
  mutate(
    time = "Baseline",
    mri_y = time0
  ) %>% 
  select(sno, group, age_visit, gender, all_core3_ew_signed_t, time, mri_y)

fig_data_t18 <- fig_data %>% 
  mutate(
    time = "Month 18",
    mri_y = time18
  ) %>% 
  select(sno, group, age_visit, gender, all_core3_ew_signed_t, time, mri_y)

fig_data_comb <- bind_rows(fig_data_t0, fig_data_t18) %>%
  mutate(
    time = factor(time, levels = c("Baseline", "Month 18")),
    group = factor(group, levels = c("HDG", "MED", "GreenMED")),
    all_core3_ew_signed_t = factor(
      all_core3_ew_signed_t,
      levels = c(1, 2, 3),
      labels = c("Low ATX(N) score", "Medium ATX(N) score", "High ATX(N) score")
    )
  ) %>%
  filter(!is.na(mri_y)) 

fig_data_comb <- fig_data_comb %>% 
  mutate(
    mri_y_resid = residuals(
      lm(mri_y ~ age_visit + gender, data = fig_data_comb, na.action = na.exclude)
    )
  )

p_box <- ggplot(
  fig_data_comb,
  aes(
    x = time,
    y = mri_y_resid,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.75
    ),
    color = "gray31",
    alpha = 0.5,
    stroke = 0.3,
    size = 1,
    shape = 21,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "Intervention group",
    values = c(
      "HDG"      = "#E4C97A",
      "MED"      = "#A7C2DD",
      "GreenMED" = "#9BB596"
    )
  ) +
  facet_wrap(
    ~ all_core3_ew_signed_t,
    nrow = 1
  ) +
  labs(
    x = NULL,
    y = paste0(mri_anno, " (residual)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.line  = element_line(linewidth = 0.3, color = "black"),
    legend.position = "top",
    
    strip.background = element_rect(
      color = "black",
      fill = "white",
      linewidth = 0.6
    ),
    
    strip.text = element_text(
      size = 10,
      color = "black"
    )
  )

ggsave(
  "results/figures/exd_fig4b_LVV.pdf",
  p_box,
  width = 6,
  height = 3
)

wide_dat <- fig_data_comb %>%
  select(sno, group, all_core3_ew_signed_t, time, mri_y_resid) %>%
  pivot_wider(names_from = time, values_from = mri_y_resid) %>%
  filter(!is.na(Baseline), !is.na(`Month 18`))

stat_df <- wide_dat %>%
  group_by(all_core3_ew_signed_t, group) %>% 
  dplyr::summarise(
    p = t.test(`Month 18`, Baseline, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    group1 = "Baseline",
    group2 = "Month 18",
    p.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )

write.csv(
  stat_df,
  "results/figures/exd_fig4b_LVV_paired_tests.csv",
  row.names = FALSE
)

### Source data

source_data_ext4b_points <- fig_data_comb %>%
  transmute(
    pair_id = as.integer(factor(sno)),
    intervention_group = as.character(group),
    ATXN_score_tertile = as.character(all_core3_ew_signed_t),
    time = as.character(time),
    LVV_residual = mri_y_resid
  ) %>%
  arrange(ATXN_score_tertile, intervention_group, pair_id, time)

source_data_ext4b_paired_tests <- stat_df %>%
  transmute(
    ATXN_score_tertile = as.character(all_core3_ew_signed_t),
    intervention_group = as.character(group),
    group1,
    group2,
    p_value = p
  )

write_xlsx(
  list(
    ext4b_points = source_data_ext4b_points,
    ext4b_paired_tests = source_data_ext4b_paired_tests
  ),
  "results/source_data/source_data_extended_data_fig4b.xlsx"
)
