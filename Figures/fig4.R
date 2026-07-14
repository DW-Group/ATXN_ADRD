library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(readxl)
library(writexl)
library(cowplot)
library(patchwork)

##---------------------------------------------------------
## Fig 4a 
##---------------------------------------------------------

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
    MIND90_cumavg,
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
    MIND90_cumavg,
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

pred_df <- data.frame(
  MIND90_cumavg = seq(
    min(plot_dat$MIND90_cumavg, na.rm = TRUE),
    max(plot_dat$MIND90_cumavg, na.rm = TRUE),
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
  aes(x = MIND90_cumavg, y = all_core3_ew_signed_resid)
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
    aes(x = MIND90_cumavg, ymin = lwr, ymax = upr),
    fill = "#6E8FA8",
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = MIND90_cumavg, y = fit),
    color = "#6E8FA8",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = MIND90_cumavg, y = lwr),
    color = "#6E8FA8",
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = MIND90_cumavg, y = upr),
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
    x = "MIND score",
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
  "results/figures/fig4a_mind_scatter_pooled.pdf",
  p_scatter_pool,
  width = 2.6,
  height = 2.8
)

### Source data

source_data_fig4a_points <- plot_dat %>%
  transmute(
    cohort = as.character(cohort),
    MIND_score = MIND90_cumavg,
    ATXN_score_residual = all_core3_ew_signed_resid
  ) %>%
  arrange(cohort, MIND_score)

source_data_fig4a_fit <- pred_df %>%
  transmute(
    MIND_score = MIND90_cumavg,
    fitted_ATXN_score_residual = fit,
    lower_95_CI = lwr,
    upper_95_CI = upr
  ) %>%
  arrange(MIND_score)

write_xlsx(
  list(
    fig4a_points = source_data_fig4a_points,
    fig4a_fit = source_data_fig4a_fit
  ),
  "results/source_data/source_data_fig4a.xlsx"
)

##---------------------------------------------------------
## Fig 4b
##---------------------------------------------------------

rm(list=ls())

### Data

load("data/NULISA_core.RData")
source("scripts/source_dietcomp.R")

fig_dat <- read_excel("results/nhs_meta_NULISA_diet.xlsx") %>%
  left_join(nhs_dietcomp_dict, by="diet") %>%
  filter(
    targetName %in% c("all_core3_ew_signed", core_list[["all"]][["core3"]], core_list[["BD"]][["core3"]]),
    LOD_flag != TRUE | is.na(LOD_flag),
    !(grepl("comp", diet) & grepl("mutual", model)),
    grepl("MIND", score)
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
        "ATX(N) score","pTau-217","pTau-181","pTau-231",
        "BD-pTau-217","BD-pTau-181","BD-pTau-231",
        "Abeta42/Abeta40","ICAM1","GDF15","NEFL","GFAP"
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

wins_top_p <- 0
wins_bot_p <- 0.02

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

### MIND heatmap

mind_dat <- fig_dat %>% 
  mutate(
    diet_anno = factor(
      diet_anno,
      level = unique(nhs_dietcomp_dict %>% filter(grepl("MIND", score)) %>% pull(diet_anno))
    )
  )

nhs_mind_dat <- mind_dat %>% 
  ungroup() %>% 
  filter(diet == "MIND90_cumavg") %>% 
  dplyr::select(targetName_anno, beta_nhs = beta, se_nhs = se, p_nhs = p, sig_nhs = sig)

mind_top <- mind_dat %>% filter(targetName_anno == "ATX(N) score")
mind_bot <- mind_dat %>% filter(targetName_anno != "ATX(N) score")

### Plot

p_mind_top <- ggplot(mind_top, aes(x = diet_anno, y = targetName_anno)) +
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

p_mind_bot <- ggplot(mind_bot, aes(x = diet_anno, y = targetName_anno)) +
  geom_tile(aes(fill = beta), color = "gray42", linewidth = 0.2) +
  geom_text(aes(label = sig), color = "black", size = 3) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
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

final_plot <- p_mind_top / plot_spacer() / p_mind_bot +
  plot_layout(heights = c(1, 0.05, length(unique(mind_bot$targetName_anno))))

ggsave("results/figures/fig4b_meta.pdf", final_plot, height = 4.5)

p_atxn_bar <- p_mind_top +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0)
  )
atxn_bar <- cowplot::get_legend(p_atxn_bar)
ggsave(
  "results/figures/fig4b_meta_colorbar_ATXN.pdf",
  cowplot::ggdraw(atxn_bar),
  width = 2,
  height = 2
)

p_comp_bar <- p_mind_bot +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0)
  )
comp_bar <- cowplot::get_legend(p_comp_bar)
ggsave(
  "results/figures/fig4b_meta_colorbar_components.pdf",
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
  "results/figures/fig4b_meta_sig_legend.pdf",
  p_sig_leg,
  width  = 1.5,
  height = 0.5
)

### Source data

source_data_fig4b_top <- mind_top %>%
  ungroup() %>%
  transmute(
    dietary_variable = as.character(diet_anno),
    target = as.character(targetName_anno),
    beta,
    p_value = p
  ) %>%
  arrange(dietary_variable, target)

source_data_fig4b_components <- mind_bot %>%
  ungroup() %>%
  transmute(
    dietary_variable = as.character(diet_anno),
    target = as.character(targetName_anno),
    beta,
    p_value = p
  ) %>%
  arrange(dietary_variable, target)

write_xlsx(
  list(
    fig4b_atxn_score = source_data_fig4b_top,
    fig4b_components = source_data_fig4b_components
  ),
  "results/source_data/source_data_fig4b.xlsx"
)

##---------------------------------------------------------
## Fig 4c
##---------------------------------------------------------

### Data

dp_mind_dat <- read_excel("results/dp_NULISA_baseline_diet.xlsx") %>%
  filter(
    targetName %in% c(
      "all_core3_ew_signed",
      core_list[["all"]][["core3"]],
      core_list[["BD"]][["core3"]]
    ),
    LOD_flag != TRUE | is.na(LOD_flag),
    diet == "MIND14",
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

merge_mind_dat <- merge(
  nhs_mind_dat,
  dp_mind_dat,
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

mind_cor <- cor.test(
  merge_mind_dat$beta_nhs,
  merge_mind_dat$beta_dp,
  method = "pearson"
)

x_lim <- range(merge_mind_dat$beta_nhs, na.rm = TRUE)
y_lim <- range(merge_mind_dat$beta_dp, na.rm = TRUE)
lim_all <- range(c(x_lim, y_lim), na.rm = TRUE)

## Fitted line

fit <- lm(beta_dp ~ beta_nhs, data = merge_mind_dat)

beta_fit <- coef(fit)[2]
p_fit    <- summary(fit)$coefficients[2, 4]
r_fit    <- cor(merge_mind_dat$beta_nhs, merge_mind_dat$beta_dp)

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
    merge_mind_dat$beta_nhs,
    merge_mind_dat$beta_dp
  ),
  na.rm = TRUE
)

### Plot

p_scatter_mind <- ggplot(
  merge_mind_dat,
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
    x = expression(beta["MIND"] ~ "in NHS"),
    y = expression(beta["MIND"] ~ "in DIRECT-PLUS")
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
  "results/figures/fig4c_mind.pdf",
  p_scatter_mind,
  width = 5.5,
  height = 3.2
)

### Source data

source_data_fig4c_points <- merge_mind_dat %>%
  transmute(
    target = as.character(targetName_anno),
    beta_MIND_NHS = beta_nhs,
    p_value_NHS = p_nhs,
    beta_MIND_DIRECT_PLUS = beta_dp,
    p_value_DIRECT_PLUS = p_dp,
    neg_log10_p_DIRECT_PLUS = neglog10_pdp
  ) %>%
  arrange(target)

source_data_fig4c_fit <- pred_df %>%
  transmute(
    beta_MIND_NHS = beta_nhs,
    fitted_beta_DIRECT_PLUS = fit,
    lower_95_CI = lwr,
    upper_95_CI = upr
  ) %>%
  arrange(beta_MIND_NHS)

write_xlsx(
  list(
    fig4c_points = source_data_fig4c_points,
    fig4c_fit = source_data_fig4c_fit
  ),
  "results/source_data/source_data_fig4c.xlsx"
)

##---------------------------------------------------------
## Fig 4d
##---------------------------------------------------------

# BioRender

##---------------------------------------------------------
## Fig 4e
##---------------------------------------------------------

rm(list=ls())

load("data/NULISA_core.RData")
source("scripts/source_col.R")
source("scripts/source_dietcomp.R")

med_res <- read_excel("results/nhs_tics_mediation.xlsx") %>% 
  filter(
    targetName == "BD_pTau_217",
    diet != "MIND90_cumavg"
  ) %>%
  left_join(nhs_dietcomp_dict, by = "diet") %>% 
  mutate(
    acme   = as.numeric(acme),
    acme_p = as.numeric(acme_p),
    acme_lcl = as.numeric(acme_lcl),
    acme_ucl = as.numeric(acme_ucl),
    sig_col = -log10(acme_p),
    sig_signed = sign(acme) * -log10(acme_p),
    diet_anno = factor(
      diet_anno,
      levels = unique(nhs_dietcomp_dict %>% filter(grepl("MIND", score)) %>% pull(diet_anno))
    )
  ) %>%
  filter(!is.na(diet_anno), !is.na(acme), !is.na(acme_p))

### Plot

p_lolipop_med <- ggplot(med_res, aes(x = diet_anno, y = acme)) +
  geom_segment(
    aes(xend = diet_anno, y = 0, yend = acme),
    linewidth = 0.3,
    alpha = 0.8,
    color = "gray45"
  ) +
  geom_point(
    aes(fill = sig_col),
    color = "gray31",
    shape = 21,
    size = 3.2,
    stroke = 0.25
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray31", linewidth = 0.4) +
  scale_fill_gradient(
    name = expression(-log[10](italic(P)[ACME])),
    low = "grey95",
    high = "#8F7BC8"
  ) +
  labs(
    x = NULL,
    y = "ACME"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8.5, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.2, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.2, color = "black"),
    axis.line.x  = element_line(linewidth = 0.2, color = "black"),
    axis.line.y  = element_line(linewidth = 0.2, color = "black"),
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
  "results/figures/fig4_mind_component_mediation_lolipop.pdf",
  p_lolipop_med,
  width = 3.8,
  height = 3.8
)

### Source data

source_data_fig4e <- med_res %>%
  transmute(
    dietary_component = as.character(diet_anno),
    mediator = "BD-pTau-217",
    ACME = acme,
    lower_95_CI = acme_lcl,
    upper_95_CI = acme_ucl,
    p_value_ACME = acme_p,
    neg_log10_p_ACME = sig_col
  ) %>%
  arrange(dietary_component)

write_xlsx(
  list(fig4e_mediation = source_data_fig4e),
  "results/source_data/source_data_fig4e.xlsx"
)

##---------------------------------------------------------
## Fig 4f
##---------------------------------------------------------

rm(list=ls())

mri_var <- "HOC_mean"

load("data/NULISA_core.RData")

forest_dat <- read_excel("results/dp_MRI_NULISA_baseline.xlsx") %>%
  filter(scale=="scaled",
         grepl("all_core3_ew_signed",targetName),
         grepl(paste0(mri_var,"_change"), MRI),
         !(group=="all"& MRI!=paste0(mri_var,"_time0")),
         gender_group=="all",
         age_group=="age_ge_50") %>% 
  mutate(
    p    = as.numeric(p),
    beta = as.numeric(beta),
    se = as.numeric(se),
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    group = factor(group, levels=c("all","HDG","MED","GreenMED"))
  )

forest_lim_min <- min(forest_dat$lower)
forest_lim_max <- max(forest_dat$upper)

heatmap_dat <- read_excel("results/dp_MRI_NULISA_baseline.xlsx") %>%
  filter(scale=="scaled",
         targetName %in% c(core_list[["all"]][["core3"]],core_list[["BD"]][["core3"]]),
         grepl(paste0(mri_var,"_change"), MRI),
         !(group=="all"& MRI!=paste0(mri_var,"_time0")),
         gender_group=="all",
         age_group=="age_ge_50") %>%
  mutate(
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      TRUE ~ NA_character_
    ),
    targetName_anno = factor(targetName_anno, level=c("pTau-217","pTau-181","pTau-231","BD-pTau-217","BD-pTau-181","BD-pTau-231","Abeta42/Abeta40",
                                                      "ICAM1","GDF15","NEFL","GFAP")),
    p    = as.numeric(p),
    beta = as.numeric(beta),
    sig = case_when(
      p < 0.01 ~ "#",
      p < 0.05 ~ "+",
      TRUE ~ ""
    ),
    group = factor(group, levels=c("all","HDG","MED","GreenMED"))
  ) %>%
  arrange(targetName_anno)

heatmap_lim_min <- min(heatmap_dat$beta)
heatmap_lim_max <- max(heatmap_dat$beta)

### Plot

p_forest1 <- ggplot(data=forest_dat %>% filter(MRI==paste0(mri_var,"_change")), aes(y=group, x=beta, xmin=lower, xmax=upper)) +
  geom_vline(xintercept=0, color='gray31', linetype='dashed', alpha=.5, linewidth =0.5) +
  geom_errorbarh(size=0.3,height=0.5, color="black")+
  geom_point(aes(fill = -log10(p)), shape=21, size=3, stroke = 0.3, position=position_dodge(width = 0.5)) +
  scale_fill_gradient(
    low  = "#334f24",   
    high = "#E3ECD5",  
    name = expression(-log[10](italic(P)))
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
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
        title=element_text(colour = "black", size = 9)) +
  scale_y_discrete(limits=rev)

p_heatmap1 <- ggplot(heatmap_dat %>% filter(MRI==paste0(mri_var,"_change")), aes(x = targetName_anno, y = group)) +
  geom_tile(aes(fill = beta), color = "gray31", linewidth = 0.3) +
  geom_text(aes(label = sig), color = "black", size = 5) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  scale_fill_gradient2(low = "#7EA6C7", mid = "white", high = "#E09A5C", midpoint = 0, limits = c(heatmap_lim_min, heatmap_lim_max)) +
  labs(fill = expression(beta["components"])) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 0, vjust = 0, color = "black", size=9),
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

ggsave(paste0("results/figures/fig4f_",mri_var,".pdf"), p_final, width = 7, height = 2.2)

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
  paste0("results/figures/fig4f_forest_",mri_var,"_legend.pdf"),
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
  paste0("results/figures/fig4f_heatmap_",mri_var,"_legend.pdf"),
  cowplot::ggdraw(heatmap_bar),
  width = 2,
  height = 2
)

### Source data

source_data_fig4f_forest <- forest_dat %>%
  filter(MRI == paste0(mri_var, "_change")) %>%
  transmute(
    intervention_group = as.character(group),
    target = "ATX(N) score",
    beta,
    lower_95_CI = lower,
    upper_95_CI = upper,
    p_value = p,
    neg_log10_p = -log10(p)
  ) %>%
  arrange(factor(intervention_group, levels = c("all", "HDG", "MED", "GreenMED")))

source_data_fig4f_heatmap <- heatmap_dat %>%
  filter(MRI == paste0(mri_var, "_change")) %>%
  transmute(
    intervention_group = as.character(group),
    target = as.character(targetName_anno),
    beta,
    p_value = p,
    significance_symbol = sig
  ) %>%
  arrange(
    factor(intervention_group, levels = c("all", "HDG", "MED", "GreenMED")),
    factor(
      target,
      levels = c(
        "pTau-217", "pTau-181", "pTau-231",
        "BD-pTau-217", "BD-pTau-181", "BD-pTau-231",
        "Abeta42/Abeta40", "ICAM1", "GDF15", "NEFL", "GFAP"
      )
    )
  )

write_xlsx(
  list(
    fig4f_forest = source_data_fig4f_forest,
    fig4f_heatmap = source_data_fig4f_heatmap
  ),
  "results/source_data/source_data_fig4f.xlsx"
)

##---------------------------------------------------------
## Fig 4g
##---------------------------------------------------------

source("scripts/source_col.R")

mri_var <- "HOC_mean"
mri_anno <- "HOC"

fig_data <- baseline_data %>% 
  filter(age >= 50) %>% 
  mutate(all_core3_ew_signed_t = ntile(all_core3_ew_signed, 3)) %>% 
  select(
    sno, group, age_visit, gender, all_core3_ew_signed_t,
    time0  = .data[[paste0(mri_var, "_time0")]],
    time18 = .data[[paste0(mri_var, "_time18")]]
  )

fig_data_t0 <- fig_data %>% 
  mutate(time = "Baseline",
         mri_y = time0) %>% 
  select(sno, group, age_visit, gender, all_core3_ew_signed_t, time, mri_y)

fig_data_t18 <- fig_data %>% 
  mutate(time = "Month 18",
         mri_y = time18) %>% 
  select(sno, group, age_visit, gender, all_core3_ew_signed_t, time, mri_y)

fig_data_comb <- bind_rows(fig_data_t0, fig_data_t18) %>%
  mutate(
    time = factor(time, levels = c("Baseline", "Month 18")),
    group = factor(group, levels = c("HDG", "MED", "GreenMED")),
    all_core3_ew_signed_t = factor(all_core3_ew_signed_t, levels = c(1, 2, 3),
                                   labels = c("Low ATX(N) score", "Medium ATX(N) score", "High ATX(N) score"))
  ) %>%
  filter(!is.na(mri_y)) 

fig_data_comb <- fig_data_comb %>% 
  mutate(mri_y_resid = residuals(lm(mri_y ~ age_visit + gender, data = fig_data_comb, na.action = na.exclude)))

### Plot

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

ggsave(paste0("results/figures/fig4g_",mri_var,".pdf"), p_box, width = 6, height = 3)

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
      TRUE          ~ "ns"
    )
  )

### Source data

source_data_fig4g_points <- fig_data_comb %>%
  mutate(pair_id = as.integer(factor(sno))) %>%
  transmute(
    pair_id,
    ATXN_score_tertile = as.character(all_core3_ew_signed_t),
    intervention_group = as.character(group),
    time = as.character(time),
    MRI_residual = mri_y_resid
  ) %>%
  arrange(
    factor(
      ATXN_score_tertile,
      levels = c("Low ATX(N) score", "Medium ATX(N) score", "High ATX(N) score")
    ),
    factor(intervention_group, levels = c("HDG", "MED", "GreenMED")),
    pair_id,
    factor(time, levels = c("Baseline", "Month 18"))
  )

source_data_fig4g_pvalues <- stat_df %>%
  transmute(
    ATXN_score_tertile = as.character(all_core3_ew_signed_t),
    intervention_group = as.character(group),
    comparison = "Month 18 vs Baseline",
    paired_t_test_p_value = p
  ) %>%
  arrange(
    factor(
      ATXN_score_tertile,
      levels = c("Low ATX(N) score", "Medium ATX(N) score", "High ATX(N) score")
    ),
    factor(intervention_group, levels = c("HDG", "MED", "GreenMED"))
  )

write_xlsx(
  list(
    fig4g_points = source_data_fig4g_points,
    fig4g_pvalues = source_data_fig4g_pvalues
  ),
  "results/source_data/source_data_fig4g.xlsx"
)
