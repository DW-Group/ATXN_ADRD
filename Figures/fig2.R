library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(readxl)
library(cowplot)
library(scam)
library(survival)
library(timeROC)
library(writexl)

##---------------------------------------------------------
## Fig 2a
##---------------------------------------------------------

rm(list=ls())

source("scripts/source_rm_markers.R")

fig_dat <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(scale == "scaled",
         LOD_flag != TRUE | is.na(LOD_flag),
         !(targetName %in% ipc_bcg_fails)) %>%
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
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH"),
    neglog10p = -log10(p),
    sig_cat = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      p < 0.05   ~ "p < 0.05",
      TRUE       ~ "NS"
    )
  )

### Plot

lab_dat <- fig_dat %>% filter(p < 0.05)
yline <- -log10(0.05)
p_dementia <- ggplot(fig_dat, aes(x = beta, y = neglog10p)) +
  geom_hline(yintercept = yline, linetype = 2, color = "gray31", linewidth = 0.4) +
  annotate(
    "text",
    x = max(fig_dat$beta, na.rm = TRUE),
    y = yline,
    label = "italic(P) == 0.05",
    parse = TRUE,
    color = "gray31",
    vjust = 1.5, hjust = 0.8, size = 3.2
  ) +
  geom_point(aes(shape = sig_cat, alpha = sig_cat, fill = sig_cat, color = sig_cat), size = 2.5, stroke = 0.4) +
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
    values = c("FDR < 0.05" = 0.95, "p < 0.05" = 0.95, "NS" = 0.5),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(
      "FDR < 0.05" = "#6E597B",  
      "p < 0.05"   = "#b0cba9",  
      "NS"         = "#BDBDBD"  
    )
  ) +
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
    x = "log(hazard ratio)",
    y = expression(-log[10](italic(P))),
    title = "Dementia risk"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size=11, color = "black"),
    axis.text.y = element_text(size=11, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    axis.title = element_text(size=12, color="black"),
    plot.title = element_text(size=12)
  )

ggsave("results/figures/fig2a.pdf", p_dementia, width = 3.5, height = 3.5)

p_leg <- ggplot(
  fig_dat %>% filter(sig_cat %in% c("FDR < 0.05", "p < 0.05")),
  aes(x = beta, y = neglog10p)
) +
  geom_point(
    aes(shape = sig_cat, fill = sig_cat, color = sig_cat),
    size = 2,
    stroke = 0.4
  ) +
  scale_shape_manual(
    values = c("FDR < 0.05" = 23, "p < 0.05" = 22)
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = "black", "p < 0.05" = "black")
  ) +
  scale_fill_manual(
    values = c("FDR < 0.05" = "#6E597B", "p < 0.05" = "#b0cba9")
  ) +
  guides(
    shape = guide_legend(
      title = NULL,
      override.aes = list(
        alpha = 1,
        fill  = c("#6E597B", "#b0cba9"),
        color = "black"
      )
    ),
    fill  = "none",
    color = "none",
    alpha = "none"
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(color = "black")
  )

legend_sig <- cowplot::get_legend(p_leg)
ggsave(
  "results/figures/fig2a_legend.pdf",
  cowplot::ggdraw(legend_sig),
  width = 1.2,
  height = 0.6
)

### Source data

source_data_fig2a <- fig_dat %>%
  transmute(
    target_name_anno = targetName_anno,
    log_hazard_ratio = beta,
    p_value = p,
    fdr = fdr
  ) %>%
  arrange(desc(p_value), target_name_anno)

write_xlsx(
  list(source_data_fig2a = source_data_fig2a),
  "results/source_data/source_data_fig2a.xlsx"
)

##---------------------------------------------------------
## Fig 2b
##---------------------------------------------------------

rm(list = ls())

source("scripts/source_rm_markers.R")

##-----------------------
## Prospective cohort
##-----------------------

pros_dat <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(
    scale == "scaled",
    LOD_flag != TRUE | is.na(LOD_flag),
    !(targetName %in% ipc_bcg_fails)
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
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH")
  ) %>%
  filter(p < 0.05) %>%
  transmute(
    targetName,
    targetName_anno,
    beta_pros = beta,
    p_pros = p,
    fdr_pros = fdr
  )

##-----------------------
## Nested case-control
##-----------------------

caco_dat <- read_excel("results/nhs_caco_dementia_NULISA.xlsx") %>%
  filter(
    scale == "scaled",
    LOD_flag != TRUE | is.na(LOD_flag),
    !(targetName %in% ipc_bcg_fails)
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
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH")
  ) %>%
  transmute(
    targetName,
    targetName_anno,
    beta_caco = beta,
    p_caco = p,
    fdr_caco = fdr
  )

##-----------------------
## Merge
##-----------------------

plot_dat <- pros_dat %>%
  inner_join(caco_dat, by = c("targetName", "targetName_anno")) %>%
  mutate(
    sig_caco = case_when(
      fdr_caco < 0.05 ~ "FDR < 0.05",
      p_caco < 0.05   ~ "p < 0.05",
      TRUE            ~ "NS"
    ),
    dir_concord = sign(beta_pros) == sign(beta_caco),
    neglog10_pcaco = ifelse(!is.na(p_caco) & p_caco > 0, -log10(p_caco), NA_real_)
  )

##-----------------------
## Plot
##-----------------------

x_lim <- range(plot_dat$beta_pros, na.rm = TRUE)
y_lim <- range(plot_dat$beta_caco, na.rm = TRUE)

x_pad <- 0.08 * diff(x_lim)
y_pad <- 0.18 * diff(y_lim)

x_lim_plot <- x_lim + c(-x_pad, x_pad)
y_lim_plot <- y_lim + c(-y_pad, y_pad)

p_beta_compare <- ggplot(plot_dat, aes(x = beta_pros, y = beta_caco)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray35", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray35", linewidth = 0.45) +
  geom_point(
    aes(fill = neglog10_pcaco),
    shape = 21,
    size = 2.6,
    stroke = 0.35,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_gradient(
    low  = "#7f9c6c",
    high = "#E3ECD5",
    name = expression(-log[10](italic(P)))
  ) +
  geom_text_repel(
    data = plot_dat %>% filter(p_pros < 0.05),
    aes(label = targetName_anno),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.18,
    segment.size = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 1
  ) +
  coord_cartesian(
    xlim = x_lim_plot,
    ylim = y_lim_plot
  ) +
  labs(
    x = "log(hazard ratio): NHS prospective cohort",
    y = "log(odds ratio): NHS nested case-control"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, color = "black"),
    legend.text  = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.title = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    plot.subtitle = element_text(size = 10, color = "black"),
    aspect.ratio = 1
  )

ggsave(
  "results/figures/fig2b_beta_compare_prospective_vs_caco.pdf",
  p_beta_compare,
  width = 5.5,
  height = 3.5
)

### Source data

source_data_fig2b <- plot_dat %>%
  transmute(
    target_name_anno = targetName_anno,
    log_hazard_ratio_NHS_prospective = beta_pros,
    p_value_NHS_prospective = p_pros,
    log_odds_ratio_NHS_nested_case_control = beta_caco,
    p_value_NHS_nested_case_control = p_caco
  ) %>%
  arrange(target_name_anno)

write_xlsx(
  list(source_data_fig2b = source_data_fig2b),
  "results/source_data/source_data_fig2b.xlsx"
)

##---------------------------------------------------------
## Fig 2c
##---------------------------------------------------------

rm(list=ls())

source("scripts/source_rm_markers.R")

fig_dat <- read_excel("results/nhs_tics_NULISA.xlsx") %>%
  filter(scale == "scaled",
         LOD_flag != TRUE | is.na(LOD_flag),
         !(targetName %in% ipc_bcg_fails),
         outcome == "av_ntotal") %>%
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
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH"),
    neglog10p = -log10(p),
    sig_cat = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      p < 0.05   ~ "p < 0.05",
      TRUE       ~ "NS"
    )
  )

### Plot

lab_dat <- fig_dat %>% filter(p < 0.05)
yline <- -log10(0.05)
p_tics <- ggplot(fig_dat, aes(x = beta, y = neglog10p)) +
  geom_hline(yintercept = yline, linetype = 2, color = "gray31", linewidth = 0.4) +
  annotate(
    "text",
    x = max(fig_dat$beta, na.rm = TRUE),
    y = yline,
    label = "italic(P) == 0.05",
    parse = TRUE,
    color = "gray31",
    vjust = -0.5, hjust = 0.8, size = 3.2
  ) +
  geom_point(aes(shape = sig_cat, alpha = sig_cat, fill = sig_cat, color = sig_cat), size = 2.5, stroke = 0.4) +
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
    values = c("FDR < 0.05" = 0.95, "p < 0.05" = 0.95, "NS" = 0.75),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(
      "FDR < 0.05" = "#6E597B",  
      "p < 0.05"   = "#b0cba9",  
      "NS"         = "#BDBDBD"  
    )
  ) +
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
    title = "Cognitive function (global cognition)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size=11, color = "black"),
    axis.text.y = element_text(size=11, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    axis.title = element_text(size=12, color="black"),
    plot.title = element_text(size=12)
  )

ggsave("results/figures/fig2c.pdf", p_tics, width = 3.5, height = 3.5)

### Source data

source_data_fig2c <- fig_dat %>%
  transmute(
    target_name_anno = targetName_anno,
    beta_coefficient = beta,
    p_value = p,
    fdr = fdr
  ) %>%
  arrange(desc(p_value), target_name_anno)

write_xlsx(
  list(source_data_fig2c = source_data_fig2c),
  "results/source_data/source_data_fig2c.xlsx"
)

##---------------------------------------------------------
## Fig 2d
##---------------------------------------------------------

rm(list = ls())

source("scripts/source_rm_markers.R")

dementia_sig <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(
    scale == "scaled",
    LOD_flag != TRUE | is.na(LOD_flag),
    !(targetName %in% ipc_bcg_fails),
    p < 0.05
  ) %>%
  pull(targetName) %>%
  unique()

##-----------------------
## Prospective cohort
##-----------------------

pros_dat <- read_excel("results/nhs_tics_NULISA.xlsx") %>%
  filter(
    scale == "scaled",
    LOD_flag != TRUE | is.na(LOD_flag),
    !(targetName %in% ipc_bcg_fails),
    outcome == "av_ntotal"
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
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH")
  ) %>%
  filter(p < 0.05) %>%
  transmute(
    targetName,
    targetName_anno,
    beta_pros = beta,
    p_pros = p,
    fdr_pros = fdr
  )

##-----------------------
## Nested case-control
##-----------------------

caco_dat <- read_excel("results/nhs_caco_tics_NULISA.xlsx") %>%
  filter(
    scale == "scaled",
    LOD_flag != TRUE | is.na(LOD_flag),
    !(targetName %in% ipc_bcg_fails),
    outcome == "av_ntotal",
    subset == "all"
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
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH")
  ) %>%
  transmute(
    targetName,
    targetName_anno,
    beta_caco = beta,
    p_caco = p,
    fdr_caco = fdr
  )

##-----------------------
## Merge
##-----------------------

plot_dat <- pros_dat %>%
  inner_join(caco_dat, by = c("targetName", "targetName_anno")) %>%
  mutate(
    sig_caco = case_when(
      fdr_caco < 0.05 ~ "FDR < 0.05",
      p_caco < 0.05   ~ "p < 0.05",
      TRUE            ~ "NS"
    ),
    dir_concord = sign(beta_pros) == sign(beta_caco),
    neglog10_pcaco = ifelse(!is.na(p_caco) & p_caco > 0, -log10(p_caco), NA_real_),
    sig_dementia = targetName %in% dementia_sig
  )

##-----------------------
## Plot
##-----------------------

x_lim <- range(plot_dat$beta_pros, na.rm = TRUE)
y_lim <- range(plot_dat$beta_caco, na.rm = TRUE)

x_pad <- 0.08 * diff(x_lim)
y_pad <- 0.12 * diff(y_lim)

x_lim_plot <- x_lim + c(-x_pad, x_pad)
y_lim_plot <- y_lim + c(-y_pad, y_pad)

p_beta_compare_tics <- ggplot(
  plot_dat,
  aes(x = beta_pros, y = beta_caco)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray35", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray35", linewidth = 0.45) +
  geom_point(
    aes(fill = neglog10_pcaco, shape = sig_dementia),
    size = 2.6,
    stroke = 0.35,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_gradient(
    low  = "#7f9c6c",
    high = "#E3ECD5",
    name = expression(-log[10](italic(P)))
  ) +
  scale_shape_manual(
    values = c(
      `FALSE` = 21,
      `TRUE`  = 22
    ),
    name = "Significant for dementia risk in NHS (P < 0.05)"
  ) +
  geom_text_repel(
    data = plot_dat %>% filter(p_pros < 0.05),
    aes(label = targetName_anno),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.18,
    segment.size = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 1
  ) +
  coord_cartesian(
    xlim = x_lim_plot,
    ylim = y_lim_plot
  ) +
  labs(
    x = expression(beta["marker"]*" NHS prospective cohort"),
    y = expression(beta["marker"]*" NHS nested case-control")
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, color = "black"),
    legend.text  = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.title = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    plot.subtitle = element_text(size = 10, color = "black"),
    aspect.ratio = 1
  )

ggsave(
  "results/figures/fig2d_beta_compare_prospective_vs_caco_tics.pdf",
  p_beta_compare_tics,
  width = 7,
  height = 3.5
)

### Source data

source_data_fig2d <- plot_dat %>%
  transmute(
    target_name_anno = targetName_anno,
    beta_NHS_prospective_global_cognition = beta_pros,
    p_value_NHS_prospective_global_cognition = p_pros,
    beta_NHS_nested_case_control_global_cognition = beta_caco,
    p_value_NHS_nested_case_control_global_cognition = p_caco
  ) %>%
  arrange(target_name_anno)

write_xlsx(
  list(source_data_fig2d = source_data_fig2d),
  "results/source_data/source_data_fig2d.xlsx"
)

##---------------------------------------------------------
## Fig 2e
##---------------------------------------------------------

rm(list=ls())

nhs_tics_dat <- read_excel("results/nhs_tics_NULISA.xlsx") %>%
  filter(scale == "scaled",
         LOD_flag != TRUE | is.na(LOD_flag),
         outcome == "av_ntotal",
         p < 0.05)

nhs_tics_all <- nhs_tics_dat

nhs_dementia_dat <- read_excel("results/nhs_dementia_NULISA.xlsx") %>%
  filter(scale == "scaled",
         LOD_flag != TRUE | is.na(LOD_flag),
         p < 0.05)

dp_dat <- read_excel("results/dp_MoCA_NULISA_cross_Total_z.xlsx") 
dp_sub <- dp_dat %>%
  filter(scale == "scaled",
         group == "all",
         gender_group == "all",
         age_group == "age_ge_50",
         age_covar == "linear",
         cog == "Total_z",
         LOD_flag != TRUE | is.na(LOD_flag)) %>%
  mutate(
    targetName_anno = case_when(
      !is.na(targetName_orig) ~ targetName_orig,
      targetName == "Ab42_div_Ab40_log2" ~ "Aβ42/Aβ40",
      targetName == "AT_comp" ~ "Amyloid-Tau score",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(targetName_anno),
         targetName %in% c(nhs_tics_all$targetName)) %>%
  mutate(
    p    = as.numeric(p),
    beta = as.numeric(beta),
    fdr  = p.adjust(p, method = "BH"),
    neglog10p = -log10(p),
    sig_cat = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      p < 0.05   ~ "p < 0.05",
      TRUE       ~ "NS"
    )
  )

nhs_fig <- nhs_tics_all %>% 
  mutate(beta_nhs = as.numeric(beta)) %>% 
  select(targetName, beta_nhs, p_nhs = p)

dp_fig <- dp_sub %>% 
  mutate(sig_dementia = ifelse(targetName %in% nhs_dementia_dat$targetName | targetName == "all_core3_ew_signed", T, F)) %>% 
  select(targetName, targetName_anno, beta_dp = beta, p_dp = p, sig_dementia)

fig_dat <- merge(dp_fig, nhs_fig, by="targetName") %>%
  mutate(
    neglog10_pdp = ifelse(!is.na(p_dp) & p_dp > 0, -log10(p_dp), NA_real_),
    sig_dementia = as.logical(sig_dementia)
  )

cor_res <- cor.test(fig_dat$beta_dp, fig_dat$beta_nhs)

### Plot

p_quad <- ggplot(fig_dat, aes(beta_nhs, beta_dp)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray35", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray35", linewidth = 0.45) +
  
  geom_point(
    aes(fill = neglog10_pdp, shape = sig_dementia),
    size = 2.6, stroke = 0.4,
    color = "black"
  ) +
  
  scale_fill_gradient(
    low  = "#7f9c6c",   
    high = "#E3ECD5",
    name = expression(-log[10](italic(P)[DIRECT-PLUS]))
  ) +
  
  ggrepel::geom_text_repel(
    data = fig_dat,
    aes(label = targetName_anno),
    size = 3.5,
    max.overlaps = Inf,        
    box.padding = 0.35,
    point.padding = 0.25,
    force = 1.2,              
    seed = 1
  ) +
  scale_shape_manual(
    values = c(
      `FALSE` = 21,  
      `TRUE`  = 22   
    ),
    name = "Significant for dementia risk in NHS (P < 0.05)"
  ) +
  
  labs(
    x = expression(beta["marker"]*" NHS global cognition z-score"),
    y = expression(beta["marker"]*" DIRECT-PLUS MoCA total z-score")
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11, color = "black"),
    axis.line  = element_line(linewidth = 0.25, color = "black"),
    axis.title = element_text(size = 11.5, color = "black"),
    axis.text  = element_text(size = 11, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black")
  )

ggsave("results/figures/fig2e.pdf", p_quad, width = 7, height = 3.5)

### Source data

source_data_fig2e <- fig_dat %>%
  transmute(
    target_name_anno = targetName_anno,
    beta_NHS_global_cognition_z_score = beta_nhs,
    p_value_NHS_global_cognition = p_nhs,
    beta_DIRECT_PLUS_MoCA_total_z_score = beta_dp,
    p_value_DIRECT_PLUS = p_dp,
    significant_for_dementia_risk_in_NHS = sig_dementia
  ) %>%
  arrange(target_name_anno)

write_xlsx(
  list(source_data_fig2e = source_data_fig2e),
  "results/source_data/source_data_fig2e.xlsx"
)

##---------------------------------------------------------
## Fig 2f
##---------------------------------------------------------

rm(list=ls())

source("scripts/source_col.R")

fig_dat <- read_excel("results/nhs_tics_NULISA_R2.xlsx") %>% 
  filter(outcome == "av_ntotal") %>%
  mutate(
    targetName_anno = case_when(
      targetName == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      targetName == "core_joint" ~ "ATX(N)",
      grepl("_",targetName) ~ gsub("_","-",targetName),
      TRUE ~ targetName
    ),
    group = case_when(targetName_anno %in% c("ATX(N)") ~ "ATX(N)",
                      targetName_anno %in% c("Abeta42/Abeta40") ~ "A (Amyloid)",
                      targetName_anno %in% c("pTau-217","pTau-181","pTau-231","BD-pTau-217","BD-pTau-181","BD-pTau-231") ~ "T (Tau)",
                      targetName_anno %in% c("ICAM1","GDF15") ~ "X (Others)",
                      targetName_anno %in% c("NEFL","GFAP") ~ "N (Neurodegeneration)"),
    group = factor(group, levels=c("ATX(N)","A (Amyloid)","T (Tau)","N (Neurodegeneration)","X (Others)")),
    delta_r2_pct = delta_r2*100) %>% 
  arrange(-delta_r2_pct)%>%
  mutate(targetName_anno = factor(targetName_anno, levels = targetName_anno),)

### Plot

p_r2 <- ggplot(fig_dat, aes(x = targetName_anno, y = delta_r2_pct, fill = group)) +
  geom_col(width = 0.9, color = "gray31", linewidth = 0.2) +
  geom_text(
  aes(
    y = delta_r2_pct + 1,   
    label = paste0(targetName_anno, "\n", sprintf("%.2f%%", delta_r2_pct))
    ),
    size = 4,
    color = "black"
  ) +
  coord_polar(start = -0.23) +
    scale_y_continuous(
      name   = expression(Delta*R^2~"(%)"),
      limits = c(-10, 25)
    ) +
  labs(
    x = NULL,
    fill = "ATX(N) domain"
  ) +
  scale_fill_manual(values = atxn_col) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_blank()
  )

ggsave("results/figures/fig2f.pdf", p_r2, width = 10, height = 10)

### Source data

source_data_fig2f <- fig_dat %>%
  transmute(
    target_name_anno = as.character(targetName_anno),
    ATXN_domain = as.character(group),
    delta_R2_percent = delta_r2_pct
  )

write_xlsx(source_data_fig2f, "results/source_data/source_data_fig2f.xlsx")

##---------------------------------------------------------
## Fig 2g
##---------------------------------------------------------

rm(list=ls())

load("data/processed_data_nhs.RData")
load("data/NULISA_core.RData")

vars_all <- c(
  "tad_final", "adcase",
  "age_final", "fhdem", "dep_antidep",
  "APOE4_carrier", "PGS002280_scaled_pool_rsq0.3",
  core_list[["all"]][["core3"]]
)

model_data_auc <- model_data %>%
  dplyr::select(all_of(vars_all)) %>%
  tidyr::drop_na()

tmp_mod1 <- coxph(
  Surv(tad_final, adcase) ~ age_final,
  ties = "breslow",
  data = model_data_auc
)

tmp_mod2 <- coxph(
  Surv(tad_final, adcase) ~ age_final + fhdem + dep_antidep,
  ties = "breslow",
  data = model_data_auc
)

tmp_mod3 <- coxph(
  Surv(tad_final, adcase) ~ age_final + fhdem + dep_antidep +
    APOE4_carrier + `PGS002280_scaled_pool_rsq0.3`,
  ties = "breslow",
  data = model_data_auc
)

tmp_mod4 <- coxph(
  as.formula(paste(
    "Surv(tad_final, adcase) ~",
    paste(
      c(
        "age_final", "fhdem", "dep_antidep",
        "APOE4_carrier", "`PGS002280_scaled_pool_rsq0.3`",
        core_list[["all"]][["core3"]]
      ),
      collapse = " + "
    )
  )),
  ties = "breslow",
  data = model_data_auc
)

## Harrell's C with 95% CI from summary(coxph)

get_cindex_ci <- function(fit, label) {
  conc <- summary(fit)$concordance
  cval <- unname(conc[1])
  se   <- unname(conc[2])
  
  data.frame(
    model = label,
    c_index = cval,
    lower = cval - 1.96 * se,
    upper = cval + 1.96 * se
  )
}

c_tab <- bind_rows(
  get_cindex_ci(tmp_mod1, "Age"),
  get_cindex_ci(tmp_mod2, "+ Family history + Depression"),
  get_cindex_ci(tmp_mod3, "+ APOE4 status + ADRD PRS"),
  get_cindex_ci(tmp_mod4, "+ ATX(N)")
)

print(c_tab)

lp1 <- predict(tmp_mod1, type = "lp")
lp2 <- predict(tmp_mod2, type = "lp")
lp3 <- predict(tmp_mod3, type = "lp")
lp4 <- predict(tmp_mod4, type = "lp")

t_grid <- as.numeric(quantile(
  model_data_auc$tad_final[model_data_auc$adcase == 1],
  probs = seq(0.1, 0.9, by = 0.1),
  na.rm = TRUE
))
t_grid <- sort(unique(t_grid))

roc1 <- timeROC(
  T = model_data_auc$tad_final,
  delta = model_data_auc$adcase,
  marker = lp1,
  cause = 1,
  weighting = "marginal",
  times = t_grid,
  iid = TRUE
)

roc2 <- timeROC(
  T = model_data_auc$tad_final,
  delta = model_data_auc$adcase,
  marker = lp2,
  cause = 1,
  weighting = "marginal",
  times = t_grid,
  iid = TRUE
)

roc3 <- timeROC(
  T = model_data_auc$tad_final,
  delta = model_data_auc$adcase,
  marker = lp3,
  cause = 1,
  weighting = "marginal",
  times = t_grid,
  iid = TRUE
)

roc4 <- timeROC(
  T = model_data_auc$tad_final,
  delta = model_data_auc$adcase,
  marker = lp4,
  cause = 1,
  weighting = "marginal",
  times = t_grid,
  iid = TRUE
)

auc_df <- data.frame(
  time = t_grid,
  Model1 = roc1$AUC,
  Model2 = roc2$AUC,
  Model3 = roc3$AUC,
  Model4 = roc4$AUC
)

auc_long <- tidyr::pivot_longer(
  auc_df,
  cols = -time,
  names_to = "model",
  values_to = "auc"
)

auc_long$model <- factor(
  auc_long$model,
  levels = c("Model1", "Model2", "Model3", "Model4"),
  labels = c(
    "Age",
    "+ Family history + Depression",
    "+ APOE4 status + ADRD PRS",
    "+ ATX(N)"
  )
)

model_cols <- c(
  "Age" = "#9DBDD2",
  "+ Family history + Depression" = "#7F9C6C",
  "+ APOE4 status + ADRD PRS" = "#DEB956",
  "+ ATX(N)" = "#C47F4E"
)

text_cols <- c(
  "Age" = "#5E7F97",
  "+ Family history + Depression" = "#4F6942",
  "+ APOE4 status + ADRD PRS" = "#A67E1E",
  "+ ATX(N)" = "#8A542B"
)

c_anno <- data.frame(
  model = factor(c_tab$model, levels = levels(auc_long$model)),
  c_index = c_tab$c_index,
  lower = c_tab$lower,
  upper = c_tab$upper
) %>%
  left_join(
    auc_long %>%
      group_by(model) %>%
      filter(time == max(time)) %>%
      ungroup(),
    by = "model"
  ) %>%
  mutate(
    label = sprintf("C = %.2f (%.2f-%.2f)", c_index, lower, upper),
    x_pos = max(auc_long$time) - diff(range(auc_long$time)) * 0.02,
    y_pos = auc
  )

### Plot

p_auc <- ggplot(auc_long, aes(x = time, y = auc, color = model)) +
  geom_line(aes(color = model), linewidth = 0.6) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.4
  ) +
  geom_text(
    data = c_anno,
    aes(
      x = x_pos,
      y = y_pos,
      label = label
    ),
    color = text_cols[as.character(c_anno$model)],
    inherit.aes = FALSE,
    hjust = 1,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = model_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Time (month)",
    y = "Time-dependent AUC"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12, color = "black"),
    axis.text  = element_text(size = 11, color = "black"),
    axis.line  = element_line(color = "black", linewidth=0.4),
    axis.ticks = element_line(color = "black", linewidth=0.4)
  )

ggsave("results/figures/fig2g_cindex.pdf", p_auc, width = 6, height = 4)

### Source data

source_data_fig2g_auc_curves <- auc_long %>%
  transmute(
    time_months = time,
    model = as.character(model),
    time_dependent_AUC = auc
  ) %>%
  arrange(model, time_months)

source_data_fig2g_c_index_labels <- c_anno %>%
  transmute(
    model = as.character(model),
    Harrell_C_index = c_index,
    lower_95_CI = lower,
    upper_95_CI = upper,
    label_shown = label
  ) %>%
  arrange(model)

write_xlsx(source_data_fig2g_auc_curves, "results/source_data/source_data_fig2g.xlsx")

##---------------------------------------------------------
## Fig 2h
##---------------------------------------------------------

rm(list=ls())

fig_dat <- read_excel("results/nhs_caco_dementia_NULISA_auc_for_fig.xlsx") %>%
  mutate(model_anno = case_when(
           model == "m1" ~ paste0("m1: Age (AUC: ", auc_only,")"),
           model == "m2" ~ paste0("m2: + Family history \n      + Depression (AUC: ", auc_only,")"),
           model == "m3" ~ paste0("m3: + APOE4 status + ADRD PRS (AUC: ", auc_only,")"),
           model == "m4" ~ paste0("m4: + ATX(N) (AUC: ", auc_only,")")
         ),
         model_anno = factor(model_anno, levels = unique(model_anno)))
fig_dat$model_anno

### Plot

p_auc <- ggplot(fig_dat, aes(x=fpr, y=tpr, color=model_anno, group=model_anno)) +
  geom_line(size=0.6) +
  geom_abline(linetype="dashed", color="grey70") +
  scale_color_manual(values = c("#9dbdd2", "#7f9c6c", "#deb956", "#c47f4e")) +
  coord_equal() +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    color = "Model"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size=12, colour="black"),
        axis.text = element_text(size=11, colour="black"),
        plot.title = element_text(size=12.5, hjust=0.5, colour="black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.ticks.x = element_line(color = "black", linewidth = 0.2),
        axis.ticks.y = element_line(color = "black", linewidth = 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.background = element_blank(),
        legend.position = c(0.65, 0.18)) +
  guides(colour = guide_legend(keyheight = 1.5))

ggsave("results/figures/fig2h_auc.pdf", p_auc, width = 3.5, height = 3.5)

rm(list = ls())

importance <- read_excel("results/nhs_caco_dementia_NULISA_importance_for_fig.xlsx")

dict_predictor <- as.data.frame(rbind(
  c("age_final", "Age", "base"),
  c("dep_antidep", "Depression", "base"),
  c("fhdem", "Family history", "base"),
  c("APOE4_carrier", "APOE4 status", "gene"),
  c("PGS002280_scaled_pool_rsq0.3", "ADRD PRS", "gene"),
  c("GDF15", "GDF15", "marker"),
  c("GFAP", "GFAP", "marker"),
  c("ICAM1", "ICAM1", "marker"),
  c("NEFL", "NEFL", "marker"),
  c("pTau_181", "pTau-181", "marker"),
  c("pTau_217", "pTau-217", "marker"),
  c("pTau_231", "pTau-231", "marker"),
  c("Ab42_div_Ab40_log2", "Abeta42/Abeta40", "marker")
))
colnames(dict_predictor) <- c("predictor", "anno_predictor", "cat_predictor")

importance_anno <- importance %>%
  left_join(dict_predictor, by = "predictor") %>%
  filter(!is.na(anno_predictor)) %>%
  arrange(overall) %>%
  mutate(
    anno_predictor = factor(anno_predictor, levels = anno_predictor)
  )

### Plot

p_importance_tile <- ggplot(importance_anno) +

  geom_tile(
    aes(x = 1, y = anno_predictor, fill = overall),
    width = 0.42,
    height = 0.95,
    color = "white",
    linewidth = 0.5
  ) +

  geom_text(
    aes(x = 1, y = anno_predictor, label = anno_predictor),
    size = 3.5,
    color = "black"
  ) +

  scale_fill_gradient(
    low = "#FBFAF7",
    high = "#C8BDB3",
    breaks = c(0.02, 0.06, 0.10),
    name = "Mean absolute\nSHAP value",
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.16, "cm"),
      barheight = unit(1.15, "cm")
    )
  ) +

  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.79, 1.21), clip = "off") +

  labs(x = NULL, y = NULL) +

  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6.8),
    legend.text = element_text(size = 6.5),
    legend.box.margin = margin(-2, 0, 0, 0),
    plot.margin = margin(1, 1, 1, 1)
  )

ggsave(
  "results/figures/fig2h_shap_importance_tile.pdf",
  p_importance_tile,
  width = 1.15,
  height = 3.7
)

### Source data

source_data_fig2h_roc_curves <- fig_dat %>%
  transmute(
    model = as.character(model_anno),
    false_positive_rate = fpr,
    true_positive_rate = tpr
  ) %>%
  arrange(model, false_positive_rate)

write_xlsx(
  list(source_data_fig2h_roc_curves = source_data_fig2h_roc_curves),
  "results/source_data/source_data_fig2h_roc.xlsx"
)

### Source data

source_data_fig2h_shap_importance <- importance_anno %>%
  transmute(
    predictor = as.character(anno_predictor),
    mean_absolute_SHAP_value = overall
  )

write_xlsx(source_data_fig2h_shap_importance, "results/source_data/source_data_fig2h_shap_importance.xlsx")

##---------------------------------------------------------
## Fig 2i
##---------------------------------------------------------

rm(list = ls())

##-------------------------------
## Load data
##-------------------------------

load("data/processed_data_dp_MRI.RData")
load("data/NULISA_core.RData")

marker_dict <- read_excel("data/NULISA_dictionary.xlsx")

##-------------------------------
## Baseline data
##-------------------------------

baseline_data <- npq_all %>% 
  filter(time == "A") %>% 
  left_join(mri_raw, by = c("sno", "group")) %>% 
  filter(
    !is.na(group),
    !is.na(HOC_mean_time0)
  ) %>% 
  mutate(
    group = case_when(
      group == "HDG" ~ "HDG",
      group == "MedDiet" ~ "MED",
      group == "MedDiet+Poly" ~ "GreenMED",
      TRUE ~ as.character(group)
    ),
    group = factor(group, levels = c("HDG", "MED", "GreenMED"))
  )

##-------------------------------
## MRI pseudotime
##-------------------------------

baseline_data <- baseline_data %>%
  mutate(
    MRI_burden_raw = -HOC_mean_time0,
    MRI_burden_z   = as.numeric(scale(MRI_burden_raw)),
    MRI_pseudotime = (MRI_burden_z - min(MRI_burden_z, na.rm = TRUE)) /
      (max(MRI_burden_z, na.rm = TRUE) - min(MRI_burden_z, na.rm = TRUE))
  )

markers_to_plot <- unique(c(core_list[["all"]][["core3"]], core_list[["BD"]][["core3"]]))
markers_to_plot <- c(markers_to_plot, "all_core3_ew_signed")
markers_to_plot <- markers_to_plot[markers_to_plot %in% colnames(baseline_data)]

##-------------------------------
## Z-score biomarkers
##-------------------------------

plot_long <- baseline_data %>%
  dplyr::select(sno, group, MRI_pseudotime, all_of(markers_to_plot)) %>%
  pivot_longer(
    cols = all_of(markers_to_plot),
    names_to = "marker",
    values_to = "value"
  ) %>%
  filter(!is.na(value), !is.na(MRI_pseudotime)) %>%
  group_by(marker) %>%
  mutate(
    value_z   = as.numeric(scale(value))
  ) %>%
  ungroup() %>%
  left_join(
    marker_dict %>% dplyr::select(targetName, targetName_anno),
    by = c("marker" = "targetName")
  ) %>%
  mutate(
    marker_lab = case_when(
      marker == "Ab42_div_Ab40_log2" ~ "Abeta42/Abeta40",
      marker == "all_core3_ew_signed" ~ "ATX(N) score",
      TRUE ~ targetName_anno
    )
  )

marker_lab_levels <- plot_long %>%
  distinct(marker, marker_lab) %>%
  mutate(marker = factor(marker, levels = markers_to_plot)) %>%
  arrange(marker) %>%
  pull(marker_lab)

plot_long <- plot_long %>%
  mutate(
    marker_lab = factor(marker_lab, levels = marker_lab_levels)
  )

##-------------------------------
## Fit smooth curves
##-------------------------------

x_grid <- seq(0, 1, length.out = 200)

fit_one_marker <- function(dat) {
  fit <- scam(
    value_z ~ s(MRI_pseudotime, bs="mpi", k = 8),
    data = dat
  )
  
  newdat <- data.frame(MRI_pseudotime = x_grid)
  pred <- predict(fit, newdata = newdat, se.fit = TRUE)
  
  tibble(
    MRI_pseudotime = x_grid,
    fit = pred$fit,
    lwr = pred$fit - 1.96 * pred$se.fit,
    upr = pred$fit + 1.96 * pred$se.fit
  )
}

pred_df <- plot_long %>%
  group_by(marker_lab) %>%
  group_modify(~ fit_one_marker(.x)) %>%
  ungroup()

bottom_df <- plot_long %>%
  distinct(sno, group, MRI_pseudotime) %>%
  mutate(
    y_strip = case_when(
      group == "HDG" ~ -0.48,
      group == "MED" ~ -0.63,
      group == "GreenMED" ~ -0.78
    )
  )

curve_cols <- c(
  "ATX(N) score" = "#C47F4E",
  
  "pTau-217"     = "#4F769E",
  "pTau-231"     = "#6F93B8",
  "pTau-181"     = "#97B3CF",
  
  "NEFL"         = "#7FA35F",
  "GFAP"         = "#AFC59A",
  
  "ICAM1"        = "#D1B25C",
  "GDF15"        = "#A98558",
  
  "Abeta42/Abeta40" = "#7F8C97",
  
  "BD-pTau-217"  = "#3F6287",
  "BD-pTau-231"  = "#5D82AA",
  "BD-pTau-181"  = "#86A5C3"
)

line_types <- c(
  "ATX(N) score" = "solid",
  
  "pTau-217"     = "solid",
  "pTau-231"     = "solid",
  "pTau-181"     = "solid",
  
  "BD-pTau-217"  = "22",
  "BD-pTau-231"  = "22",
  "BD-pTau-181"  = "22",
  
  "NEFL"         = "solid",
  "GFAP"         = "solid",
  "ICAM1"        = "solid",
  "GDF15"        = "solid",
  "Abeta42/Abeta40"    = "solid"
)

marker_order <- c(
  "ATX(N) score",
  "pTau-217",
  "pTau-231",
  "pTau-181",
  "BD-pTau-217",
  "BD-pTau-231",
  "BD-pTau-181",
  "Abeta42/Abeta40",
  "ICAM1",
  "GDF15",
  "NEFL",
  "GFAP"
)

curve_cols <- curve_cols[names(curve_cols) %in% levels(plot_long$marker_lab)]

plot_long$marker_lab <- factor(plot_long$marker_lab, levels = marker_order)
pred_df$marker_lab   <- factor(pred_df$marker_lab, levels = marker_order)

##-------------------------------
## Plot
##-------------------------------

p_pseudotime <- ggplot() +
  geom_line(
    data = pred_df %>% filter(marker_lab != "ATX(N) score"),
    aes(
      x = MRI_pseudotime,
      y = fit,
      color = marker_lab,
      linetype = marker_lab
    ),
    linewidth = 0.28,
    alpha = 0.7
  ) +
  geom_line(
    data = pred_df %>% filter(marker_lab == "ATX(N) score"),
    aes(
      x = MRI_pseudotime,
      y = fit,
      color = marker_lab,
      linetype = marker_lab
    ),
    linewidth = 0.5,
    alpha = 1
  ) +
  scale_color_manual(
    values = curve_cols,
    breaks = marker_order
  ) +
  
  scale_linetype_manual(
    values = line_types,
    breaks = marker_order
  ) +
  geom_rug(
    data = plot_long %>% distinct(sno, MRI_pseudotime),
    aes(x = MRI_pseudotime),
    sides = "b",
    inherit.aes = FALSE,
    alpha = 0.18,
    linewidth = 0.18,
    color = "grey45"
  ) +
  scale_fill_manual(values = curve_cols) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00")
  ) +
  labs(
    x = "MRI structural burden",
    y = "NULISA protein quantification (z-score)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(size = 8.5, color = "black"),
    axis.title = element_text(size = 8.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black"),
    axis.line  = element_line(linewidth = 0.25, color = "black"),
    legend.position = c(0.01, 1.05),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.3, "cm"),
    legend.key.width = grid::unit(0.55, "cm"),
    plot.margin = margin(8, 8, 8, 8),
    panel.grid = element_blank()
  ) +
  guides(
    linetype = "none",  
    color = guide_legend(
      override.aes = list(linewidth = 0.6)
    )
  )

ggsave(
  "results/figures/fig2i_mri_pseudotime_atxn.pdf",
  p_pseudotime,
  width = 2.3,
  height = 3.3
)

### Source data

source_data_fig2i_curves <- pred_df %>%
  transmute(
    marker = as.character(marker_lab),
    MRI_structural_burden = MRI_pseudotime,
    fitted_z_score = fit
  ) %>%
  arrange(marker, MRI_structural_burden)

source_data_fig2i_rug <- plot_long %>%
  distinct(sno, MRI_pseudotime) %>%
  transmute(
    MRI_structural_burden = MRI_pseudotime
  ) %>%
  arrange(MRI_structural_burden)

write_xlsx(
  list(
    source_data_fig2i_curves = source_data_fig2i_curves,
    source_data_fig2i_rug = source_data_fig2i_rug
  ),
  "results/source_data/source_data_fig2i.xlsx"
)

