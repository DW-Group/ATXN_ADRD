library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(readxl)
library(ggridges)
library(ggrepel)
library(scales)
library(cowplot)
library(patchwork)
library(vegan)
library(ggnewscale)
library(scam)
library(writexl)

##---------------------------------------------------------
## Fig 1a
##---------------------------------------------------------

# BioRender

##---------------------------------------------------------
## Fig 1b
##---------------------------------------------------------

rm(list=ls())

##-----------------------------
## Data
##-----------------------------

load("data/data_NULISA_nhs.RData")
merged_data <- fread("data/data_dementia_tics.csv") 
marker_dict_anno <- read_excel("data/NULISA_dictionary.xlsx")
source("scripts/source_col.R")

prot_select <- prot_list[!grepl("raw", prot_list)]

tmp <- npq_all %>%
  filter(!qc_flag,
         id %in% (merged_data %>% filter(!is.na(adcase)) %>% pull(id))) %>% 
  dplyr::select(all_of(prot_select))

marker_median <- tmp %>%
  dplyr::summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "marker",
    values_to = "median"
  )
summary(marker_median$median)

marker_cv_tech <- fread("results/DP_IPC_BCG_marker_technical_CV.csv") %>%
  transmute(
    marker = targetName,
    cv_tech = technical_cv
  )
summary(marker_cv_tech$cv_tech)

nhs_data <- marker_median %>%
  left_join(marker_cv, by = "marker") %>%
  left_join(marker_cv_tech, by = "marker") %>%
  left_join(lod_counts, by = c("marker" = "targetName")) %>% 
  left_join(marker_dict_anno, by = c("marker" = "targetName")) %>%
  mutate(
    prop_below_lod = ifelse(marker %in% c("APOE", "APOE4", "CRP"), NA_real_, prop_below_lod),
    hallmark_anno = factor(hallmark_anno, levels = names(hallmark_col)),
  ) %>%
  arrange(hallmark_anno, uniprotID, commonName, targetName_anno) %>%
  mutate(
    targetName_anno = factor(targetName_anno, levels = unique(targetName_anno))
  )

r_max  <- max(nhs_data$median, na.rm = TRUE)

r_hole <- r_max * 0.6    
r_ring <- r_max * 0.06
gap_r  <- r_max * 0.06
r0     <- r_hole + r_ring + gap_r

axis_gap  <- 6.0
bar_width <- 1.00

nhs_plot <- nhs_data %>%
  mutate(
    x_raw = as.numeric(targetName_anno),
    x_num = x_raw + axis_gap
  )

n_bar <- max(nhs_plot$x_raw)
max_x <- n_bar + axis_gap

x_axis <- 0.75 + axis_gap / 2

start_angle <- -pi/2 - 2*pi*(x_axis - 0.5) / max_x

label_df <- nhs_plot %>%
  mutate(
    angle0 = 90 - 360 * (x_num - 0.5) / max_x,
    angle  = ifelse(angle0 < -90, angle0 + 280, angle0 + 100),
    hjust  = ifelse(angle0 < -90, 1, 0)
  )

r_label <- r0 + r_max * 1.20

r_outer <- r0 + r_max
y_lim   <- r_outer * 1.20

breaks_bar <- pretty_breaks(n = 5)(c(0, r_max))
breaks_rad <- r0 + breaks_bar

tick_len <- 0.50
lab_off  <- 0.70

axis_df <- tibble(x = x_axis, y0 = r0, y1 = r_outer)

tick_df <- tibble(
  x   = x_axis,
  x1  = x_axis + tick_len,
  y   = breaks_rad,
  lab = sprintf("%.2f", breaks_bar)
)

ring_df <- tibble(y = breaks_rad)

### Plot

p_circle <- ggplot(nhs_plot, aes(x = x_num)) +
  
  geom_hline(
    data = ring_df,
    aes(yintercept = y),
    color = "grey75",
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  
  geom_tile(
    aes(y = r_hole + r_ring/2, fill = hallmark_anno),
    height = r_ring,
    width  = bar_width,
    show.legend = F
  ) +
  scale_fill_manual(values = hallmark_col, name = "Hallmarks") +
  ggnewscale::new_scale_fill() +
  
  geom_segment(
    aes(xend = x_num, y = r0, yend = r0 + median, color = hallmark_anno),
    linewidth = 0.3,
    color = "gray75"
  ) +
  geom_point(
    aes(y = r0 + median, size = cv_tech, fill = prop_below_lod * 100),
    shape = 21,
    stroke = 0.45,
    color = "gray31"
  ) +
  geom_text(
    data = tick_df,
    aes(x = x_axis + lab_off, y = y+5, label = lab),
    inherit.aes = FALSE,
    size = 3.2,
    hjust = 0,
    color = "gray31"
  ) +
  geom_text(
    data = label_df,
    aes(x = x_num, y = r_label, label = as.character(targetName_anno), color=hallmark_anno,
        angle = angle, hjust = hjust),
    inherit.aes = FALSE,
    size = 2.6,
    color = "black",
    show.legend = F
  ) +
  
  scale_size_continuous(
    range = c(0.6, 6),
    name = "Coefficient of variation"
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "#F0F0F0",  
      "#F2C66D",  
      "#B65C5C", 
      "#3B0F70"   
    ),
    values = scales::rescale(c(0, 0.2, 0.6, 1)),
    name = "% < LOD",
    na.value = "white",
    guide = guide_colourbar(
      barwidth  = unit(35, "mm"),
      barheight = unit(3, "mm")
    )
  ) +
  scale_color_manual(values = hallmark_col) +
  
  coord_polar(theta = "x", start = start_angle, clip = "off") +
  scale_x_continuous(limits = c(0.5, max_x + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, y_lim * 1.02), expand = c(0, 0)) +
  
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.key.size = unit(0.9, "lines"),
    plot.margin = margin(10, 35, 10, 35)
  )

ggsave("results/figures/fig1b.pdf",
       p_circle, width = 7.5, height = 7.5)

legend_df <- data.frame(
  hallmark_anno = factor(names(hallmark_col), levels = names(hallmark_col)),
  x = 1
)
p_hallmark_legend <- ggplot(legend_df, aes(x = x, fill = hallmark_anno)) +
  geom_bar() +
  scale_fill_manual(values = hallmark_col, name = "Hallmarks") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  ) +
  guides(
    fill = guide_legend(ncol = 1, byrow = TRUE)
  )

leg <- cowplot::get_legend(p_hallmark_legend)

pdf("results/figures/fig1b_legend.pdf", width = 6.5, height = 10)
grid::grid.newpage()
grid::grid.draw(leg)
dev.off()

### Source data

source_data_fig1b <- nhs_plot %>%
  transmute(
    target_name_anno = as.character(targetName_anno),
    hallmark_anno = as.character(hallmark_anno),
    median = median,
    technical_cv = cv_tech,
    percent_below_lod = prop_below_lod * 100
  ) %>%
  arrange(plot_order)

write_xlsx(
  list(source_data_fig1b = source_data_fig1b),
  "results/source_data/source_data_fig1b.xlsx"
)

##---------------------------------------------------------
## Fig 1c
##---------------------------------------------------------

rm(list=ls())

##-----------------------------
## Age
##-----------------------------

source("scripts/source_col.R")
load("data/processed_data_nhs.RData"); age_nhs_pro <- model_data %>% dplyr::select(age=age_final) %>% mutate(cohort="NHS (prospective)")
load("data/processed_data_nhs_caco.RData"); age_nhs_caco <- model_data %>% dplyr::select(age=age_final) %>% mutate(cohort="NHS (case-control)")
load("data/processed_data_dp.RData"); age_dp <- merged_data %>% distinct(sno, .keep_all = TRUE) %>% mutate(cohort="DIRECT-PLUS") %>% dplyr::select(cohort, age=age_baseline)

age_dat <- as.data.frame(rbind(age_nhs_pro, age_nhs_caco, age_dp)) %>%
  mutate(cohort = factor(cohort, levels = c("NHS (prospective)", "NHS (case-control)", "DIRECT-PLUS"))) %>%
  filter(!is.na(age))
table(age_dat$cohort)

p_age <- ggplot(age_dat, aes(x = age, y = cohort, fill = cohort)) +
  geom_density_ridges(
    scale = 1.3, rel_min_height = 0.01,
    color = "black", linewidth = 0.2, alpha = 0.85
  ) +
  scale_fill_manual(values = study_col) +
  scale_x_continuous(name = "Age at baseline (years)", expand = c(0, 0)) +
  scale_y_discrete(name = "", limits = rev) +
  guides(fill = "none") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.25),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(3, "pt"),
    axis.text.x = element_text(color = "black", size = 9),
    axis.title.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.title  = element_text(color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none",
    legend.text = element_text(colour = "black", size = 10),
    legend.key.size = unit(3, "mm")
  ) 
  
ggsave("results/figures/fig1c_age.pdf", p_age, width = 2.5, height = 2)

##-----------------------------
## Age and NULISA
##-----------------------------

rm(list = ls())

source("scripts/source_rm_markers.R")

age_res_dat <- read_excel("results/dp_age_NULISA_boot.xlsx") %>%
  filter(!(targetName %in% ipc_bcg_fails))

nhs_dat <- read_excel("results/nhs_age_NULISA_boot.xlsx") %>%
  filter(!(targetName %in% ipc_bcg_fails)) %>%
  mutate(targetName_anno = gsub("β", "b", targetName_orig))

caco_dat <- read_excel("results/nhs_caco_age_NULISA_boot.xlsx") %>%
  filter(!(targetName %in% ipc_bcg_fails)) %>%
  mutate(targetName_anno = gsub("β", "b", targetName_orig))

table(age_res_dat$p_age10 < 0.05)
table(age_res_dat$`beta_lcl.2.5%` * age_res_dat$`beta_ucl.97.5%` > 0)

top_dat <- age_res_dat %>%
  filter(`beta_lcl.2.5%` * `beta_ucl.97.5%` > 0) %>%
  mutate(targetName_anno = gsub("β", "b", targetName_orig))

top_dat <- top_dat %>%
  mutate(
    dir = ifelse(beta_age10 > 0, "Increase", "Decrease")
  )

nhs_dat  <- nhs_dat  %>% filter(targetName %in% top_dat$targetName)
caco_dat <- caco_dat %>% filter(targetName %in% top_dat$targetName)

### Plot

p_r2_bar <- ggplot(
  top_dat,
  aes(x = reorder(targetName_anno, r2), y = r2, fill = dir)
) +
  geom_col(
    width = 0.7,
    alpha = 0.9
  ) +
  geom_errorbar(
    aes(ymin = `r2_lcl.2.5%`, ymax = `r2_ucl.97.5%`),
    width = 0.4,
    linewidth = 0.4,
    color = "black"
  ) +
  geom_point(size = 0.5, color = "black", legend=F) +
  scale_fill_manual(
    values = c(
      "Increase" = "#B24C4C",
      "Decrease" = "#6C83B5",
      "NHS | Prospective cohort" = "white",
      "NHS | Nested case-control" = "black"
    ),
    breaks = c("Increase", "Decrease", "NHS | Prospective cohort", "NHS | Nested case-control"),
    name = "Association with age \nin DIRECT-PLUS"
  )+
  ## NHS prospective 
  geom_point(
    data = nhs_dat %>% mutate(cohort = "NHS | Prospective cohort"),
    aes(x = targetName_anno, y = r2, fill = cohort),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.6,
    color = "black"
  ) +
  
  ## NHS caco 
  geom_point(
    data = caco_dat %>% mutate(cohort = "NHS | Nested case-control"),
    aes(x = targetName_anno, y = r2, fill = cohort),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.6,
    color = "black"
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Variance explained by age (" * italic(R)^2 * ")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 12.5, color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13.5, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.title = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.position = "top"
  )

ggsave(
  "results/figures/fig1c_age_nhs.pdf",
  p_r2_bar,
  width = 5,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

### Source data

target_order_fig1c <- top_dat %>%
  arrange(r2) %>%
  transmute(
    targetName,
    targetName_anno
  )

source_data_fig1c_age_nulisa <- bind_rows(
  top_dat %>%
    left_join(target_order_fig1c, by = c("targetName", "targetName_anno")) %>%
    transmute(
      target_name_anno = targetName_anno,
      cohort = "DIRECT-PLUS",
      association_with_age = dir,
      r2 = r2,
      r2_lcl = `r2_lcl.2.5%`,
      r2_ucl = `r2_ucl.97.5%`
    ),
  
  nhs_dat %>%
    left_join(target_order_fig1c, by = c("targetName", "targetName_anno")) %>%
    transmute(
      target_name_anno = targetName_anno,
      cohort = "NHS | Prospective cohort",
      association_with_age = NA_character_,
      r2 = r2,
      r2_lcl = NA_real_,
      r2_ucl = NA_real_
    ),
  
  caco_dat %>%
    left_join(target_order_fig1c, by = c("targetName", "targetName_anno")) %>%
    transmute(
      target_name_anno = targetName_anno,
      cohort = "NHS | Nested case-control",
      association_with_age = NA_character_,
      r2 = r2,
      r2_lcl = NA_real_,
      r2_ucl = NA_real_
    )
) %>%
  arrange(cohort)

source_data_fig1c_age_distribution <- age_dat %>%
  transmute(
    cohort = as.character(cohort),
    age_years = age
  ) %>%
  arrange(cohort, age_years)

write_xlsx(
  list(
    source_data_fig1c_age_distribution = source_data_fig1c_age_distribution,
    source_data_fig1c_age_nulisa = source_data_fig1c_age_nulisa
  ),
  "results/source_data/source_data_fig1c.xlsx"
)

##---------------------------------------------------------
## Fig 1d
##---------------------------------------------------------

rm(list = ls())

load("data/processed_data_dp.RData")
marker_dict <- read_excel("data/NULISA_dictionary.xlsx")

baseline_data <- npq_all %>% 
  filter(time == "A") %>% 
  left_join(allgeno, by = "sno", suffix = c("", "_genetic")) %>% 
  filter(
    !is.na(group),
    !is.na(age_baseline)
  )

age_res_dat <- read_excel("results/dp_age_NULISA_boot.xlsx")

markers_to_plot <- age_res_dat %>% 
  arrange(-r2) %>% 
  slice_head(n = 5) %>% 
  pull(targetName)
markers_to_plot <- c(markers_to_plot, "pTau_217")

plot_long <- baseline_data %>%
  dplyr::select(sno, age_baseline, all_of(markers_to_plot)) %>%
  pivot_longer(
    cols = all_of(markers_to_plot),
    names_to = "marker",
    values_to = "value"
  ) %>%
  group_by(marker) %>%
  mutate(
    value_z = as.numeric(scale(value))
  ) %>%
  ungroup() %>%
  filter(!is.na(value), !is.na(age_baseline)) %>%
  left_join(
    marker_dict %>% dplyr::select(targetName, targetName_anno),
    by = c("marker" = "targetName")
  ) %>%
  mutate(
    marker_lab = case_when(
      marker == "Ab42_div_Ab40_log2" ~ "Aβ42/Aβ40",
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

x_grid <- seq(
  min(plot_long$age_baseline, na.rm = TRUE),
  max(plot_long$age_baseline, na.rm = TRUE),
  length.out = 200
)

fit_one_marker <- function(dat) {
  fit <- scam(
    value_z ~ s(age_baseline, k = 8),
    data = dat
  )
  
  newdat <- data.frame(age_baseline = x_grid)
  pred <- predict(fit, newdata = newdat, se.fit = TRUE)
  
  tibble(
    age_baseline = x_grid,
    fit = pred$fit,
    lwr = pred$fit - 1.96 * pred$se.fit,
    upr = pred$fit + 1.96 * pred$se.fit
  )
}

pred_df <- plot_long %>%
  group_by(marker_lab) %>%
  group_modify(~ fit_one_marker(.x)) %>%
  ungroup()

turn_df <- pred_df %>%
  filter(marker_lab %in% c("NEFL", "GFAP")) %>%
  arrange(marker_lab, age_baseline) %>%
  group_by(marker_lab) %>%
  mutate(
    slope = c(NA, diff(fit) / diff(age_baseline))
  ) %>%
  filter(!is.na(slope)) %>%
  summarise(
    turn_age_max = age_baseline[which.max(slope)],
    turn_y_max   = fit[which.max(slope)],
    turn_age_min = age_baseline[which.min(slope)],
    turn_y_min   = fit[which.min(slope)],
    .groups = "drop"
  )
turn_df_long <- turn_df %>%
  tidyr::pivot_longer(
    cols = -marker_lab,
    names_to = c(".value", "type"),
    names_pattern = "turn_(.*)_(.*)"
  )

curve_cols <- c(
  "NEFL"     = "#2C7FB8",  
  "GFAP"     = "#B44B4B", 
  "CHI3L1"   = "#7B6FD0", 
  "CCL11"    = "#D28A2E", 
  "GDF15"    = "#1F8A70",  
  "pTau-217" = "#7A9A3A"   
  
)

curve_cols <- curve_cols[names(curve_cols) %in% levels(plot_long$marker_lab)]

### Plot

p_agetraj <- ggplot() +
  geom_line(
    data = pred_df,
    aes(x = age_baseline, y = fit, color = marker_lab),
    linewidth = 0.3,
    alpha = 0.85
  ) +
  geom_segment(
    data = turn_df_long,
    aes(
      x = age,
      xend = age,
      y = y,
      yend = min(pred_df$lwr, pred_df$fit, na.rm = TRUE),
      color = marker_lab
    ),
    linetype = "42",
    linewidth = 0.15,
    show.legend = FALSE
  ) +
  geom_point(
    data = turn_df_long,
    aes(x = age, y = y, color = marker_lab),
    size = 0.4,
    show.legend = FALSE
  ) +
  geom_rug(
    data = plot_long %>% distinct(sno, age_baseline),
    aes(x = age_baseline),
    sides = "b",
    inherit.aes = FALSE,
    alpha = 0.2,
    linewidth = 0.20
  ) +
  scale_color_manual(
    values = curve_cols,
    breaks = names(curve_cols)
  ) +
  scale_fill_manual(
    values = curve_cols
  ) +
  scale_x_continuous(
    breaks = pretty(baseline_data$age_baseline, n = 5)
  ) +
  coord_cartesian(ylim = c(-1.2, NA)) +
  labs(
    x = "Age (years)",
    y = "NULISA Protein\nQuantification (z-score)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black"),
    axis.line  = element_line(linewidth = 0.25, color = "black"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.key.width = grid::unit(0.7, "cm"),
    plot.margin = margin(8, 8, 8, 8),
    panel.grid = element_blank()
  )

ggsave(
  "results/figures/fig1d_age_traj.pdf",
  p_agetraj,
  width = 5.8,
  height = 1.8
)

### Source data

source_data_fig1d_curves <- pred_df %>%
  transmute(
    marker = as.character(marker_lab),
    age_years = age_baseline,
    fitted_z_score = fit
  ) %>%
  arrange(marker, age_years)

source_data_fig1d_turnpoints <- turn_df_long %>%
  transmute(
    marker = as.character(marker_lab),
    turnpoint_type = type,
    age_years = age,
    fitted_z_score = y
  ) %>%
  arrange(marker, turnpoint_type)

write_xlsx(
  list(
    source_data_fig1d_curves = source_data_fig1d_curves,
    source_data_fig1d_turnpoints = source_data_fig1d_turnpoints
  ),
  "results/source_data/source_data_fig1d.xlsx"
)

##---------------------------------------------------------
## Fig 1e 
##---------------------------------------------------------

rm(list = ls())

load("data/processed_data_nhs.RData")
load("data/NULISA_core.RData")
source("scripts/source_col.R")
source("scripts/source_dietcomp.R")

protall_list <- prot_list[
  !grepl("raw", prot_list) &
    !(prot_list %in% c(marker_dict %>%
                        filter(LOD_flag | ipc_bcg_flag) %>%
                        pull(targetName)))
]

diet_vars_use <- nhs_dietcomp_dict %>%
  filter(score %in% c("MIND", "aMED", "MIND component")) %>%
  pull(diet) %>%
  unique()

prs_var_use <- c("PGS002280_scaled_pool_rsq0.3", "PGS000334_scaled_pool_rsq0.3")

pca_dat <- npq_all %>%
  dplyr::select(id, all_of(protall_list)) %>%
  left_join(
    model_data %>%
      dplyr::select(id, APOE4_carrier),
    by = "id"
  ) %>%
  filter(
    if_all(all_of(protall_list), ~ !is.na(.)),
    !is.na(APOE4_carrier)
  ) %>%
  mutate(
    APOE4_anno = factor(
      ifelse(APOE4_carrier == 1, "Carrier", "Noncarrier"),
      levels = c("Carrier", "Noncarrier")
    )
  )

X_pca <- pca_dat %>%
  dplyr::select(all_of(protall_list))

pca_fit <- prcomp(X_pca, center = TRUE, scale. = TRUE)

diet_dat <- model_data %>%
  dplyr::select(id, APOE4_carrier, any_of(prs_var_use), any_of(diet_vars_use)) %>%
  distinct()

pc_scores <- as.data.frame(predict(pca_fit, newdata = X_pca)) %>%
  bind_cols(
    pca_dat %>% dplyr::select(id, APOE4_carrier, APOE4_anno)
  ) %>%
  left_join(diet_dat, by = c("id", "APOE4_carrier"))

pc_scores <- pc_scores %>%
  mutate(across(any_of(c(prs_var_use, diet_vars_use)), ~ suppressWarnings(as.numeric(.))))

pc_var <- data.frame(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  var_explained = (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
)

arrow_vars <- nhs_dietcomp_dict %>%
  filter(score %in% c("MIND", "aMED", "MIND component")) %>%
  pull(diet) %>%
  unique()

pco_corr_list <- list()
i <- 1

for (var in arrow_vars) {
  
  if (!var %in% colnames(pc_scores)) next
  
  x  <- pc_scores[[var]]
  y1 <- pc_scores$PC1
  y2 <- pc_scores$PC2
  
  if (sum(!is.na(x)) < 10) next
  if (length(unique(x[!is.na(x)])) < 2) next
  
  tmp_cor_pc1 <- tryCatch(
    cor.test(x, y1, use = "complete.obs"),
    error = function(e) NULL
  )
  
  tmp_cor_pc2 <- tryCatch(
    cor.test(x, y2, use = "complete.obs"),
    error = function(e) NULL
  )
  
  if (is.null(tmp_cor_pc1) || is.null(tmp_cor_pc2)) next
  
  pco_corr_list[[i]] <- c(
    var,
    sum(!is.na(x)),
    tmp_cor_pc1$estimate, tmp_cor_pc1$p.value,
    tmp_cor_pc2$estimate, tmp_cor_pc2$p.value
  )
  
  i <- i + 1
}

pco_corr_dat <- as.data.frame(do.call("rbind", pco_corr_list))
colnames(pco_corr_dat) <- c("diet", "n", "corr_pc1", "p_pc1", "corr_pc2", "p_pc2")

range_pc1 <- diff(range(pc_scores$PC1, na.rm = TRUE))
range_pc2 <- diff(range(pc_scores$PC2, na.rm = TRUE))
scale_factor <- 1.5 * max(range_pc1, range_pc2)

pco_corr_dat <- pco_corr_dat %>%
  mutate(
    across(c(n, corr_pc1, p_pc1, corr_pc2, p_pc2), as.numeric)
  ) %>%
  left_join(nhs_dietcomp_dict, by = "diet") %>%
  mutate(
    arrow_x = corr_pc1 * scale_factor,
    arrow_y = corr_pc2 * scale_factor,
    arrow_col = diet_col
  ) %>%
  filter(score %in% c("MIND", "aMED", "MIND component"))

##-----------------------------
## PCA plot
##-----------------------------

p_pc_apoe <- ggplot(pc_scores, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(
      fill  = PGS002280_scaled_pool_rsq0.3,
      shape = APOE4_anno
    ),
    size = 3,
    stroke = 0.1,
    color = "gray31",
    alpha = 0.85
  ) +
  scale_shape_manual(
    values = c(Carrier = 22, Noncarrier = 21),
    name = "APOE4 status"
  ) +
  scale_fill_gradient2(
    low      = "#7c9559",
    mid      = "white",
    high     = "#4a75a0",
    midpoint = 0,
    name     = "ADRD PRS",
    guide = guide_colorbar(
      direction = "horizontal",
      barwidth  = unit(35, "mm"),
      barheight = unit(3, "mm")
    )
  ) +
  xlab(paste0("Proteomic PC1 (", round(pc_var$var_explained[1] * 100, 1), "%)")) +
  ylab(paste0("Proteomic PC2 (", round(pc_var$var_explained[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 13, colour = "black"),
    axis.text.y        = element_text(size = 13, colour = "black"),
    axis.title.x       = element_text(size = 14, colour = "black"),
    axis.title.y       = element_text(size = 14, colour = "black"),
    axis.ticks         = element_line(colour = "black", linewidth = 0.5),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.title         = element_text(size = 14, colour = "black", hjust = 0.5),
    legend.title       = element_text(size = 11, colour = "black"),
    legend.text        = element_text(size = 10, colour = "black"),
    legend.position    = c(0.66, 0.8),
    legend.background  = element_blank(),
    legend.box         = "vertical"
  ) +
  ggnewscale::new_scale_color() +
  geom_segment(
    data = pco_corr_dat,
    aes(
      x = 0, y = 0,
      xend = arrow_x,
      yend = arrow_y,
      color = arrow_col
    ),
    arrow = arrow(length = unit(0.08, "inches")),
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  scale_color_identity(guide = "none") +
  ggnewscale::new_scale_fill() +
  ggrepel::geom_label_repel(
    data = pco_corr_dat,
    aes(
      x = arrow_x,
      y = arrow_y,
      label = diet_anno,
      color = arrow_col,
      segment.color = arrow_col
    ),
    fill = scales::alpha("white", 0.8), 
    size             = 4.5,
    fontface         = "plain",
    box.padding      = 0.28,
    point.padding    = 0.18,
    label.size       = 0.3,
    label.r          = unit(0.10, "lines"),
    segment.size     = 0.35,
    segment.alpha    = 0.85,
    max.overlaps     = Inf,
    force            = 2.2,
    show.legend      = FALSE
  ) +
  scale_color_identity(guide = "none") +
  scale_fill_identity() +
  scale_x_continuous(
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = expression("Correlation " * italic(r) * " with proteomic PC1")
    ),
    expand = expansion(mult = 0.10)
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = expression("Correlation " * italic(r) * " with proteomic PC2")
    ),
    expand = expansion(mult = 0.10)
  )

ggsave(
  "results/figures/fig1_pca_apoe.pdf",
  p_pc_apoe,
  width = 5,
  height = 5,
  dpi = 600,
  bg = "white"
)

### Source data

source_data_fig1e_points <- pc_scores %>%
  transmute(
    PC1 = PC1,
    PC2 = PC2,
    APOE4_status = as.character(APOE4_anno),
    ADRD_PRS = PGS002280_scaled_pool_rsq0.3
  ) %>%
  arrange(APOE4_status, PC1, PC2)

source_data_fig1e_diet_arrows <- pco_corr_dat %>%
  transmute(
    diet_label = diet_anno,
    correlation_with_PC1 = corr_pc1,
    correlation_with_PC2 = corr_pc2,
    plotted_arrow_x = arrow_x,
    plotted_arrow_y = arrow_y
  ) %>%
  arrange(diet_label)

write_xlsx(
  list(
    source_data_fig1e_points = source_data_fig1e_points,
    source_data_fig1e_diet_arrows = source_data_fig1e_diet_arrows
  ),
  "results/source_data/source_data_fig1e.xlsx"
)
