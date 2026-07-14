library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(readxl)
library(vegan)
library(ggnewscale)
library(purrr)
library(writexl)
library(ggrepel)
library(patchwork)

##---------------------------------------------------------
## Fig 3a 
##---------------------------------------------------------

rm(list=ls())

load("data/processed_data_nhs_2nd_calibrated.RData")

marker_dict %>% filter(LOD_flag)

summary(bld_data$bld_yr_diff)

protall_list <- prot_list[
  !grepl("raw", prot_list) &
    !(prot_list %in% (marker_dict %>% filter(LOD_flag | ipc_bcg_flag) %>% pull(targetName)))
]

fig_dat <- npq_all %>% 
  select(id, 
         all_of(paste0(protall_list,"_1")), 
         all_of(paste0(protall_list,"_2")))

X1 <- fig_dat %>% select(all_of(paste0(protall_list, "_1")))
X2 <- fig_dat %>% select(all_of(paste0(protall_list, "_2")))

colnames(X1) <- protall_list
colnames(X2) <- protall_list

X_pca <- bind_rows(X1, X2)
stopifnot(!anyNA(X_pca))

pca_fit <- prcomp(X_pca, center = TRUE, scale. = TRUE)

pc_1 <- as.data.frame(predict(pca_fit, newdata = X1)) %>%
  mutate(id = fig_dat$id, time_anno = "1st blood")

pc_2 <- as.data.frame(predict(pca_fit, newdata = X2)) %>%
  mutate(id = fig_dat$id, time_anno = "2nd blood")

pc_long <- bind_rows(pc_1, pc_2) %>%
  mutate(time_anno = factor(time_anno, levels = c("1st blood", "2nd blood")))

pc_var <- data.frame(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  var_explained = (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
)

pc_loadings <- as.data.frame(pca_fit$rotation) %>%
  mutate(marker = rownames(.))

delta_mat <- as.matrix(X1) - as.matrix(X2)

colnames(delta_mat) <- protall_list

dist_mat <- dist(X_pca, method = "euclidean")

meta <- data.frame(
  id = rep(fig_dat$id, 2),
  time   = rep(c("1st blood","2nd blood"), each = nrow(fig_dat))
)

adonis_res <- adonis2(
  dist_mat ~ time,
  data = meta,
  permutations = 9999,
  strata = meta$id
)

permanova_label <- paste0(
  "PERMANOVA ",
  "R² = ", sprintf("%.1f", adonis_res$R2[1] * 100), "%\n",
  "P = ", sprintf("%.5f", adonis_res$`Pr(>F)`[1])
)

### Plot

p_pc_2nd <- ggplot(pc_long, aes(x = PC1, y = PC2)) +
  
  geom_line(
    aes(group = id),
    color = "gray31",
    alpha = 0.4,
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  
  geom_point(
    aes(fill = time_anno),
    shape = 21,
    size = 3.4,
    stroke = 0.4,
    color = "black",
    alpha = 0.8
  ) +
  
  scale_fill_manual(
    values = c(`1st blood` = "#6FAED6", `2nd blood` = "#D8A45B"),
    name = ""
  ) +
  
  xlab(paste0("Proteomic PC1 (", round(pc_var$var_explained[1] * 100, 1), "%)")) +
  ylab(paste0("Proteomic PC2 (", round(pc_var$var_explained[2] * 100, 1), "%)")) +
  ggtitle("Global proteomic changes in NHS") +
  
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
  stat_ellipse(
    aes(fill = time_anno),
    geom = "polygon",
    type = "norm",
    alpha = 0.10,
    color = NA,
    show.legend = FALSE
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

ggsave("results/figures/fig3_pca_nhs_2nd.pdf", p_pc_2nd, width = 5, height = 5)

### Source data

source_data_fig3a_points <- pc_long %>%
  mutate(pair_id = as.integer(factor(id))) %>%
  transmute(
    pair_id,
    time = as.character(time_anno),
    PC1,
    PC2
  ) %>%
  arrange(pair_id, time)

write_xlsx(source_data_fig3a_points, "results/source_data/source_data_fig3a.xlsx")

##---------------------------------------------------------
## Fig 3b 
##---------------------------------------------------------

rm(list=ls())

load("data/processed_data_nhs_2nd_calibrated.RData")
source("scripts/source_col.R")

marker_dict_anno <- read_excel("data/NULISA_dictionary.xlsx") %>% 
  select(targetName, hallmark_anno)

fig_dat <- read_excel("results/nhs_2nd_NULISA_paired_tests.xlsx") %>% 
  left_join(marker_dict, by = c("marker" = "targetName")) %>% 
  left_join(marker_dict_anno, by = c("marker" = "targetName")) %>% 
  filter(
    APOE4_group == "all",
    LOD_flag != TRUE,
    ipc_bcg_flag != TRUE,
    p_w < 0.05
  ) %>% 
  mutate(
    hallmark_anno = factor(hallmark_anno, levels = names(hallmark_col)),
    targetName_anno = if_else(grepl("_", marker), gsub("_", "-", marker), marker),
    change_abs = abs(median_change_yr),
    change_dir = if_else(median_change_yr >= 0, "up", "down"),
    signed_logp = sign(median_change_yr) * -log10(p_w)
  ) %>%
  arrange(hallmark_anno, targetName_anno) %>%
  mutate(
    targetName_anno = factor(targetName_anno, levels = targetName_anno)
  )

cap_val <- quantile(fig_dat$change_abs, 0.95, na.rm = TRUE)

fig_dat <- fig_dat %>%
  mutate(
    change_abs_plot = pmin(change_abs, cap_val)
  )

max_y <- max(fig_dat$change_abs_plot, na.rm = TRUE)

ring_center <- -0.03 * max_y
ring_height <- 0.015 * max_y

inner_hole_limit <- -0.2 * max_y

axis_gap  <- 4.5
bar_width <- 0.90

fig_plot <- fig_dat %>%
  mutate(
    x_raw = as.numeric(targetName_anno),
    x_num = x_raw + axis_gap
  )

n_bar <- max(fig_plot$x_raw)
max_x <- n_bar + axis_gap

x_axis <- 0.5 + axis_gap / 2

start_angle <- -pi/2 - 2*pi*(x_axis - 0.5) / max_x

ring_dat <- fig_plot %>%
  mutate(ring_y = ring_center)

tick_vals <- c(0, 0.33 * max_y, 0.66 * max_y, max_y)
tick_labs <- c(
  "0",
  sprintf("%.2f", 0.33 * max_y),
  sprintf("%.2f", 0.66 * max_y),
  sprintf("%.2f", max_y)
)

tick_df <- tibble(
  x = x_axis,
  y = tick_vals,
  lab = tick_labs
)

label_df <- fig_plot %>%
  mutate(
    angle_raw = 90 - 360 * (x_num - 0.5) / max_x - start_angle * 180 / pi,
    angle = ifelse(angle_raw < -90, angle_raw + 180, angle_raw),
    hjust = ifelse(angle_raw < -90, 1, 0)
  )

label_gap <- 0.02 * max_y

### Plot

p_change <- ggplot(fig_plot, aes(x = x_num, y = change_abs_plot)) +
  
  geom_hline(
    yintercept = tick_vals,
    color = "grey88",
    linewidth = 0.3
  ) +
  
  geom_col(
    aes(fill = signed_logp),
    width = bar_width,
    color = NA,
    linewidth = 0.2
  ) +
  
  scale_fill_gradient2(
    low  = "#7A97C9",
    mid  = "#F5F5F5",
    high = "#D77A61",
    midpoint = 0,
    name = expression("sign(change) × " * -log[10](P))
    ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    data = ring_dat,
    aes(
      x = x_num,
      y = ring_y,
      fill = hallmark_anno
    ),
    width = 1.0,
    height = ring_height,
    color = NA,
    alpha = 0.75,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = hallmark_col,
    name = "Hallmarks"
  ) +
  
  geom_text(
    data = label_df,
    aes(
      x = x_num,
      y = change_abs_plot + label_gap, 
      label = targetName_anno,
      angle = angle,
      hjust = hjust
    ),
    size = 2.5,             
    color = "gray31",
    fontface = "plain",     
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = tick_df,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.5,
    color = "grey35"
  ) +
  
  coord_polar(theta = "x", start = start_angle, clip = "off") +
  
  scale_x_continuous(
    limits = c(0.5, max_x + 0.5),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(inner_hole_limit, 1.15 * max_y),
    expand = c(0, 0)
  ) +
  
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(), 
    plot.margin = margin(10, 30, 10, 30)
  )

ggsave(
  "results/figures/fig3b_change_sig.pdf",
  p_change,
  width = 15,
  height = 15
)

### Source data

source_data_fig3b <- fig_plot %>%
  transmute(
    target_name_anno = as.character(targetName_anno),
    hallmark_anno = as.character(hallmark_anno),
    annualized_median_change = median_change_yr,
    plotted_absolute_annualized_median_change = change_abs_plot,
    wilcoxon_p_value = p_w
  ) 

write_xlsx(
  list(source_data_fig3b = source_data_fig3b),
  "results/source_data/source_data_fig3b.xlsx"
)

##---------------------------------------------------------
## Fig 3c 
##---------------------------------------------------------

rm(list=ls())

load("data/processed_data_nhs_2nd_calibrated.RData")
source("scripts/source_col.R")
source("scripts/source_rm_markers.R")

res_paired <- read_excel("results/nhs_2nd_NULISA_paired_tests.xlsx") %>% 
  left_join(marker_dict, by=c("marker"="targetName")) %>% 
  filter(APOE4_group == "all",
         LOD_flag == F,
         ipc_bcg_flag == F,
         p_w < 0.05)

gene_res_dat <- read_excel("results/nhs_2nd_NULISA_gene.xlsx")
gene_res_dat <- gene_res_dat %>%
  filter(p < 0.05,
         LOD_flag == F,
         !(targetName %in% ipc_bcg_fails),
         prot_type == "chg_yr")

### Overall

marker_show <- c("pTau_217", "pTau_181", "GFAP", "CRH", "TREM2")

marker_lab <- c(
  pTau_217 = "pTau-217",
  pTau_181 = "pTau-181",
  GFAP     = "GFAP",
  CRH      = "CRH",
  TREM2    = "TREM2"
)

marker_levels_overall <- c("pTau-217", "pTau-181", "GFAP", "CRH", "TREM2")

time_col <- c(
  `1st blood` = "#6FAED6",
  `2nd blood` = "#D8A45B"
)

fig_dat <- npq_all %>%
  select(id, all_of(paste0(marker_show, "_1")), all_of(paste0(marker_show, "_2"))) %>%
  pivot_longer(
    cols = -id,
    names_to = c("marker", "time_code"),
    names_pattern = "^(.*)_(1|2)$",
    values_to = "value"
  ) %>%
  mutate(
    time_anno = case_when(
      time_code == "1" ~ "1st blood",
      time_code == "2" ~ "2nd blood"
    ),
    time_anno = factor(time_anno, levels = c("1st blood", "2nd blood")),
    marker = dplyr::recode(marker, !!!marker_lab),
    marker = factor(marker, levels = marker_levels_overall)
  ) %>%
  filter(!is.na(value))

stat_df <- map_dfr(marker_show, function(m) {
  x1 <- npq_all[[paste0(m, "_1")]]
  x2 <- npq_all[[paste0(m, "_2")]]
  keep <- !is.na(x1) & !is.na(x2)
  wt <- wilcox.test(x2[keep], x1[keep], paired = TRUE, exact = FALSE)
  
  tibble(
    marker_raw = m,
    marker = marker_lab[[m]],
    n = sum(keep),
    p_w = wt$p.value,
    p_lab = paste0("P = ", formatC(wt$p.value, format = "e", digits = 2))
  )
}) %>%
  mutate(marker = factor(marker, levels = marker_levels_overall))

y_pos <- fig_dat %>%
  group_by(marker) %>%
  dplyr::summarise(
    y = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  )

stat_df <- stat_df %>%
  left_join(y_pos, by = "marker")

### Plot

p_marker_change <- ggplot(fig_dat, aes(x = time_anno, y = value)) +
  geom_line(
    aes(group = id),
    color = "grey70",
    alpha = 0.6,
    linewidth = 0.2,
    na.rm = TRUE
  ) +
  geom_boxplot(
    aes(fill = time_anno),
    width = 0.42,
    outlier.shape = NA,
    color = "grey35",
    linewidth = 0.35,
    alpha = 0.35
  ) +
  geom_point(
    aes(fill = time_anno),
    shape = 21,
    size = 1.35,
    stroke = 0.22,
    color = "grey25",
    alpha = 0.70,
    position = position_jitter(width = 0.04, height = 0, seed = 1),
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.45,
    color = "black",
    fatten = 0
  ) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    linewidth = 0.75,
    color = "grey31"
  ) +
  stat_summary(
    aes(group = 1, fill = time_anno),
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.0,
    stroke = 0.7,
    color = "grey31"
  ) +
  geom_text(
    data = stat_df,
    aes(x = 1.5, y = y, label = p_lab),
    inherit.aes = FALSE,
    size = 3.0,
    color = "black"
  ) +
  facet_wrap(~ marker, nrow = 2, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = time_col, name = NULL) +
  labs(
    x = NULL,
    y = "NULISA Protein Quantification"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, color = "black"),
    legend.position = "none",
    panel.spacing = unit(0.9, "lines")
  )

### Genetic callout

callout_map <- tibble::tribble(
  ~marker_raw,    ~marker_plot,   ~prs_var,                         ~prs_label,
  "BD_pTau_231",  "BD-pTau-231",  "PGS000334_scaled_pool_rsq0.3",   "ADRD PRS (with APOE4)"
)

marker_levels_genetic <- callout_map$marker_plot

prs_col <- c(
  "Low"    = "#9FA8A3",
  "Medium" = "#D2A876",
  "High"   = "#B65C5C"
)

fig_dat_g <- purrr::map_dfr(seq_len(nrow(callout_map)), function(i) {
  
  marker_i  <- callout_map$marker_raw[i]
  prs_i     <- callout_map$prs_var[i]
  label_i   <- callout_map$marker_plot[i]
  prs_lab_i <- callout_map$prs_label[i]
  
  model_data %>%
    dplyr::select(
      id,
      all_of(prs_i),
      all_of(paste0(marker_i, "_1")),
      all_of(paste0(marker_i, "_2"))
    ) %>%
    dplyr::mutate(
      PRS_tertile_num = dplyr::ntile(.data[[prs_i]], 3),
      PRS_tertile = factor(
        PRS_tertile_num,
        levels = c(1, 2, 3),
        labels = c("Low", "Medium", "High")
      )
    ) %>%
    pivot_longer(
      cols = all_of(c(paste0(marker_i, "_1"), paste0(marker_i, "_2"))),
      names_to = c("marker_raw_tmp", "time_code"),
      names_pattern = "^(.*)_(1|2)$",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      marker = label_i,
      prs_var = prs_i,
      prs_label = prs_lab_i,
      time_anno = case_when(
        time_code == "1" ~ "1st blood",
        time_code == "2" ~ "2nd blood"
      ),
      time_anno = factor(time_anno, levels = c("1st blood", "2nd blood")),
      PRS_tertile = factor(PRS_tertile, levels = c("Low", "Medium", "High"))
    ) %>%
    dplyr::select(id, marker, prs_var, prs_label, PRS_tertile, time_anno, value)
}) %>%
  dplyr::filter(!is.na(value), !is.na(PRS_tertile)) %>%
  dplyr::mutate(marker = factor(marker, levels = marker_levels_genetic))

stat_df_g <- callout_map %>%
  left_join(
    gene_res_dat,
    by = c("marker_raw" = "targetName", "prs_var" = "gene")
  ) %>%
  mutate(
    marker = factor(marker_plot, levels = marker_levels_genetic),
    p_lab = paste0("P = ", formatC(p, format = "e", digits = 2))
  ) %>%
  select(marker, marker_raw, prs_var, prs_label, p, p_lab)

y_pos_g <- fig_dat_g %>%
  group_by(marker) %>%
  dplyr::summarise(
    y = max(value, na.rm = TRUE) + 0.12 * diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  )

stat_df_g <- stat_df_g %>%
  left_join(y_pos_g, by = "marker") %>%
  mutate(x = 1.5)

offset_map <- c(Low = -0.2, Medium = 0, High = 0.2)

fig_dat_g2 <- fig_dat_g %>%
  mutate(
    time_num = ifelse(time_anno == "1st blood", 1, 2),
    x_plot = time_num + offset_map[as.character(PRS_tertile)]
  )

### Plot

p_callout <- ggplot() +
  geom_boxplot(
    data = fig_dat_g2,
    aes(x = x_plot, y = value, group = interaction(time_anno, PRS_tertile), fill = PRS_tertile),
    width = 0.16,
    outlier.shape = NA,
    color = "grey35",
    linewidth = 0.3,
    alpha = 0.35
  ) +
  geom_point(
    data = fig_dat_g2,
    aes(x = x_plot, y = value, fill = PRS_tertile),
    shape = 21,
    size = 1.25,
    stroke = 0.2,
    color = "grey25",
    alpha = 0.65,
    position = position_jitter(width = 0.02, height = 0, seed = 1),
    show.legend = FALSE
  ) +
  stat_summary(
    data = fig_dat_g2,
    aes(x = x_plot, y = value, group = interaction(time_anno, PRS_tertile)),
    fun = median,
    geom = "crossbar",
    width = 0.16,
    linewidth = 0.4,
    color = "black",
    fatten = 0
  ) +
  stat_summary(
    data = fig_dat_g2,
    aes(x = x_plot, y = value, group = PRS_tertile, color = PRS_tertile),
    fun = median,
    geom = "line",
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  stat_summary(
    data = fig_dat_g2,
    aes(x = x_plot, y = value, group = interaction(time_anno, PRS_tertile), fill = PRS_tertile),
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.1,
    stroke = 0.5,
    color = "grey30",
    show.legend = FALSE
  ) +
  geom_text(
    data = stat_df_g,
    aes(x = x, y = y, label = p_lab),
    size = 3.0,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ marker, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("1st blood", "2nd blood")
  ) +
  scale_fill_manual(values = prs_col, name = "ADRD PRS tertile") +
  scale_color_manual(values = prs_col, guide = "none") +
  labs(
    x = NULL,
    y = "NULISA Protein Quantification"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9.5),
    panel.spacing = unit(0.9, "lines")
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.7))
  )

### Combine

p_combined <- p_marker_change + p_callout +
  patchwork::plot_layout(widths = c(2, 0.85))

ggsave(
  "results/figures/fig3c_marker_change_and_prs.pdf",
  p_combined,
  width = 6.5,
  height = 5.5
)

### Source data

source_data_fig3c_marker_change_points <- fig_dat %>%
  mutate(pair_id = as.integer(factor(id))) %>%
  transmute(
    pair_id,
    marker = as.character(marker),
    time = as.character(time_anno),
    NULISA_protein_quantification = value
  ) %>%
  arrange(marker, pair_id, time)

source_data_fig3c_prs_points <- fig_dat_g2 %>%
  mutate(pair_id = as.integer(factor(id))) %>%
  transmute(
    pair_id,
    marker = as.character(marker),
    PRS_tertile = as.character(PRS_tertile),
    time = as.character(time_anno),
    x_position = x_plot,
    NULISA_protein_quantification = value
  ) %>%
  arrange(marker, PRS_tertile, pair_id, time)

write_xlsx(
  list(
    fig3c_marker_points = source_data_fig3c_marker_change_points,
    fig3c_prs_points    = source_data_fig3c_prs_points
  ),
  "results/source_data/source_data_fig3c.xlsx"
)

##---------------------------------------------------------
## Fig 3d 
##---------------------------------------------------------

rm(list=ls())

load("data/processed_data_nhs_2nd_calibrated.RData")
load("data/NULISA_core.RData")
source("scripts/source_col.R")
source("scripts/source_rm_markers.R")

manual_label_markers <- c(
  "pTau_217", "pTau_181", "pTau_231",
  "BD_pTau_217", "BD_pTau_181", "BD_pTau_231",
  "GFAP", "APOE4",
  "IL10", "AGRN"
)

res_paired <- read_excel("results/nhs_2nd_NULISA_paired_tests.xlsx") %>% 
  filter(APOE4_group == "all") %>%
  filter(!(marker %in% ipc_bcg_fails))

dementia_res <- read_excel("results/nhs_2nd_dementia_NULISA.xlsx") %>% 
  filter(
    prot_type == "chg_yr",
    scale == "scaled",
    LOD_flag == FALSE,
    !(targetName %in% ipc_bcg_fails)
  ) %>% 
  left_join(res_paired, by = c("targetName" = "marker")) %>% 
  mutate(
    beta = as.numeric(beta),
    logp_change = -log10(p_w),
    logp_dementia = -log10(p)
  )

label_df <- dementia_res %>%
  filter(targetName %in% manual_label_markers) %>%
  mutate(
    label = if_else(grepl("_", targetName), gsub("_", "-", targetName), targetName)
  )

### Plot

p_scatter <- ggplot(dementia_res, aes(x = beta, y = logp_dementia)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "gray31", linewidth = 0.4) +
  annotate(
    "text",
    x = max(dementia_res$beta, na.rm = TRUE),
    y = -log10(0.05),
    label = "italic(P) == 0.05",
    parse = TRUE,
    color = "gray31",
    vjust = -0.5, hjust = 0.8, size = 3.2
  ) +
  geom_point(
    aes(size = logp_change, fill = median_change_yr),
    shape = 21,
    color = "gray31",
    alpha = 0.8,
    stroke = 0.3
  ) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "#F5F5F5",
    high = "#d73027",
    midpoint = 0,
    name = "Change per year"
  ) +
  scale_size_continuous(name = expression(-log[10](P[change]))) +
  geom_text_repel(
    data = label_df,
    aes(label = label),
    segment.size = 0.2,
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.15,
    min.segment.length = 0,
    seed = 1
  ) +
  labs(
    x = "log(Hazard ratio of dementia risk)",
    y = expression(-log[10](italic(P)))
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    panel.border = element_blank(),
    legend.title = element_text(size = 11, color = "black"),
    legend.text  = element_text(size = 11, color = "black"),
    legend.position = "right",
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_line(color = "black", linewidth = 0.2),
    axis.title = element_text(size = 11, color = "black")
  )

ggsave(
  "results/figures/fig3d_dementia_scatter.pdf",
  p_scatter,
  width = 5,
  height = 5
)

### Source data

source_data_fig3d <- dementia_res %>%
  transmute(
    target_name_anno = if_else(grepl("_", targetName), gsub("_", "-", targetName), targetName),
    log_hazard_ratio_dementia = beta,
    p_value_dementia = p,
    neg_log10_p_dementia = logp_dementia,
    annualized_median_change = median_change_yr,
    p_value_change = p_w,
    neg_log10_p_change = logp_change
  ) %>%
  arrange(desc(neg_log10_p_dementia), target_name_anno)

write_xlsx(
  list(source_data_fig3d = source_data_fig3d),
  "results/source_data/source_data_fig3d.xlsx"
)

##---------------------------------------------------------
## Fig 3e 
##---------------------------------------------------------

rm(list = ls())

load("data/processed_data_dp.RData")
load("data/NULISA_core.RData")

marker_dict %>% filter(LOD_flag)

protall_list <- prot_list[
  !grepl("raw", prot_list) &
    !(prot_list %in% (marker_dict %>% filter(LOD_flag | ipc_bcg_flag) %>% pull(targetName)))
]

fig_dat <- npq_all %>%
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

pc_loadings <- as.data.frame(pca_fit$rotation) %>%
  mutate(marker = rownames(.))

### PERMANOVA for overall time effect

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

p_pc_dp_AD <- ggplot(pc_long, aes(x = PC1, y = PC2)) +
  
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
  ggtitle("Global proteomic changes in DIRECT-PLUS") +
  
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
  "results/figures/fig3e_pca_dp_baseline_fit.pdf",
  p_pc_dp_AD,
  width = 5,
  height = 5
)

### Source data

source_data_fig3e_points <- pc_long %>%
  mutate(pair_id = as.integer(factor(sno))) %>%
  transmute(
    pair_id,
    time = as.character(time_anno),
    PC1,
    PC2
  ) %>%
  arrange(pair_id, time)

write_xlsx(
  list(
    source_data_fig3e_points = source_data_fig3e_points
  ),
  "results/source_data/source_data_fig3e.xlsx"
)

##---------------------------------------------------------
## Fig 3f
##---------------------------------------------------------

rm(list = ls())

load("data/processed_data_dp.RData")
source("scripts/source_rm_markers.R")

marker_dict_anno <- read_excel("data/NULISA_dictionary.xlsx") %>%
  dplyr::select(targetName, targetName_anno)

marker_qc <- marker_dict %>%
  dplyr::select(targetName, targetName_orig, LOD_flag, ipc_bcg_flag)

time_col <- c(
  "Baseline" = "#8DB7A8",  
  "FIT"      = "#D7A07A" 
)

markers_to_plot <- c(
  "pTau_217",
  "FCN2",
  "MME",
  "TARDBP",
  "VEGFA",
  "GDF15"
)

marker_label_order <- c(
  "pTau-217",
  "FCN2",
  "MME",
  "TDP43",
  "VEGFA",
  "GDF15"
)

qc_check_fig3f <- marker_qc %>%
  dplyr::filter(targetName %in% markers_to_plot) %>%
  dplyr::select(targetName, targetName_orig, LOD_flag, ipc_bcg_flag)

print(qc_check_fig3f)

stopifnot(all(qc_check_fig3f$LOD_flag != TRUE))
stopifnot(all(qc_check_fig3f$ipc_bcg_flag != TRUE))
stopifnot(!any(markers_to_plot %in% ipc_bcg_fails))

marker_levels <- marker_label_order

fig_dat <- npq_all %>%
  filter(
    !is.na(group),
    time_anno %in% c("Baseline", "Year 5")
  ) %>%
  dplyr::select(sno, time_anno, all_of(markers_to_plot)) %>%
  pivot_longer(
    cols = all_of(markers_to_plot),
    names_to = "marker",
    values_to = "value"
  ) %>%
  left_join(marker_dict_anno, by = c("marker" = "targetName")) %>%
  mutate(
    time_anno = case_when(
      time_anno == "Baseline" ~ "Baseline",
      time_anno == "Year 5" ~ "FIT",
      TRUE ~ time_anno
    ),
    time_anno = factor(time_anno, levels = c("Baseline", "FIT")),
    marker_lab = case_when(
      marker == "pTau_217" ~ "pTau-217",
      marker == "FCN2" ~ "FCN2",
      marker == "MME" ~ "MME",
      marker == "TARDBP" ~ "TDP43",
      marker == "VEGFA" ~ "VEGFA",
      marker == "GDF15" ~ "GDF15",
      TRUE ~ marker
    ),
    marker_lab = factor(marker_lab, levels = marker_label_order)
  ) %>%
  filter(!is.na(value))

paired_ids <- fig_dat %>%
  dplyr::distinct(sno, marker_lab, time_anno) %>%
  dplyr::count(sno, marker_lab) %>%
  filter(n == 2)

fig_dat <- fig_dat %>%
  semi_join(paired_ids, by = c("sno", "marker_lab"))

stat_df <- fig_dat %>%
  dplyr::select(sno, marker_lab, time_anno, value) %>%
  pivot_wider(names_from = time_anno, values_from = value) %>%
  filter(!is.na(Baseline), !is.na(FIT)) %>%
  group_by(marker_lab) %>%
  dplyr::summarise(
    n = n(),
    p_w = wilcox.test(FIT, Baseline, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_lab = paste0("P = ", formatC(p_w, format = "e", digits = 2))
  )

y_pos <- fig_dat %>%
  group_by(marker_lab) %>%
  dplyr::summarise(
    y = max(value, na.rm = TRUE) + 0.08 * diff(range(value, na.rm = TRUE)),
    .groups = "drop"
  )

stat_df <- stat_df %>%
  left_join(y_pos, by = "marker_lab")

### Plot

p_marker_change_dp <- ggplot(fig_dat, aes(x = time_anno, y = value)) +
  geom_line(
    aes(group = sno),
    color = "grey70",
    alpha = 0.5,
    linewidth = 0.1,
    na.rm = TRUE
  ) +
  geom_boxplot(
    aes(fill = time_anno),
    width = 0.42,
    outlier.shape = NA,
    color = "grey35",
    linewidth = 0.25,
    alpha = 0.35
  ) +
  geom_point(
    aes(fill = time_anno),
    shape = 21,
    size = 1.4,
    stroke = 0.1,
    color = "grey25",
    alpha = 0.65,
    position = position_jitter(width = 0.04, height = 0, seed = 1),
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.45,
    color = "black",
    fatten = 0
  ) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    linewidth = 0.75,
    color = "grey31"
  ) +
  stat_summary(
    aes(group = 1, fill = time_anno),
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.1,
    stroke = 0.7,
    color = "grey31"
  ) +
  geom_text(
    data = stat_df,
    aes(x = 1.5, y = y, label = p_lab),
    inherit.aes = FALSE,
    size = 3,
    color = "black"
  ) +
  facet_wrap(~ marker_lab, scales = "free_y", nrow = 2, ncol = 3) +
  scale_fill_manual(values = time_col, name = NULL) +
  labs(
    x = NULL,
    y = "NULISA Protein Quantification"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  )

ggsave(
  "results/figures/fig3f_dp_marker_change_baseline_fit.pdf",
  p_marker_change_dp,
  width = 5.5,
  height = 4.2
)

### Source data

source_data_fig3f_points <- fig_dat %>%
  mutate(pair_id = as.integer(factor(sno))) %>%
  transmute(
    pair_id,
    marker = as.character(marker_lab),
    time = as.character(time_anno),
    NULISA_protein_quantification = value
  ) %>%
  arrange(marker, pair_id, time)

write_xlsx(
  list(
    source_data_fig3f_points = source_data_fig3f_points
  ),
  "results/source_data/source_data_fig3f.xlsx"
)