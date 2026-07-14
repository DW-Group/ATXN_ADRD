library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggrepel)
library(readxl)
library(writexl)
library(cowplot)
library(patchwork)
library(pheatmap)

##---------------------------------------------------------
## Supp Fig 1
##---------------------------------------------------------

rm(list=ls())

load("data/processed_data_nhs.RData")

marker_dict$targetName_orig <- gsub("β","b",marker_dict$targetName_orig)

cor_long <- read_excel("results/nhs_NULISA_spearman_corr.xlsx") %>% 
  filter(!(var1 %in% c(scoreall_list)),
         !(var2 %in% c(scoreall_list)),
         !grepl("_log2|AT_comp",var1),
         !grepl("_log2|AT_comp",var2)) %>% 
  left_join(marker_dict %>% select(var1=targetName, var1_anno=targetName_orig), by="var1") %>% 
  left_join(marker_dict %>% select(var2=targetName, var2_anno=targetName_orig), by="var2")
length(unique(cor_long$var1))

cor_mat <- cor_long %>%
  select(var1_anno, var2_anno, r) %>%
  pivot_wider(names_from = var2_anno, values_from = r) %>%
  column_to_rownames("var1_anno") %>%
  as.matrix()
cor_mat[is.na(cor_mat)] <- 0

p_mat <- cor_long %>%
  select(var1_anno, var2_anno, p) %>%
  pivot_wider(names_from = var2_anno, values_from = p) %>%
  column_to_rownames("var1_anno") %>%
  as.matrix()
p_mat[is.na(p_mat)] <- 1
diag(p_mat) <- 0

pdf("results/figures/nhs_NULISA_spearman_corr.pdf", width = 16, height = 16)
pheatmap(
  cor_mat,
  color = colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100),
  breaks = seq(-1, 1, length.out = 101),
  clustering_method = "complete",
  border_color = NA,
  display_numbers = ifelse(p_mat < 0.05, "+", ""),
  number_color = "black",
  fontsize = 9.5,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  main = ""
)
dev.off()

### Source data

source_data_suppfig1_correlations <- cor_long %>%
  dplyr::transmute(
    protein_1 = var1_anno,
    protein_2 = var2_anno,
    spearman_r = r,
    p_value = p
  ) %>%
  dplyr::arrange(protein_1, protein_2)

write_xlsx(source_data_suppfig1_correlations, "results/source_data/source_data_supplementary_fig1.xlsx")

##---------------------------------------------------------
## Supp Fig 2
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
  filter(score %in% c("MIND", "aMED", "aMED component")) %>%
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

### Correlation arrows for diet

arrow_vars <- nhs_dietcomp_dict %>%
  filter(score %in% c("MIND", "aMED", "aMED component")) %>%
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

### Use one consistent scaling factor

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
  filter(score %in% c("MIND", "aMED", "aMED component"))

### PCA plot

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
  "results/figures/supp_fig2_pca_apoe.pdf",
  p_pc_apoe,
  width = 5,
  height = 5,
  dpi = 600,
  bg = "white"
)

### Source data

source_data_suppfig2_points <- pc_scores %>%
  dplyr::transmute(
    PC1,
    PC2,
    APOE4_status = as.character(APOE4_anno),
    ADRD_PRS = PGS002280_scaled_pool_rsq0.3
  ) 

source_data_suppfig2_arrows <- pco_corr_dat %>%
  dplyr::transmute(
    diet_label = diet_anno,
    correlation_with_PC1 = corr_pc1,
    correlation_with_PC2 = corr_pc2,
    arrow_x,
    arrow_y
  )

write_xlsx(
  list(
    suppfig2_points = source_data_suppfig2_points,
    suppfig2_arrows = source_data_suppfig2_arrows
  ),
  "results/source_data/source_data_supplementary_fig2.xlsx"
)

##---------------------------------------------------------
## Supp Fig 3
##---------------------------------------------------------

# plot_qq_pqtl_34.py

##---------------------------------------------------------
## Supp Fig 4
##---------------------------------------------------------

rm(list = ls())

##------------------
## Read selected TREM2-MS4A locus and GTEx brain QTL hits
##------------------

trem2_locus_file <- "results/hotspot_TREM2_MS4A/TREM2_MS4A_selected_locus.csv"
qtl_file <- "results/gtex_v8_brain_qtl_region_scan/gtex_v8_brain_eqtl_sqtl_region_hits_detail.csv"

stopifnot(file.exists(trem2_locus_file))
stopifnot(file.exists(qtl_file))

trem2_ms4a_loci <- fread(trem2_locus_file)
qtl_hits <- fread(qtl_file)

trem2_qtl <- qtl_hits %>%
  filter(locus_id %in% trem2_ms4a_loci$locus_id) %>%
  filter(
    !is.na(gene_name),
    gene_name != "",
    !is.na(tissue),
    !is.na(qtl_type),
    !is.na(pval_nominal)
  ) %>%
  mutate(
    qtl_type = ifelse(qtl_type == "eQTL", "eQTL", "sQTL"),
    neglog10_p = -log10(pval_nominal)
  )

stopifnot(nrow(trem2_qtl) > 0)

##------------------
## Brain tissue labels
##------------------

brain_tissues_order <- c(
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Substantia_nigra"
)

brain_tissue_labels <- c(
  "Brain_Amygdala" = "Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24" = "Anterior cingulate cortex",
  "Brain_Caudate_basal_ganglia" = "Caudate",
  "Brain_Cerebellar_Hemisphere" = "Cerebellar hemisphere",
  "Brain_Cerebellum" = "Cerebellum",
  "Brain_Cortex" = "Cerebral cortex",
  "Brain_Frontal_Cortex_BA9" = "Frontal cortex",
  "Brain_Hippocampus" = "Hippocampus",
  "Brain_Hypothalamus" = "Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia" = "Nucleus accumbens",
  "Brain_Putamen_basal_ganglia" = "Putamen",
  "Brain_Spinal_cord_cervical_c-1" = "Spinal cord",
  "Brain_Substantia_nigra" = "Substantia nigra"
)

trem2_qtl <- trem2_qtl %>%
  filter(tissue %in% brain_tissues_order)

##------------------
## Define genes to plot
##------------------

gene_track_genes <- c(
  "GIF", "OOSP2", "MS4A3", "MS4A2", "MS4A6A", "MS4A4E",
  "MS4A4A", "MS4A7", "MS4A14", "MS4A5", "MS4A1", "MS4A12",
  "MS4A13", "OR4D9", "OSBP", "PATL1", "OR10V1", "STX3",
  "MRPL16", "TCN1", "OOSP1", "MS4A6E", "ZP1", "MS4A8",
  "MS4A18", "MS4A15", "MS4A10", "CCDC86", "PRPF19",
  "CD6", "CD5", "PTGDR2", "TMEM109", "TMEM132A", "SLC15A3"
)

all_qtl_genes <- trem2_qtl %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(gene_name) %>%
  arrange(gene_name) %>%
  pull(gene_name)

focus_genes <- c(
  gene_track_genes,
  setdiff(all_qtl_genes, gene_track_genes)
)

##------------------
## Build matrix
##------------------

matrix_dat <- trem2_qtl %>%
  dplyr::group_by(gene_name, tissue, qtl_type) %>%
  dplyr::summarise(
    min_p = min(pval_nominal, na.rm = TRUE),
    max_neglog10_p = max(neglog10_p, na.rm = TRUE),
    top_variant = variant_id[which.min(pval_nominal)][1],
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = qtl_type,
    values_from = c(min_p, max_neglog10_p, top_variant),
    values_fill = NA
  ) %>%
  mutate(
    regulatory_class = case_when(
      !is.na(min_p_eQTL) & !is.na(min_p_sQTL) ~ "Both",
      !is.na(min_p_eQTL) ~ "eQTL",
      !is.na(min_p_sQTL) ~ "sQTL",
      TRUE ~ "None"
    ),
    score = pmax(
      max_neglog10_p_eQTL,
      max_neglog10_p_sQTL,
      na.rm = TRUE
    )
  )

matrix_grid <- expand.grid(
  gene_name = focus_genes,
  tissue = brain_tissues_order,
  stringsAsFactors = FALSE
)

matrix_dat_full <- matrix_grid %>%
  left_join(matrix_dat, by = c("gene_name", "tissue")) %>%
  mutate(
    regulatory_class = ifelse(is.na(regulatory_class), "None", regulatory_class),
    score = ifelse(is.na(score) | is.infinite(score), 0, score),
    gene_name = factor(gene_name, levels = focus_genes),
    tissue = factor(tissue, levels = rev(brain_tissues_order))
  )

stopifnot(nrow(matrix_dat_full) == length(focus_genes) * length(brain_tissues_order))

##------------------
## Plot
##------------------

qtl_cols <- c(
  "eQTL" = "#C8BEDD",
  "sQTL" = "#D8B28A",
  "Both" = "#A8B8B6"
)

plot_dat <- matrix_dat_full %>%
  filter(regulatory_class != "None", score > 0) %>%
  mutate(
    regulatory_class = factor(regulatory_class, levels = c("eQTL", "sQTL", "Both"))
  )

n_genes <- length(focus_genes)

p_supp_fig4 <- ggplot(
  matrix_dat_full,
  aes(x = gene_name, y = tissue)
) +
  geom_point(
    shape = 21,
    size = 0.9,
    fill = "#EEE9E1",
    color = "#D8D0C5",
    stroke = 0.18
  ) +
  geom_point(
    data = plot_dat,
    aes(fill = regulatory_class, size = score),
    shape = 21,
    color = "#3A3A3A",
    stroke = 0.22,
    alpha = 0.92
  ) +
  scale_fill_manual(
    values = qtl_cols,
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.2, 5.5),
    name = expression(max~ -log[10](P))
  ) +
  scale_x_discrete(
    expand = expansion(add = c(0.8, 0.8))
  ) +
  scale_y_discrete(
    labels = brain_tissue_labels,
    expand = expansion(add = c(0.8, 0.8))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(
      size = 7.2,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 8.5,
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "#E6E6E6",
      linewidth = 0.25
    ),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.key = element_blank(),
    legend.key.size = unit(5, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  "results/figures/supp_fig4_TREM2_MS4A_GTEx_brain_QTL_all_locus_genes.pdf",
  p_supp_fig4,
  width = 12,
  height = 4.2
)

### Source data

source_data_suppfig4_matrix <- matrix_dat_full %>%
  dplyr::transmute(
    gene = as.character(gene_name),
    brain_tissue_label = brain_tissue_labels[as.character(tissue)],
    regulatory_class,
    score,
    min_p_eQTL,
    min_p_sQTL,
    max_neglog10_p_eQTL,
    max_neglog10_p_sQTL
  ) %>%
  dplyr::arrange(gene, brain_tissue)

write_xlsx(
  list(
    suppfig4_matrix = source_data_suppfig4_matrix
  ),
  "results/source_data/source_data_supplementary_fig4.xlsx"
)
