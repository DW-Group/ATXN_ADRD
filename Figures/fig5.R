library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(ggVennDiagram)
library(eulerr)
library(locuszoomr)
library(readxl)
library(cowplot)
library(scales)
library(TwoSampleMR)
library(writexl)
library(stringr)
library(tidyr)
library(patchwork)

##---------------------------------------------------------
## Fig 5a
##---------------------------------------------------------

# plot_mht_gwaslab_top1.py

### Source data

infile <- "results/mht_extract/combined_mht_sumstats_top1_per_snp.txt"
out_file <- "results/source_data/source_data_fig5a.xlsx"

df <- fread(infile) %>%
  mutate(
    Chr  = as.numeric(Chr),
    Pos  = as.numeric(Pos),
    Pval = as.numeric(Pval)
  ) %>%
  filter(
    !is.na(Chr),
    !is.na(Pos),
    !is.na(Pval),
    Pval > 0,
    Pval <= 1
  ) %>%
  mutate(
    Chr = as.integer(Chr),
    Pos = as.integer(Pos),
    rsid = ifelse(is.na(rsid) | rsid == "", ".", as.character(rsid)),
    Trait = ifelse(is.na(Trait) | Trait == "", "NA", as.character(Trait)),
    SNPID = paste0(Chr, ":", Pos),
    LABEL = paste0(Trait, " | ", rsid),
    neg_log10_p = -log10(Pval)
  ) %>%
  arrange(Pval)

sig_df <- df %>%
  dplyr::filter(Pval < 5e-8) %>%
  dplyr::arrange(Trait, Chr, Pos, Pval)

locus_window <- 500000

sig_df <- sig_df %>%
  dplyr::group_by(Trait) %>%
  dplyr::arrange(Chr, Pos, Pval) %>%
  dplyr::mutate(
    new_locus = dplyr::case_when(
      dplyr::row_number() == 1 ~ 1L,
      Chr != dplyr::lag(Chr) ~ 1L,
      Pos - dplyr::lag(Pos) > locus_window ~ 1L,
      TRUE ~ 0L
    ),
    locus_id = cumsum(new_locus)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-new_locus)

source_data_fig5a_annotated_loci <- sig_df %>%
  arrange(Trait, locus_id, Pval) %>%
  group_by(Trait, locus_id) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    protein = Trait,
    rsid,
    SNPID,
    chr = Chr,
    position = Pos,
    p_value = Pval,
    neg_log10_p
  ) %>%
  arrange(p_value)

write_xlsx(source_data_fig5a_annotated_loci, out_file)

##---------------------------------------------------------
## Fig 5b
##---------------------------------------------------------

rm(list=ls())

source("scripts/source_rm_markers.R")
lod_dict <- read_excel("data/LOD_dictionary.xlsx")
marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx") %>% 
  left_join(lod_dict, by="targetName") %>% 
  mutate(ipc_bcg_flag = ifelse(targetName %in% ipc_bcg_fails, T, F))
qc_prot <- unique(marker_dict_final %>% filter(!LOD_flag_all & !ipc_bcg_flag) %>% pull(targetName))

protall <- fread("results/prot_for_mht.txt", header=F)$V1
lambdagc <- read_excel("results/ldsc/ldsc_h2_res.xlsx") %>% 
  filter(protein %in% qc_prot)

p_lambda <- lambdagc %>%
  ggplot(aes(x = lambda_gc)) +
  geom_histogram(
    fill = "gray80",
    color = "black",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = median(lambdagc$lambda_gc, na.rm = TRUE),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  labs(
    x = expression(lambda[GC]),
    y = "Number of plasma protein GWASs"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black", size=10),
    axis.title = element_text(color = "black", size=11),
    panel.grid = element_blank()
  )

ggsave("results/figures/fig5b_lambda_gc_hist.pdf", p_lambda, width = 2.5, height = 3)

### Source data

source_data_fig5b_lambda <- lambdagc %>%
  transmute(
    protein,
    lambda_GC = lambda_gc
  ) %>%
  arrange(protein)

write_xlsx(
  list(
    fig5b_lambda_gc = source_data_fig5b_lambda
  ),
  "results/source_data/source_data_fig5b.xlsx"
)

##---------------------------------------------------------
## Fig 5c
##---------------------------------------------------------

rm(list=ls())

### Load data

source("scripts/source_rm_markers.R")
lod_dict <- read_excel("data/LOD_dictionary.xlsx")
marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx") %>% 
  left_join(lod_dict, by="targetName") %>% 
  mutate(ipc_bcg_flag = ifelse(targetName %in% ipc_bcg_fails, T, F))
qc_prot <- unique(marker_dict_final %>% filter(!LOD_flag_all & !ipc_bcg_flag) %>% pull(targetName))

gwas_hits_all <- fread("results/pqtl_allsig_res.csv") %>% 
  filter(Trait %in% qc_prot,
         Pval < 5e-8) %>% 
  left_join(marker_dict_final, by=c("Trait"="targetName"), suffix=c("","_dict")) %>% 
  mutate(type = ifelse(chr==Chr & Pos > cis_start & Pos < cis_end, "cis", "trans"))

cis_res_dat <- fread("results/pqtl_cis_res.csv")
table(cis_res_dat$sig_Bonf)

### Count

cis_prot <- cis_res_dat %>%
  filter(sig_Bonf == TRUE) %>%
  pull(Trait) %>%
  unique()

trans_prot <- gwas_hits_all %>%
  filter(type == "trans") %>%
  pull(Trait) %>%
  unique()

all_prot <- unique(c(cis_prot, trans_prot))

length(cis_prot)
length(trans_prot)
length(intersect(cis_prot, trans_prot))

cis_only   <- length(setdiff(cis_prot, trans_prot))
trans_only <- length(setdiff(trans_prot, cis_prot))
both       <- length(intersect(cis_prot, trans_prot))

cis_col   <- "#9FB2C6"   
both_col  <- "#C7CCD2"  
trans_col <- "#D8D1C7" 

plot_df <- data.frame(
  category = factor(
    c("cis-QTL", "Both", "trans-QTL"),
    levels = c("trans-QTL", "Both", "cis-QTL")
  ),
  n = c(cis_only, both, trans_only)
)

total_n <- sum(plot_df$n)

plot_df$label <- plot_df$n 
plot_df$cum_n <- cumsum(plot_df$n)
plot_df$mid_n <- plot_df$cum_n - plot_df$n / 2

p_bar <- ggplot(plot_df, aes(x = n, y = "", fill = category)) +
  geom_col(
    color = "black",
    linewidth = 0.4,
    height = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "cis-QTL" = cis_col,
      "Both" = both_col,
      "trans-QTL" = trans_col
    )
  ) +
  geom_text(
    aes(label = label),
    color="black",
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  
  geom_text(
    data = plot_df,
    aes(x = mid_n, y = 1.5, label = category),
    inherit.aes = FALSE,
    vjust = -1,
    size = 6
  ) +
  
  annotate(
    "text",
    x = sum(plot_df$n) / 2,
    y = 0,
    label = "Number of plasma proteins with QTL",
    size = 6.5
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.margin = margin(15, 5, 20, 5) 
  )

ggsave(
  "results/figures/fig5c_cis_trans_bar.pdf",
  p_bar,
  width = 5,
  height = 1.2
)

### Source data

source_data_fig5c_proteins <- tibble(
  protein = all_prot
) %>%
  mutate(
    QTL_category = case_when(
      protein %in% cis_prot & protein %in% trans_prot ~ "Both",
      protein %in% cis_prot ~ "cis-QTL",
      protein %in% trans_prot ~ "trans-QTL",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(QTL_category, protein)

write_xlsx(
  list(
    fig5c_proteins = source_data_fig5c_proteins
  ),
  "results/source_data/source_data_fig5c.xlsx"
)

##---------------------------------------------------------
## Fig 5d
##---------------------------------------------------------

rm(list=ls())

### Extract hotspot region

locus_sum <- fread("results/pqtl_hotspot.csv")

hotspot_row <- locus_sum %>%
  filter(hotspot_class2 == "trans hotspot") %>%
  slice(1)

hotspot_prot_list <- hotspot_row$trans_proteins %>%
  str_split(";\\s*") %>%
  .[[1]] %>%
  unique()

hotspot_row$locus_start
hotspot_row$locus_end

hotspot_chr <- hotspot_row$Chr
hotspot_start <- hotspot_row$locus_start
hotspot_end <- hotspot_row$locus_end

gwas_hits_region <- fread("results/pqtl_allsig_res_anno.csv") %>%
  filter(
    Trait %in% hotspot_prot_list,
    type == "trans",
    Chr == hotspot_chr,
    Pos >= hotspot_start,
    Pos <= hotspot_end
  )

lead_snp_shared <- gwas_hits_region %>%
  group_by(rsid) %>%
  summarise(
    n_protein = n_distinct(Trait),
    best_p = min(Pval)
  ) %>%
  arrange(desc(n_protein), best_p) %>%
  slice(1)

### Get LD info

lead_snp <- "rs66530140"
lead_chr <- 4

ld_raw <- LDproxy(
  snp = lead_snp,
  pop = "CEU",
  r2d = "r2",
  token = ld_token,
  genome_build = "grch37",
  win_size = 1000000   
)

ld_df <- ld_raw %>%
  as.data.frame() %>%
  transmute(
    rsid = RS_Number,
    r2 = as.numeric(R2)
  ) %>% 
  distinct()

length(ld_df$rsid); length(unique(ld_df$rsid))

ld_df <- ld_df %>%
  group_by(rsid) %>%
  summarise(r2 = max(r2, na.rm = TRUE), .groups = "drop") %>% 
  filter(rsid != ".")

lead_pos <- unique(gwas_hits_region %>% filter(rsid == "rs66530140") %>% pull(Pos))
window <- 5e5
start <- lead_pos - window
end   <- lead_pos + window

res_list <- list()
for (tmp_prot in hotspot_prot_list){
  
  file <- paste0(
    "results/trans_hotspot_region/",
    tmp_prot,
    ".trans_hotspot_chr4_186152506_188152506.txt"
  )
  
  tmp <- fread(file) %>%
    mutate(rsid = as.character(rsid)) %>%
    left_join(ld_df, by = "rsid") %>% 
    filter(
      Chr == lead_chr,
      Pos >= start,
      Pos <= end
    ) %>%
    mutate(protein = tmp_prot)
  
  res_list[[tmp_prot]] <- tmp
}

trans_all_anno <- as.data.frame(do.call("rbind", res_list))
head(trans_all_anno)
table(is.na(trans_all_anno$r2))

### Plot

# trans_hotspot.R

### Source data

load("results/pqtl_hotspot.RData")
source_data_fig5d_region <- trans_all_anno %>%
  transmute(
    protein,
    chr = Chr,
    position = Pos,
    rsid,
    p_value = Pval,
    neg_log10_p = -log10(Pval),
    LD_r2_to_lead = r2,
    lead_snp = lead_snp
  ) %>%
  arrange(protein, position)

write_xlsx(
  list(
    fig5d_region_points = source_data_fig5d_region
  ),
  "results/source_data/source_data_fig5d.xlsx"
)

##---------------------------------------------------------
## Fig 5e 
##---------------------------------------------------------

rm(list=ls())

### Load dictionary

source("scripts/source_rm_markers.R")
lod_dict <- read_excel("data/LOD_dictionary.xlsx")
marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx") %>% 
  left_join(lod_dict, by="targetName") %>% 
  mutate(ipc_bcg_flag = ifelse(targetName %in% ipc_bcg_fails, T, F))
qc_prot <- unique(marker_dict_final %>% filter(!LOD_flag_all & !ipc_bcg_flag) %>% pull(targetName))

### Trans

coloc_sig_all <- read_excel("results/colocalization/coloc_sig_res.xlsx") %>% 
  mutate(
    protein = sub("_rs.*$", "", targetName),
    rsid    = sub("^.*_(rs[0-9]+)$", "\\1", targetName)
  ) %>%
  filter(protein %in% qc_prot)

region_list <- lapply(unique(coloc_sig_all$protein), function(prot) {
  reg_path <- file.path("results/remeta_NULISA/sig_region/", paste0(prot, ".sig_region.txt"))
  tmp_region <- fread(reg_path)
  tmp_region$protein <- prot
  tmp_region
})
tmp_region_all <- rbindlist(region_list, fill = TRUE) %>% 
  select(rsid, Chr, Pos) %>% 
  filter(rsid %in% unique(coloc_sig_all$rsid)) %>% 
  distinct()

coloc_sig_anno <- coloc_sig_all %>% 
  left_join(tmp_region_all, by="rsid") %>% 
  left_join(marker_dict_final, by=c("protein"="targetName"), suffix=c("","_dict")) %>% 
  mutate(type = ifelse(chr==Chr & Pos > cis_start & Pos < cis_end, "cis", "trans"))

table(coloc_sig_anno$type)
coloc_trans_anno <- coloc_sig_anno %>% 
  filter(type == "trans") %>% 
  select(type, protein, outcome,  
         PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, PP3_plus_PP4, PP4_div_PP3_plus_PP4)

### Cis

coloc_cis_all <- read_excel("results/colocalization/coloc_cis_res.xlsx") %>% 
  filter(targetName %in% qc_prot)
tmp <- as.data.frame(table(coloc_cis_all$targetName))

coloc_cis_anno <- coloc_cis_all %>% 
  mutate(type = "cis",
         protein = targetName) %>% 
  select(type, protein, outcome, 
         PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, PP3_plus_PP4, PP4_div_PP3_plus_PP4)

### Scatter plot

coloc_anno <- as.data.frame(rbind(coloc_cis_anno, coloc_trans_anno))

fig_dat <- coloc_anno %>%
  mutate(
    outcome_label = recode(
      outcome,
      "bellenguez_ad"   = "AD/ADRD (Bellenguez et al.)",
      "wightman_ad"     = "Clinical AD (Wightman et al.)",
      "wightmanukb_ad"  = "AD + proxy (Wightman et al.)",
      "lee"             = "Cognitive performance (Lee et al.)"
    ),
    type_label = recode(
      type,
      "cis"   = "cis-pQTL",
      "trans" = "trans-pQTL"
    ),
    score = PP4_div_PP3_plus_PP4 + PP.H3.abf
  )

label_df <- fig_dat %>%
  filter(PP4_div_PP3_plus_PP4 >= 0.6 | PP.H3.abf >= 0.6) %>% 
  distinct(protein, outcome, .keep_all = TRUE) %>%
  mutate(label = protein)

outcome_col <- c(
  "AD/ADRD (Bellenguez et al.)"            = "#5F7F82",
  "Clinical AD (Wightman et al.)" = "#8FA8B8",
  "AD + proxy (Wightman et al.)"     = "#C4D2C8",
  "Cognitive performance (Lee et al.)" = "#B89C6B"
)

p_coloc <- ggplot(
  fig_dat,
  aes(
    x = PP4_div_PP3_plus_PP4,
    y = PP.H3.abf
  )
) +
  geom_point(
    aes(fill = outcome_label, shape = type_label),
    color = "#F2F0EF",
    size = 3,
    stroke = 0.25,
    alpha = 0.85
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = label),
    size = 3.0,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.35,
    segment.color = "grey",
    segment.size = 0.2,
    force = 2,                
    force_pull = 0.5,
    max.iter = 5000
  ) +
  scale_x_continuous(
    name = "Conditional PP.H4",
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(add = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    name = "PP.H3",
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(add = c(0.05, 0.05))
  ) +
  scale_fill_manual(
    values = outcome_col,
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "cis-pQTL" = 21,
      "trans-pQTL" = 23
    ),
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        color = "gray31",   
        size = 3
      )
    ),
    shape = guide_legend(
      override.aes = list(
        fill = "grey80",
        color = "gray31",
        size = 3
      )
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
     
    legend.position = "right",
    legend.background = element_blank()
  ) +
  coord_equal()

ggsave("results/figures/fig5e_coloc_scatter.pdf", p_coloc, width = 7, height = 7)

### Source data

source_data_fig5e_points <- fig_dat %>%
  transmute(
    protein,
    outcome = outcome_label,
    QTL_type = type_label,
    conditional_PP_H4 = PP4_div_PP3_plus_PP4,
    PP_H3 = PP.H3.abf,
    PP_H0 = PP.H0.abf,
    PP_H1 = PP.H1.abf,
    PP_H2 = PP.H2.abf,
    PP_H4 = PP.H4.abf,
    PP3_plus_PP4
  ) %>%
  arrange(QTL_type, outcome, protein)

write_xlsx(
  list(
    fig5e_points = source_data_fig5e_points
  ),
  "results/source_data/source_data_fig5e.xlsx"
)

##---------------------------------------------------------
## Fig 5f Coloc callout
##---------------------------------------------------------

# coloc_callout.R 

##---------------------------------------------------------
## Fig 5g
##---------------------------------------------------------

rm(list = ls())

##-------------------------------
## TREM2-MS4A locus
##-------------------------------

trans_loci <- fread("results/tissue/pqtl_query_loci_from_existing_coloc_files.csv") %>% 
  filter(pqtl_class=="trans")

trem2_loci <- trans_loci %>%
  dplyr::filter(protein == "TREM2") %>%
  as.data.frame() %>%
  dplyr::arrange(lead_p)

trem2_loci %>%
  dplyr::select(
    protein,
    protein_gene,
    lead_rsid,
    lead_chr,
    lead_pos,
    lead_beta,
    lead_p,
    region_file,
    locus_id
  )

trem2_ms4a_loci <- trem2_loci %>%
  dplyr::filter(lead_chr == 11) %>%
  dplyr::arrange(lead_p)

trem2_ms4a_loci <- trem2_ms4a_loci[1, ]

##-------------------------------
## GTEx brain eQTL/sQTL hits for selected locus
##-------------------------------

qtl_hits <- fread(
  "results/gtex_v8_brain_qtl_region_scan/gtex_v8_brain_eqtl_sqtl_region_hits_detail.csv"
)

trem2_qtl <- qtl_hits %>%
  dplyr::filter(locus_id %in% trem2_ms4a_loci$locus_id) %>%
  dplyr::mutate(
    neglog10_p = -log10(pval_nominal)
  )

##-------------------------------
## Gene-level summary
##-------------------------------

trem2_qtl_summary_gene <- trem2_qtl %>%
  dplyr::filter(!is.na(gene_name), gene_name != "") %>%
  dplyr::group_by(qtl_type, gene_name) %>%
  dplyr::summarise(
    n_tissues = dplyr::n_distinct(tissue),
    min_p = min(pval_nominal, na.rm = TRUE),
    max_neglog10_p = max(neglog10_p, na.rm = TRUE),
    top_tissue = tissue[which.min(pval_nominal)][1],
    top_variant = variant_id[which.min(pval_nominal)][1],
    .groups = "drop"
  ) %>%
  dplyr::arrange(qtl_type, min_p)

##-------------------------------
## Brain tissue labels
##-------------------------------

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
  "Brain_Anterior_cingulate_cortex_BA24" = "ACC",
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

##-------------------------------
## Genes in MS4A locus gene track
##-------------------------------

gene_track_genes <- c(
  "GIF", "OOSP2", "MS4A3", "MS4A2", "MS4A6A", "MS4A4E",
  "MS4A4A", "MS4A7", "MS4A14", "MS4A5", "MS4A1", "MS4A12",
  "MS4A13", "OR4D9", "OSBP", "PATL1", "OR10V1", "STX3",
  "MRPL16", "TCN1", "OOSP1", "MS4A6E", "ZP1", "MS4A8",
  "MS4A18", "MS4A15", "MS4A10", "CCDC86", "PRPF19",
  "CD6", "CD5", "PTGDR2", "TMEM109", "TMEM132A", "SLC15A3"
)

gene_order_dat <- trem2_qtl %>%
  dplyr::filter(
    !is.na(gene_name),
    gene_name != "",
    gene_name %in% gene_track_genes
  ) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(
    min_p = min(pval_nominal, na.rm = TRUE),
    max_score = max(neglog10_p, na.rm = TRUE),
    n_tissues = dplyr::n_distinct(tissue),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    gene_name = factor(gene_name, levels = gene_track_genes)
  ) %>%
  dplyr::arrange(gene_name) %>%
  dplyr::mutate(
    gene_name = as.character(gene_name)
  )

focus_ms4a_genes <- gene_order_dat %>%
  dplyr::pull(gene_name)

focus_ms4a_genes

trem2_qtl_focus <- trem2_qtl %>%
  dplyr::filter(gene_name %in% focus_ms4a_genes)

trem2_qtl_focus_summary <- trem2_qtl_summary_gene %>%
  dplyr::filter(gene_name %in% focus_ms4a_genes)

##-------------------------------
## Build matrix
##-------------------------------

matrix_dat <- trem2_qtl_focus %>%
  dplyr::filter(!is.na(gene_name), gene_name != "") %>%
  dplyr::group_by(gene_name, tissue, qtl_type) %>%
  dplyr::summarise(
    min_p = min(pval_nominal, na.rm = TRUE),
    max_neglog10_p = max(neglog10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    qtl_type = ifelse(qtl_type == "eQTL", "eQTL", "sQTL")
  ) %>%
  tidyr::pivot_wider(
    names_from = qtl_type,
    values_from = c(min_p, max_neglog10_p),
    values_fill = NA
  ) %>%
  dplyr::mutate(
    regulatory_class = dplyr::case_when(
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
  gene_name = focus_ms4a_genes,
  tissue = brain_tissues_order,
  stringsAsFactors = FALSE
)

matrix_dat_full <- matrix_grid %>%
  dplyr::left_join(matrix_dat, by = c("gene_name", "tissue")) %>%
  dplyr::mutate(
    regulatory_class = ifelse(is.na(regulatory_class), "None", regulatory_class),
    score = ifelse(is.na(score) | is.infinite(score), 0, score),
    gene_name = factor(gene_name, levels = focus_ms4a_genes),
    tissue = factor(tissue, levels = rev(brain_tissues_order))
  )

##-------------------------------
## Plot matrix
##-------------------------------

qtl_cols <- c(
  "eQTL" = "#C8BEDD",
  "sQTL" = "#D8B28A",
  "Both" = "#A8B8B6"
)

plot_dat <- matrix_dat_full %>%
  dplyr::filter(regulatory_class != "None", score > 0) %>%
  dplyr::mutate(
    regulatory_class = factor(regulatory_class, levels = c("eQTL", "sQTL", "Both"))
  )

p_trem2_ms4a_matrix <- ggplot(
  matrix_dat_full,
  aes(x = gene_name, y = tissue)
) +
  geom_point(
    shape = 21,
    size = 1.05,
    fill = "#EEE9E1",
    color = "#D8D0C5",
    stroke = 0.20
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
    range = c(1.4, 5.0),
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
      size = 8,
      color = "#2F2F2F",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 8,
      color = "#2F2F2F"
    ),
    panel.border = element_rect(
      color = "#3A3A3A",
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
    legend.position = "left",
    legend.title = element_text(
      size = 7,
      color = "#2F2F2F"
    ),
    legend.text = element_text(
      size = 7,
      color = "#2F2F2F"
    ),
    legend.key = element_blank(),
    legend.key.size = unit(5, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(
  "results/figures/fig5g_TREM2_MS4A_brain_QTL.pdf",
  p_trem2_ms4a_matrix,
  width = 5.2,
  height = 2.8
)

### Source data

source_data_fig5g_matrix <- matrix_dat_full %>%
  transmute(
    gene = as.character(gene_name),
    brain_tissue_label = brain_tissue_labels[as.character(tissue)],
    regulatory_class,
    score = score,
    min_p_eQTL,
    min_p_sQTL,
    max_neglog10_p_eQTL,
    max_neglog10_p_sQTL
  ) %>%
  arrange(gene, brain_tissue)

write_xlsx(
  list(
    fig5g_matrix = source_data_fig5g_matrix
  ),
  "results/source_data/source_data_fig5g.xlsx"
)

##---------------------------------------------------------
## Fig 5h
##---------------------------------------------------------

rm(list = ls())

mr_res <- read_excel("results/mr/mr_res.xlsx")
anno_res <- fread("results/cis_region_sig_in_outcome_gwas.csv") %>% 
  mutate(
    sig_wightman_bonf    = ifelse(minp_wightman <= 0.05 / nsnp_wightman, TRUE, FALSE),
    sig_wightmanukb_bonf = ifelse(minp_wightmanukb <= 0.05 / nsnp_wightmanukb, TRUE, FALSE),
    sig_bellenguez_bonf  = ifelse(minp_bellenguez <= 0.05 / nsnp_bellenguez, TRUE, FALSE),
    sig_lee_bonf         = ifelse(minp_lee <= 0.05 / nsnp_lee, TRUE, FALSE),
    sig_any_bonf         = sig_wightman_bonf | sig_wightmanukb_bonf | sig_bellenguez_bonf | sig_lee_bonf
  )

mr_anno <- mr_res %>%
  left_join(
    anno_res %>% select(targetName, sig_any, sig_any_bonf),
    by = c("id.exposure" = "targetName")
  ) %>% 
  mutate(
    outcome_label = recode(
      outcome,
      "bellenguez_ad"  = "AD/ADRD (Bellenguez et al.)",
      "wightman_ad"    = "Clinical AD (Wightman et al.)",
      "wightmanukb_ad" = "AD + proxy (Wightman et al.)",
      "lee_cp"         = "Cognitive performance (Lee et al.)"
    ),
    sig = case_when(
      pval < 0.01 ~ "#",
      pval < 0.05 ~ "+",
      TRUE ~ ""
    )
  )

q_low  <- quantile(mr_anno$b, 0.1, na.rm = TRUE)
q_high <- quantile(mr_anno$b, 0.9, na.rm = TRUE)

mr_anno <- mr_anno %>%
  mutate(
    b_winsor = pmin(pmax(b, q_low), q_high)
  )

protein_levels <- unique(mr_anno$id.exposure)

mr_anno <- mr_anno %>%
  mutate(
    id.exposure   = factor(id.exposure, levels = rev(protein_levels)),
    outcome_label = factor(
      outcome_label,
      levels = c(
        "AD/ADRD (Bellenguez et al.)",
        "Clinical AD (Wightman et al.)",
        "AD + proxy (Wightman et al.)",
        "Cognitive performance (Lee et al.)"
      )
    )
  )

prot_anno <- mr_anno %>%
  distinct(id.exposure, sig_any_bonf)

p_mr_heat <- ggplot(mr_anno, aes(x = outcome_label, y = id.exposure)) +
  geom_tile(aes(fill = b_winsor), color = "white", linewidth = 0.2) +
  geom_text(aes(label = sig), color = "black", size = 2.5) +
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#9FB2C6",
    mid = "#F7F5F2",
    high = "#8C6A5A",
    midpoint = 0,
    limits = c(q_low, q_high),
    oob = squish
  ) +
  labs(fill = expression(beta["ATX(N) marker"])) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 6
    ),
    axis.text.y = element_text(color = "black", size = 6),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "left",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    panel.border = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  coord_fixed(ratio = 1)

p_side <- ggplot(prot_anno, aes(x = 1, y = id.exposure)) +  
  geom_tile(aes(fill = sig_any_bonf), width = 0.95, height = 0.95, linewidth = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(
    name = "AD/Cognition signal\nin cis region",
    values = c(
      "TRUE"  = "#7D8C72",
      "FALSE" = "#D6D2CB"
    ),
    breaks = c("TRUE", "FALSE"),  
    labels = c("TRUE" = "Yes", "FALSE" = "No"),
    guide = "none"
  ) +
  theme_void() +   
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

p_mr_heat_final <- p_mr_heat + p_side + 
  plot_layout(widths = c(1, 0.01)) 

ggsave(
  "results/figures/fig5_mr_heatmap.pdf",
  p_mr_heat_final,
  width = 5,
  height = 4
)

p_side_leg <- p_side +
  guides(
    fill = guide_legend(
      keywidth = unit(2, "mm"),
      keyheight = unit(4, "mm")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6),
    legend.key.size = unit(3, "mm")
  )

legend_side <- get_legend(p_side_leg)

ggsave(
  "results/figures/fig5_side_legend.pdf",
  legend_side,
  width = 3,
  height = 3
)

### Source data

source_data_fig5h_heatmap <- mr_anno %>%
  transmute(
    protein = as.character(id.exposure),
    outcome = as.character(outcome_label),
    MR_beta = b,
    MR_beta_winsorized_for_plot = b_winsor,
    p_value = pval
  ) %>%
  arrange(protein, outcome)

write_xlsx(
  list(
    fig5h_heatmap = source_data_fig5h_heatmap
  ),
  "results/source_data/source_data_fig5h.xlsx"
)

##---------------------------------------------------------
## Fig 5i
##---------------------------------------------------------

rm(list = ls())

load("data/mr.RData")

##-------------------------------
## Harmonize MR data 
##-------------------------------

harmonized_list <- list()
for (pheno in names(outcome_list)) {
  harmonized_list[[pheno]] <- harmonise_data(
    exposure_dat = instru_pqtl_list[[pheno]],
    outcome_dat  = outcome_list[[pheno]]
  )
}

CNTN2_inst <- harmonized_list[["CNTN2"]] %>%
  filter(outcome == "bellenguez_ad", mr_keep) %>%
  mutate(
    panel = "CNTN2",
    outcome_label = "AD (Bellenguez et al.)",
    exposure_label = "CNTN2"
  )

MSLN_inst <- harmonized_list[["MSLN"]] %>%
  filter(outcome == "lee_cp", mr_keep) %>%
  mutate(
    panel = "MSLN",
    outcome_label = "Cognitive performance",
    exposure_label = "MSLN"
  )

inst_dat <- bind_rows(CNTN2_inst, MSLN_inst) %>%
  mutate(
    SNP = as.character(SNP),
    p_exp = 2 * pnorm(-abs(beta.exposure / se.exposure)),
    xmin = beta.exposure - 1.96 * se.exposure,
    xmax = beta.exposure + 1.96 * se.exposure,
    ymin = beta.outcome  - 1.96 * se.outcome,
    ymax = beta.outcome  + 1.96 * se.outcome
  )

##-------------------------------
## Define index SNP 
##-------------------------------

index_snp_df <- inst_dat %>%
  group_by(panel) %>%
  slice_min(order_by = p_exp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(panel, index_snp = SNP)

print(index_snp_df)

##-------------------------------
## Helper function
##-------------------------------

read_coloc_scatter <- function(prot, out, panel, outcome_label) {
  
  f <- paste0("data/colocalization/coloc_cis/coloc_cis_", out, "_", prot, ".tsv")
  
  fread(f) %>%
    filter(
      is.finite(beta_exp),
      is.finite(varbeta_exp),
      is.finite(beta_out),
      is.finite(varbeta_out),
      !is.na(snp)
    ) %>%
    mutate(
      snp = as.character(snp),
      se_exp = sqrt(varbeta_exp),
      se_out = sqrt(varbeta_out),
      panel = panel,
      outcome_label = outcome_label,
      exposure_label = prot
    ) %>%
    distinct(snp, .keep_all = TRUE)
}

get_ld_to_index <- function(index_snp,
                            coloc_snps,
                            token,
                            pop = "EUR",
                            genome_build = "grch37",
                            win_size = 1000000,
                            outfile) {
    
  ld_raw <- LDproxy(
    snp = index_snp,
    pop = pop,
    r2d = "r2",
    token = token,
    genome_build = genome_build,
    win_size = win_size
  )
  
  ld_raw <- as.data.frame(ld_raw)
  names(ld_raw) <- gsub("\\s+", "_", names(ld_raw))
  
  rs_col <- intersect(
    c("RS_Number", "RS_Number_", "RS_number", "rsid", "RSID"),
    names(ld_raw)
  )[1]
  
  r2_col <- intersect(c("R2", "r2"), names(ld_raw))[1]

  ld_out <- ld_raw %>%
    transmute(
      snp1 = index_snp,
      snp2 = as.character(.data[[rs_col]]),
      r2   = suppressWarnings(as.numeric(.data[[r2_col]]))
    ) %>%
    filter(!is.na(snp2)) %>%
    distinct() %>%
    filter(snp2 %in% coloc_snps)
  
  if (!(index_snp %in% ld_out$snp2)) {
    ld_out <- bind_rows(
      data.frame(snp1 = index_snp, snp2 = index_snp, r2 = 1),
      ld_out
    )
  }
  
  missing_snps <- setdiff(coloc_snps, ld_out$snp2)
  if (length(missing_snps) > 0) {
    ld_out <- bind_rows(
      ld_out,
      data.frame(snp1 = index_snp, snp2 = missing_snps, r2 = NA_real_)
    )
  }
  
  fwrite(ld_out, outfile, sep = "\t")
  ld_out
}

read_ld_to_index <- function(ld_file, index_snp, panel_name) {
  
  ld <- fread(ld_file)
  names(ld) <- tolower(names(ld))
  
  if (!"r2" %in% names(ld)) {
    if ("r2_corr" %in% names(ld)) {
      ld$r2 <- ld$r2_corr
    } else if ("rsquare" %in% names(ld)) {
      ld$r2 <- ld$rsquare
    } else {
      stop("No r2 column found in ", ld_file)
    }
  }
  
  snp_col1 <- names(ld)[grepl("^snp|^rsid|^var1|^snp_a|^snp1$", names(ld))][1]
  snp_col2 <- names(ld)[grepl("^snp|^rsid|^var2|^snp_b|^snp2$", names(ld))][2]
  
  ld %>%
    mutate(
      snp1 = as.character(.data[[snp_col1]]),
      snp2 = as.character(.data[[snp_col2]])
    ) %>%
    filter(snp1 == index_snp | snp2 == index_snp) %>%
    mutate(
      snp = ifelse(snp1 == index_snp, snp2, snp1),
      panel = panel_name
    ) %>%
    select(panel, snp, r2) %>%
    distinct() %>%
    bind_rows(data.frame(panel = panel_name, snp = index_snp, r2 = 1)) %>%
    distinct()
}

##-------------------------------
## Read regional data
##-------------------------------

CNTN2_reg <- read_coloc_scatter(
  prot = "CNTN2",
  out = "bellenguez_ad",
  panel = "CNTN2",
  outcome_label = "AD/ADRD (Bellenguez et al.)"
)

MSLN_reg <- read_coloc_scatter(
  prot = "MSLN",
  out = "lee",
  panel = "MSLN",
  outcome_label = "Cognitive performance"
)

##-------------------------------
## Generate LD files
##-------------------------------

CNTN2_index <- index_snp_df %>%
  filter(panel == "CNTN2") %>%
  pull(index_snp)

MSLN_index <- index_snp_df %>%
  filter(panel == "MSLN") %>%
  pull(index_snp)

CNTN2_ld_file <- "results/CNTN2_bellenguez_ad_ld.tsv"
MSLN_ld_file  <- "results/MSLN_lee_ld.tsv"

##-------------------------------
## Read LD for plotting
##-------------------------------

CNTN2_ld <- read_ld_to_index(
  ld_file = CNTN2_ld_file,
  index_snp = CNTN2_index,
  panel_name = "CNTN2"
)

MSLN_ld <- read_ld_to_index(
  ld_file = MSLN_ld_file,
  index_snp = MSLN_index,
  panel_name = "MSLN"
)

ld_dat <- bind_rows(CNTN2_ld, MSLN_ld)

##-------------------------------
## Merge LD into regional data
##-------------------------------

reg_dat <- bind_rows(CNTN2_reg, MSLN_reg)

plot_reg <- reg_dat %>%
  left_join(ld_dat, by = c("panel", "snp")) %>%
  mutate(
    ld_bin = case_when(
      is.na(r2)      ~ "0-0.2",
      r2 < 0.2       ~ "0-0.2",
      r2 < 0.4       ~ "0.2-0.4",
      r2 < 0.6       ~ "0.4-0.6",
      r2 < 0.8       ~ "0.6-0.8",
      TRUE           ~ "0.8-1.0"
    )
  ) %>%
  left_join(
    inst_dat %>%
      select(panel, SNP) %>%
      distinct() %>%
      mutate(is_instrument = TRUE),
    by = c("panel", "snp" = "SNP")
  ) %>%
  mutate(is_instrument = ifelse(is.na(is_instrument), FALSE, is_instrument))

##-------------------------------
## Plotting function
##-------------------------------

ld_cols <- c(
  "0-0.2"   = "#E8E6E3",
  "0.2-0.4" = "#9FB9B9",
  "0.4-0.6" = "#6FA8A6",
  "0.6-0.8" = "#E6A86C",
  "0.8-1.0" = "#C65A5A"
)

plot_one_panel <- function(panel_name, width = 2.2, height = 2.3) {
  
  reg_sub <- plot_reg %>%
    filter(
      panel == panel_name,
      is.finite(beta_exp),
      is.finite(beta_out)
    )
  
  inst_sub <- inst_dat %>%
    filter(
      panel == panel_name,
      is.finite(beta.exposure),
      is.finite(beta.outcome),
      is.finite(xmin),
      is.finite(xmax),
      is.finite(ymin),
      is.finite(ymax)
    )
  
  exposure_label <- unique(inst_sub$exposure_label)
  outcome_label  <- unique(inst_sub$outcome_label)
  
  x_range <- diff(range(c(reg_sub$beta_exp, inst_sub$xmin, inst_sub$xmax), na.rm = TRUE))
  y_range <- diff(range(c(reg_sub$beta_out, inst_sub$ymin, inst_sub$ymax), na.rm = TRUE))
  
  cap_x <- 0.015 * x_range
  cap_y <- 0.015 * y_range
  
  p <- ggplot() +
    
    geom_point(
      data = reg_sub %>% filter(!is_instrument, is.na(r2) | r2 < 0.2),
      aes(x = beta_exp, y = beta_out),
      shape = 21,
      size = 1.9,
      stroke = 0.15,         
      fill = "#E6E6E6",      
      color = "#BFBFBF",      
      alpha = 0.35
    ) +
    
    geom_hline(yintercept = 0, linetype = 2, color = "gray60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray60", linewidth = 0.3) +
    
    geom_point(
      data = reg_sub %>%
        filter(!is_instrument, !is.na(r2), r2 >= 0.2) %>%
        arrange(r2),  
      aes(x = beta_exp, y = beta_out, fill = ld_bin),
      shape = 21,
      size = 1.9,
      stroke = 0.2,
      color = "gray31",
      alpha = 0.8
    ) +
    
    geom_segment(
      data = inst_sub,
      aes(x = xmin, xend = xmax, y = beta.outcome, yend = beta.outcome),
      linewidth = 0.2,
      color = "black"
    ) +
    geom_segment(
      data = inst_sub,
      aes(x = xmin, xend = xmin,
          y = beta.outcome - cap_y, yend = beta.outcome + cap_y),
      linewidth = 0.2,
      color = "black"
    ) +
    geom_segment(
      data = inst_sub,
      aes(x = xmax, xend = xmax,
          y = beta.outcome - cap_y, yend = beta.outcome + cap_y),
      linewidth = 0.2,
      color = "black"
    ) +
    
    geom_segment(
      data = inst_sub,
      aes(x = beta.exposure, xend = beta.exposure, y = ymin, yend = ymax),
      linewidth = 0.2,
      color = "black"
    ) +
    geom_segment(
      data = inst_sub,
      aes(x = beta.exposure - cap_x, xend = beta.exposure + cap_x,
          y = ymin, yend = ymin),
      linewidth = 0.2,
      color = "black"
    ) +
    geom_segment(
      data = inst_sub,
      aes(x = beta.exposure - cap_x, xend = beta.exposure + cap_x,
          y = ymax, yend = ymax),
      linewidth = 0.2,
      color = "black"
    ) +
    
    geom_point(
      data = inst_sub,
      aes(x = beta.exposure, y = beta.outcome),
      shape = 23,
      size = 2.2,
      stroke = 0.3,
      fill = "#9A8FBF",
      color = "black"
    ) +
    
    geom_text_repel(
      data = inst_sub,
      aes(x = beta.exposure, y = beta.outcome, label = SNP),
      size = 3,
      box.padding = 0.2,
      point.padding = 0.15,
      segment.color = "gray50",
      min.segment.length = 0,
      seed = 1
    ) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    scale_fill_manual(values = ld_cols, drop = FALSE) +
    labs(
      x = bquote(beta[.(exposure_label)]),
      y = bquote(beta[.(outcome_label)])
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 8.5, color = "black"),
      axis.title = element_text(size = 10),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
  
  out_file <- paste0("results/figures/fig5_mr_callout_regional_scatter_", panel_name, ".pdf")
  ggsave(out_file, p, width = width, height = height)
  
  p
}

##-------------------------------
## Make figures
##-------------------------------

p_CNTN2 <- plot_one_panel("CNTN2")
p_MSLN  <- plot_one_panel("MSLN")

### Source data

source_data_fig5i_regional_points <- plot_reg %>%
  transmute(
    panel,
    SNP = snp,
    exposure = exposure_label,
    outcome = outcome_label,
    beta_exposure = beta_exp,
    beta_outcome = beta_out,
    LD_bin = ld_bin
  ) %>%
  arrange(panel, SNP)

write_xlsx(
  list(
    fig5i_regional_points = source_data_fig5i_regional_points
  ),
  "results/source_data/source_data_fig5i.xlsx"
)

