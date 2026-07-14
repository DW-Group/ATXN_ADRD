library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(cowplot)
library(writexl)

## ---------------------------------------------------------
## Setup
## ---------------------------------------------------------

out_dir <- "figures/fig5f_coloc_callout"
source_dir <- "results/source_data"
ld_dir <- "data/colocalization/ld_cache"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

ld_cols <- c(
  "0-0.2"   = "#E8E6E3",
  "0.2-0.4" = "#9FB9B9",
  "0.4-0.6" = "#6FA8A6",
  "0.6-0.8" = "#E6A86C",
  "0.8-1.0" = "#C65A5A"
)

fig5f_cfg <- tribble(
  ~protein,       ~outcome,     ~outcome_label,                         ~lead_snp,      ~qtl_mode, ~coloc_pattern,
  "BD_pTau_181",  "lee",        "Cognitive performance (Lee et al.)",    "rs242948",     "cis",    "pp4",
  "TREM2",        "bellenguez", "AD/ADRD (Bellenguez et al.)",           "rs72924659",   "trans",  "pp4",
  "SFTPD",        "bellenguez", "AD/ADRD (Bellenguez et al.)",           NA_character_,  "cis",    "pp3"
)

## ---------------------------------------------------------
## Cached LD helper
## ---------------------------------------------------------

ld_file <- function(lead_snp) {
  file.path(ld_dir, paste0("ld_", lead_snp, ".tsv"))
}

read_ld <- function(lead_snp) {
  f <- ld_file(lead_snp)
  if (!file.exists(f)) stop("Missing LD file: ", f)
  
  fread(f) %>%
    transmute(
      rsid = as.character(rsid),
      r2 = as.numeric(r2)
    ) %>%
    filter(!is.na(rsid), rsid != ".") %>%
    group_by(rsid) %>%
    summarise(r2 = max(r2, na.rm = TRUE), .groups = "drop")
}

## ---------------------------------------------------------
## Add LD bins for plotting
## ---------------------------------------------------------

add_ld <- function(dt) {
  dt %>%
    mutate(
      logp = -log10(p),
      pos_mb = pos_b37 / 1e6,
      ld_bin = case_when(
        is.na(r2) ~ "0-0.2",
        r2 >= 0.8 ~ "0.8-1.0",
        r2 >= 0.6 ~ "0.6-0.8",
        r2 >= 0.4 ~ "0.4-0.6",
        r2 >= 0.2 ~ "0.2-0.4",
        TRUE ~ "0-0.2"
      ),
      ld_bin = factor(
        ld_bin,
        levels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
      )
    ) %>%
    arrange(ld_bin, pos_b37, rsid)
}

## ---------------------------------------------------------
## Data reader
## ---------------------------------------------------------

read_callout_data <- function(protein, outcome, lead_snp = NA_character_, qtl_mode = "cis") {
  
  if (qtl_mode == "trans") {
    file_protein <- file.path(
      "data/colocalization",
      paste0("coloc_NULISA_sig_", protein, "_", lead_snp, ".tsv")
    )
  } else {
    file_protein <- file.path(
      "data/colocalization",
      paste0("coloc_NULISA_", protein, ".tsv")
    )
  }
  
  protein_dt <- fread(file_protein) %>%
    transmute(
      rsid = as.character(rsid),
      chr = as.integer(Chr),
      pos_b37 = as.numeric(Pos),
      p = as.numeric(Pval)
    ) %>%
    filter(
      !is.na(rsid),
      rsid != ".",
      !is.na(chr),
      !is.na(pos_b37),
      !is.na(p),
      p > 0,
      p <= 1
    )
  
  if (outcome == "lee") {
    
    if (qtl_mode == "trans") {
      file_outcome <- file.path(
        "data/colocalization",
        paste0("coloc_lee_sig_", protein, "_", lead_snp, ".tsv")
      )
    } else {
      file_outcome <- file.path(
        "data/colocalization",
        paste0("coloc_lee_", protein, "_harmonized.tsv")
      )
    }
    
    outcome_dt <- fread(file_outcome) %>%
      transmute(
        rsid = as.character(MarkerName),
        p = as.numeric(Pval)
      )
  }
  
  if (outcome == "bellenguez") {
    
    if (qtl_mode == "trans") {
      file_outcome <- file.path(
        "data/colocalization",
        paste0("coloc_bellenguez_sig_", protein, "_", lead_snp, ".tsv")
      )
    } else {
      file_outcome <- file.path(
        "data/colocalization",
        paste0("coloc_bellenguez_ad_", protein, "_harmonized.tsv")
      )
    }
    
    outcome_dt <- fread(file_outcome) %>%
      transmute(
        rsid = as.character(variant_id),
        p = as.numeric(p_value)
      )
  }
  
  outcome_dt <- outcome_dt %>%
    filter(
      !is.na(rsid),
      rsid != ".",
      !is.na(p),
      p > 0,
      p <= 1
    )
  
  overlap_rsid <- intersect(protein_dt$rsid, outcome_dt$rsid)
  if (length(overlap_rsid) == 0) stop("No overlapping variants for ", protein, " / ", outcome)
  
  variant_anno <- protein_dt %>%
    filter(rsid %in% overlap_rsid) %>%
    select(rsid, chr, pos_b37) %>%
    distinct()
  
  protein_dt <- protein_dt %>%
    filter(rsid %in% overlap_rsid) %>%
    select(rsid, chr, pos_b37, p)
  
  outcome_dt <- outcome_dt %>%
    filter(rsid %in% overlap_rsid) %>%
    select(rsid, p) %>%
    left_join(variant_anno, by = "rsid") %>%
    select(rsid, chr, pos_b37, p)
  
  list(
    protein_dt = protein_dt,
    outcome_dt = outcome_dt
  )
}

get_ld_anchors <- function(protein_dt, outcome_dt, lead_snp = NA_character_, coloc_pattern = "pp4") {
  
  if (coloc_pattern == "pp3") {
    ld_outcome <- outcome_dt %>%
      arrange(p) %>%
      slice(1) %>%
      pull(rsid)
    
    ld_protein <- protein_dt %>%
      arrange(p) %>%
      slice(1) %>%
      pull(rsid)
  } else {
    ld_outcome <- lead_snp
    ld_protein <- lead_snp
  }
  
  tibble(
    LD_anchor_outcome_GWAS = ld_outcome,
    LD_anchor_protein_pQTL = ld_protein
  )
}

## ---------------------------------------------------------
## Panel plot
## ---------------------------------------------------------

make_panel <- function(dt, title_label, label_df, xlim, ylim, show_x = FALSE) {
  
  p <- ggplot(dt, aes(x = pos_mb, y = logp)) +
    geom_point(
      data = dt %>% filter(ld_bin == "0-0.2"),
      aes(fill = ld_bin),
      shape = 21,
      size = 2.5,
      color = "gray70",
      stroke = 0.15,
      alpha = 0.7
    ) +
    geom_point(
      data = dt %>% filter(ld_bin != "0-0.2"),
      aes(fill = ld_bin),
      shape = 21,
      size = 2.5,
      color = "gray31",
      stroke = 0.25,
      alpha = 0.8
    ) +
    geom_point(
      data = label_df,
      aes(x = pos_mb, y = logp),
      shape = 23,
      fill = "#9A8FBF",
      color = "black",
      size = 2.8,
      stroke = 0.4
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = rsid),
      size = 4,
      color = "black",
      box.padding = 0.3,
      point.padding = 0.2,
      force = 1.2,
      segment.color = "grey50"
    ) +
    scale_fill_manual(values = ld_cols) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(
      title = title_label,
      x = "Position (Mb)",
      y = expression(-log[10](P))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 12),
      plot.title = element_text(color = "black", size = 12, hjust = 0),
      legend.position = "none"
    )
  
  if (!show_x) {
    p <- p +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  p
}

## ---------------------------------------------------------
## Generate one callout panel
## ---------------------------------------------------------

make_callout <- function(protein, outcome, outcome_label,
                         lead_snp = NA_character_,
                         qtl_mode = "cis",
                         coloc_pattern = "pp4") {
  
  dat_list <- read_callout_data(
    protein = protein,
    outcome = outcome,
    lead_snp = lead_snp,
    qtl_mode = qtl_mode
  )
  
  protein_dt <- dat_list$protein_dt
  outcome_dt <- dat_list$outcome_dt
  
  chr_use <- unique(protein_dt$chr)
  if (length(chr_use) != 1) stop("More than one chromosome found for ", protein)
  
  anchor_df <- get_ld_anchors(
    protein_dt = protein_dt,
    outcome_dt = outcome_dt,
    lead_snp = lead_snp,
    coloc_pattern = coloc_pattern
  )
  
  ld_outcome <- anchor_df$LD_anchor_outcome_GWAS
  ld_protein <- anchor_df$LD_anchor_protein_pQTL
  
  outcome_dt <- outcome_dt %>%
    left_join(read_ld(ld_outcome), by = "rsid") %>%
    add_ld()
  
  protein_dt <- protein_dt %>%
    left_join(read_ld(ld_protein), by = "rsid") %>%
    add_ld()
  
  xlim <- range(protein_dt$pos_b37, na.rm = TRUE) / 1e6
  ylim_outcome <- c(0, max(outcome_dt$logp, na.rm = TRUE) * 1.05)
  ylim_protein <- c(0, max(protein_dt$logp, na.rm = TRUE) * 1.08)
  
  label_outcome <- outcome_dt %>%
    filter(rsid == ld_outcome)
  
  label_protein <- protein_dt %>%
    filter(rsid == ld_protein)
  
  p_outcome <- make_panel(
    dt = outcome_dt,
    title_label = outcome_label,
    label_df = label_outcome,
    xlim = xlim,
    ylim = ylim_outcome,
    show_x = FALSE
  )
  
  p_protein <- make_panel(
    dt = protein_dt,
    title_label = protein,
    label_df = label_protein,
    xlim = xlim,
    ylim = ylim_protein,
    show_x = FALSE
  )
  
  loc <- locus(
    data = protein_dt,
    seqname = chr_use,
    xrange = range(protein_dt$pos_b37, na.rm = TRUE),
    ens_db = "EnsDb.Hsapiens.v75",
    chrom = "chr",
    pos = "pos_b37",
    p = "p",
    labs = "rsid"
  )
  
  p_gene <- gg_genetracks(
    loc,
    filter_gene_biotype = "protein_coding",
    xticks = TRUE,
    xlab = paste0("Chr ", chr_use, " position (Mb)"),
    gene_col = "#7C6F78",
    exon_col = "#4A3C43",
    cex.text = 0.8
  ) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.text.x = element_text(color = "black", size = 11),
      axis.title.x = element_text(color = "black", size = 12)
    )
  
  p_final <- p_outcome / p_protein / p_gene +
    plot_layout(heights = c(0.6, 0.6, 1))
  
  fig_file <- file.path(
    out_dir,
    paste0("fig5f_", protein, "_", outcome, ".pdf")
  )
  
  ggsave(
    filename = fig_file,
    plot = p_final,
    width = 3.5,
    height = 9.2
  )
  
  source_points <- bind_rows(
    outcome_dt %>%
      mutate(
        protein = protein,
        outcome = outcome,
        outcome_label = outcome_label,
        track = "Outcome GWAS",
        lead_snp = ifelse(is.na(lead_snp), "", lead_snp),
        LD_anchor_snp = ld_outcome,
        is_labeled_snp = rsid == ld_outcome
      ),
    protein_dt %>%
      mutate(
        protein = protein,
        outcome = outcome,
        outcome_label = outcome_label,
        track = "Protein pQTL",
        lead_snp = ifelse(is.na(lead_snp), "", lead_snp),
        LD_anchor_snp = ld_protein,
        is_labeled_snp = rsid == ld_protein
      )
  ) %>%
    transmute(
      protein,
      outcome,
      outcome_label,
      track,
      chr,
      position = pos_b37,
      position_Mb = pos_mb,
      rsid,
      p_value = p,
      neg_log10_p = logp,
      lead_snp,
      LD_anchor_snp,
      LD_r2_to_anchor = r2,
      LD_bin = as.character(ld_bin),
      is_labeled_snp
    ) %>%
    arrange(protein, outcome, track, position, rsid)
  
  list(
    points = source_points
  )
}

## ---------------------------------------------------------
## Run plotting
## ---------------------------------------------------------

callout_res <- lapply(seq_len(nrow(fig5f_cfg)), function(i) {
  make_callout(
    protein = fig5f_cfg$protein[i],
    outcome = fig5f_cfg$outcome[i],
    outcome_label = fig5f_cfg$outcome_label[i],
    lead_snp = fig5f_cfg$lead_snp[i],
    qtl_mode = fig5f_cfg$qtl_mode[i],
    coloc_pattern = fig5f_cfg$coloc_pattern[i]
  )
})

source_data_fig5f_points <- bind_rows(lapply(callout_res, `[[`, "points"))

## ---------------------------------------------------------
## LD legend
## ---------------------------------------------------------

p_legend <- ggplot(
  data.frame(
    x = 1,
    y = 1,
    ld_bin = factor(names(ld_cols), levels = names(ld_cols))
  ),
  aes(x = x, y = y, fill = ld_bin)
) +
  geom_point(
    shape = 21,
    size = 3,
    color = "gray31",
    stroke = 0.3
  ) +
  scale_fill_manual(
    values = ld_cols,
    name = expression(r^2)
  ) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

legend <- cowplot::get_legend(p_legend)

pdf(file.path(out_dir, "fig5f_ld_legend.pdf"), width = 2.5, height = 3)
grid::grid.draw(legend)
dev.off()

## ---------------------------------------------------------
## Source data
## ---------------------------------------------------------

write_xlsx(
  list(
    fig5f_points = source_data_fig5f_points
  ),
  file.path(source_dir, "source_data_fig5f.xlsx")
)