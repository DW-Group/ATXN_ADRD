library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(patchwork)
library(cowplot)
library(writexl)

load("data/colocalization/trans/pqtl_hotspot.RData")

### Stacked locus

plot_dat <- trans_all_anno %>%
  mutate(
    logp = -log10(Pval),
    pos_mb = Pos / 1e6,
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
  arrange(ld_bin)

xrange <- range(plot_dat$logp, na.rm = TRUE); xrange[2] <- xrange[2] + 0.5

ld_cols <- c(
  "0-0.2"   = "#EEF2F1",
  "0.2-0.4" = "#AFC6C4",
  "0.4-0.6" = "#8DAFB0",
  "0.6-0.8" = "#6F9397",
  "0.8-1.0" = "#4F6F78"
)

lead_snp <- "rs66530140"

label_df <- plot_dat %>%
  dplyr::filter(rsid == lead_snp) %>%
  mutate(
    p_lab = formatC(Pval, format = "e", digits = 2),
    p_lab = sub("e\\+?", "e", p_lab),
    label = paste0(protein, "\nP = ", p_lab)
  )

yrange_bp <- range(plot_dat$Pos, na.rm = TRUE)
yrange_mb <- yrange_bp / 1e6

prot_order <- c("NPY", "POSTN", "SFTPD", "TAFA5", "VEGFD")
prot_order <- prot_order[prot_order %in% unique(plot_dat$protein)]

make_one_panel <- function(prot, show_y = TRUE, show_x = TRUE, show_label = TRUE) {
  dat_i <- plot_dat %>% dplyr::filter(protein == prot)
  lab_i <- label_df %>% dplyr::filter(protein == prot)
  
  p_i <- ggplot(dat_i, aes(x = logp, y = pos_mb)) +
    geom_point(
      aes(color = ld_bin),
      size = 1.5,
      alpha = 0.9
    ) +
    scale_color_manual(values = ld_cols, name = expression(r^2)) +
    scale_x_continuous(
      name = expression(-log[10](P)),
      limits = xrange,
      breaks = function(x) {
        b <- pretty(x)
        b[b < max(b)]
      },
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "Chr 4 position (Mb)",
      limits = yrange_mb,
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.line.y = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(4, "pt"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 13),
      plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
      legend.position = "none",
      legend.title = element_text(size = 11.5),
      legend.text = element_text(size = 11.5)
    )
  
  if (show_label && nrow(lab_i) > 0) {
    p_i <- p_i +
      ggrepel::geom_text_repel(
        data = lab_i,
        aes(label = label),
        size = 4,
        color = "#4F6F78",
        box.padding = 0.5,
        point.padding = 0.6,
        force = 2,
        force_pull = 0.3,
        segment.color = "grey50"
      )
  }
  
  if (!show_y) {
    p_i <- p_i +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )
  }
  
  if (!show_x) {
    p_i <- p_i +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  p_i
}

p_list <- lapply(seq_along(prot_order), function(i) {
  make_one_panel(
    prot = prot_order[i],
    show_y = (i == 1),  
    show_x = TRUE,
    show_label = TRUE
  )
})

p_legend <- ggplot(
  plot_dat %>% dplyr::filter(protein == prot_order[1]),
  aes(x = logp, y = pos_mb, color = ld_bin)
) +
  geom_point(size = 1.5, alpha = 0.9) +
  scale_color_manual(values = ld_cols, name = expression(r^2)) +
  theme_void() +
  theme(legend.position = "left",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

legend <- patchwork::wrap_elements(
  cowplot::get_legend(p_legend)
)

p_main <- wrap_plots(plotlist = p_list, nrow = 1)

p_final <- (legend | p_main) + plot_layout(widths = c(1.2, 10))

ggsave(
  "figures/fig5_trans_hotspot.pdf",
  p_final,
  width = 11,
  height = 3.8
)

### Annotate gene

chr_use <- unique(plot_dat$Chr)
stopifnot(length(chr_use) == 1)

gene_dat <- plot_dat %>%
  dplyr::filter(protein == unique(protein)[1]) %>%
  mutate(
    Chr = as.integer(Chr),
    Pos = as.numeric(Pos),
    Pval = as.numeric(Pval),
    rsid = as.character(rsid)
  )

loc <- locus(
  data = as.data.frame(gene_dat),
  seqname = chr_use,
  xrange = yrange_bp,
  ens_db = "EnsDb.Hsapiens.v75",
  chrom = "Chr",
  pos = "Pos",
  p = "Pval",
  labs = "rsid"
)

p_gene <- gg_genetracks(
  loc,
  filter_gene_biotype = "protein_coding",
  xticks = TRUE,
  xlab = "Chr 4 position (Mb)",
  maxrows = 10,
  cex.text = 0.9,
  text_pos = "top",
  italics = FALSE,
  gene_col = scales::alpha("#7C6F78", 0.55),
  exon_col = "#4A3C43",
  exon_border = NA
) +
  scale_x_continuous(
    limits = yrange_mb,
    expand = c(0, 0)
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black", size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y = element_line(color = "black", linewidth = 0.4),
    
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12)
  )

ggsave(
  "figures/fig5_trans_hotspot_gene_track.pdf",
  p_gene,
  width = 3.5,
  height = 2.2
) 
