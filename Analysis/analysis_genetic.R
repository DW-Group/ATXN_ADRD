library(dplyr)
library(data.table)
library(readxl)
library(writexl)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

rm(list=ls())

##---------------------------------------------------------
## Load dictionary
##---------------------------------------------------------

lod_dict <- read_excel("data/LOD_dictionary.xlsx")
marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx") %>% 
  left_join(lod_dict, by="targetName")
qc_prot <- unique(marker_dict_final %>% filter(!LOD_flag_all) %>% pull(targetName))

##---------------------------------------------------------
## Subset to cis and trans loci
##---------------------------------------------------------

### cis proteins with Bonferroni-significant cis pQTL

cis_res_dat <- fread("results/pqtl_cis_res.csv")

cis_prot <- cis_res_dat %>%
  filter(sig_Bonf == TRUE) %>%
  pull(Trait) %>%
  unique()

length(cis_prot)

### All genome-wide significant pQTL hits, annotated as cis/trans

gwas_hits_all <- fread("results/pqtl_allsig_res.csv") %>% 
  filter(Trait %in% qc_prot,
         Pval < 5e-8) %>% 
  left_join(marker_dict_final, by = c("Trait" = "targetName"), suffix = c("", "_dict")) %>% 
  mutate(type = ifelse(chr == Chr & Pos > cis_start & Pos < cis_end, "cis", "trans"))

### trans proteins

trans_prot <- gwas_hits_all %>%
  filter(type == "trans") %>%
  pull(Trait) %>%
  unique()

length(trans_prot)

### Proteins with any pQTL evidence

all_pqtl_prot <- unique(c(cis_prot, trans_prot))

length(all_pqtl_prot)
length(intersect(cis_prot, trans_prot))

### cis loci

cis_loci <- tibble(
  protein = cis_prot
) %>%
  left_join(
    marker_dict_final,
    by = c("protein" = "targetName")
  ) %>%
  transmute(
    protein,
    pqtl_class = "cis",
    protein_signal_group = case_when(
      protein %in% cis_prot & protein %in% trans_prot ~ "cis and trans",
      protein %in% cis_prot ~ "cis only",
      protein %in% trans_prot ~ "trans only",
      TRUE ~ "other"
    ),
    locus_id = paste0(protein, "_cis"),
    region_file = paste0(
      "data/colocalization/cis/NULISA_meta/coloc_NULISA_",
      targetName_unique,
      ".tsv"
    ),
    protein_gene = geneName
  )

### trans loci

sig_files_existing <- list.files(
  "data/colocalization/sig",
  pattern = "^coloc_NULISA_sig_.*\\.tsv$",
  full.names = FALSE
)

sig_loci_existing <- tibble(
  file_name = sig_files_existing,
  region_file = file.path("data/colocalization/sig", sig_files_existing),
  prot_rsid = gsub("^coloc_NULISA_sig_|\\.tsv$", "", sig_files_existing)
) %>%
  tidyr::separate(
    prot_rsid,
    into = c("protein", "lead_rsid"),
    sep = "_rs",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(
    lead_rsid = paste0("rs", lead_rsid)
  )

sig_loci_anno <- sig_loci_existing %>%
  left_join(
    gwas_hits_all %>%
      select(
        protein = Trait,
        lead_rsid = rsid,
        Chr,
        Pos,
        Pval,
        Effect,
        type,
        targetName_unique
      ) %>%
      distinct(protein, lead_rsid, .keep_all = TRUE),
    by = c("protein", "lead_rsid")
  )

trans_loci <- sig_loci_anno %>%
  filter(type == "trans") %>%
  left_join(
    marker_dict_final,
    by = c("protein" = "targetName")
  ) %>%
  transmute(
    protein,
    pqtl_class = "trans",
    protein_signal_group = case_when(
      protein %in% cis_prot & protein %in% trans_prot ~ "cis and trans",
      protein %in% cis_prot ~ "cis only",
      protein %in% trans_prot ~ "trans only",
      TRUE ~ "other"
    ),
    locus_id = paste0(protein, "_trans_", lead_rsid),
    lead_rsid,
    lead_chr = Chr,
    lead_pos = Pos,
    lead_p = Pval,
    lead_beta = Effect,
    region_file,
    protein_gene = geneName
  )

### Combine query loci

query_loci <- bind_rows(cis_loci, trans_loci) %>%
  distinct(locus_id, .keep_all = TRUE)

dim(query_loci)[1]
# 41

fwrite(query_loci, "results/tissue/pqtl_query_loci_from_existing_coloc_files.csv")

##---------------------------------------------------------
## Build pQTL query regions from existing coloc files
##---------------------------------------------------------

region_list <- list()

for (i in seq_len(nrow(query_loci))) {
  
  message("Reading region file ", i, "/", nrow(query_loci), ": ", query_loci$locus_id[i])
  
  tmp <- fread(query_loci$region_file[i]) %>%
    filter(
      !is.na(rsid),
      rsid != ".",
      !is.na(Chr),
      !is.na(Pos)
    )
  
  tmp_region <- tmp %>%
    summarise(
      region_chr_b37 = unique(Chr)[1],
      region_start_b37 = min(Pos, na.rm = TRUE),
      region_end_b37 = max(Pos, na.rm = TRUE),
      n_region_snps = n(),
      n_region_rsids = n_distinct(rsid)
    ) %>%
    mutate(
      protein = query_loci$protein[i],
      pqtl_class = query_loci$pqtl_class[i],
      protein_signal_group = query_loci$protein_signal_group[i],
      locus_id = query_loci$locus_id[i],
      region_file = query_loci$region_file[i],
      protein_gene = query_loci$protein_gene[i]
    )
  
  region_list[[i]] <- tmp_region
}

query_regions_b37 <- bind_rows(region_list) %>%
  mutate(
    region_chr_b37 = paste0("chr", region_chr_b37)
  ) %>%
  select(
    protein,
    pqtl_class,
    protein_signal_group,
    locus_id,
    protein_gene,
    region_file,
    region_chr_b37,
    region_start_b37,
    region_end_b37,
    n_region_snps,
    n_region_rsids
  )

fwrite(query_regions_b37, "results/tissue/pqtl_query_regions_b37.csv")

##---------------------------------------------------------
## Liftover b37 query regions to b38
##---------------------------------------------------------

query_regions <- fread("results/tissue/pqtl_query_regions_b37_b38.csv")

##-------------------------------
## GTEx directories
##-------------------------------

eqtl_dir <- file.path("gtex/v8/GTEx_Analysis_v8_eQTL")
sqtl_dir <- file.path("gtex/v8/GTEx_Analysis_v8_sQTL")

out_dir <- "results/gtex_v8_brain_qtl_region_scan"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

brain_tissues <- c(
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

##---------------------------------------------------------
## GTEx brain eQTL region scan
##---------------------------------------------------------

eqtl_pair_files <- list.files(
  eqtl_dir,
  pattern = "\\.v8\\.signif_variant_gene_pairs\\.txt\\.gz$",
  full.names = TRUE
)

eqtl_pair_files <- eqtl_pair_files[
  gsub("\\.v8\\.signif_variant_gene_pairs\\.txt\\.gz$", "", basename(eqtl_pair_files)) %in% brain_tissues
]

egene_files <- list.files(
  eqtl_dir,
  pattern = "\\.v8\\.egenes\\.txt\\.gz$",
  full.names = TRUE
)

egene_files <- egene_files[
  gsub("\\.v8\\.egenes\\.txt\\.gz$", "", basename(egene_files)) %in% brain_tissues
]

egene_anno_list <- list()

for (f in egene_files) {
  
  tissue <- gsub("\\.v8\\.egenes\\.txt\\.gz$", "", basename(f))
  
  tmp_anno <- fread(f, select = c("gene_id", "gene_name")) %>%
    distinct(gene_id, gene_name) %>%
    mutate(tissue = tissue)
  
  egene_anno_list[[length(egene_anno_list) + 1]] <- tmp_anno
}

egene_anno <- bind_rows(egene_anno_list)

eqtl_hits_list <- list()

query_dt <- as.data.table(query_regions)
query_dt[, region_start_b38 := as.numeric(region_start_b38)]
query_dt[, region_end_b38 := as.numeric(region_end_b38)]
setkey(query_dt, region_chr_b38, region_start_b38, region_end_b38)

for (f in eqtl_pair_files) {
  
  tissue <- gsub("\\.v8\\.signif_variant_gene_pairs\\.txt\\.gz$", "", basename(f))
  message("Scanning eQTL: ", tissue)
  
  qtl <- fread(f)
  
  qtl <- qtl %>%
    mutate(
      tissue = tissue,
      qtl_type = "eQTL",
      variant_chr_b38 = sub("^(chr[0-9XY]+)_.*$", "\\1", variant_id),
      variant_pos_b38 = as.numeric(sub("^chr[0-9XY]+_([0-9]+)_.*$", "\\1", variant_id))
    ) %>%
    left_join(
      egene_anno %>% filter(tissue == !!tissue),
      by = c("gene_id", "tissue")
    )
  
  qtl_dt <- as.data.table(qtl)
  qtl_dt[, qtl_start := variant_pos_b38]
  qtl_dt[, qtl_end := variant_pos_b38]
  setkey(qtl_dt, variant_chr_b38, qtl_start, qtl_end)
  
  hits <- foverlaps(
    qtl_dt,
    query_dt,
    by.x = c("variant_chr_b38", "qtl_start", "qtl_end"),
    by.y = c("region_chr_b38", "region_start_b38", "region_end_b38"),
    type = "within",
    nomatch = 0
  )
  
  if (nrow(hits) > 0) {
    eqtl_hits_list[[length(eqtl_hits_list) + 1]] <- as.data.frame(hits)
  }
}

eqtl_hits <- bind_rows(eqtl_hits_list) %>%
  mutate(
    qtl_type = "eQTL",
    is_target_gene = gene_name == protein_gene,
    neglog10_p = -log10(pval_nominal)
  )

fwrite(eqtl_hits, file.path(out_dir, "gtex_v8_brain_eqtl_region_hits_detail.csv"))

##---------------------------------------------------------
## GTEx brain sQTL region scan
##---------------------------------------------------------

sqtl_pair_files <- list.files(
  sqtl_dir,
  pattern = "\\.v8\\.sqtl_signifpairs\\.txt\\.gz$",
  full.names = TRUE
)

sqtl_pair_files <- sqtl_pair_files[
  gsub("\\.v8\\.sqtl_signifpairs\\.txt\\.gz$", "", basename(sqtl_pair_files)) %in% brain_tissues
]

sgene_files <- list.files(
  sqtl_dir,
  pattern = "\\.v8\\.sgenes\\.txt\\.gz$",
  full.names = TRUE
)

sgene_files <- sgene_files[
  gsub("\\.v8\\.sgenes\\.txt\\.gz$", "", basename(sgene_files)) %in% brain_tissues
]

sgene_anno_list <- list()

for (f in sgene_files) {
  
  tissue <- gsub("\\.v8\\.sgenes\\.txt\\.gz$", "", basename(f))
  
  tmp_anno <- fread(f, select = c("gene_id", "gene_name")) %>%
    distinct(gene_id, gene_name) %>%
    mutate(tissue = tissue)
  
  sgene_anno_list[[length(sgene_anno_list) + 1]] <- tmp_anno
}

sgene_anno <- bind_rows(sgene_anno_list)

sqtl_hits_list <- list()

for (f in sqtl_pair_files) {
  
  tissue <- gsub("\\.v8\\.sqtl_signifpairs\\.txt\\.gz$", "", basename(f))
  message("Scanning sQTL: ", tissue)
  
  qtl <- fread(f)
  
  qtl <- qtl %>%
    mutate(
      gene_id = sub(".*:", "", phenotype_id),
      tissue = tissue,
      qtl_type = "sQTL",
      variant_chr_b38 = sub("^(chr[0-9XY]+)_.*$", "\\1", variant_id),
      variant_pos_b38 = as.numeric(sub("^chr[0-9XY]+_([0-9]+)_.*$", "\\1", variant_id))
    ) %>%
    left_join(
      sgene_anno %>% filter(tissue == !!tissue),
      by = c("gene_id", "tissue")
    )
  
  qtl_dt <- as.data.table(qtl)
  qtl_dt[, qtl_start := variant_pos_b38]
  qtl_dt[, qtl_end := variant_pos_b38]
  setkey(qtl_dt, variant_chr_b38, qtl_start, qtl_end)
  
  hits <- foverlaps(
    qtl_dt,
    query_dt,
    by.x = c("variant_chr_b38", "qtl_start", "qtl_end"),
    by.y = c("region_chr_b38", "region_start_b38", "region_end_b38"),
    type = "within",
    nomatch = 0
  )
  
  if (nrow(hits) > 0) {
    sqtl_hits_list[[length(sqtl_hits_list) + 1]] <- as.data.frame(hits)
  }
}

sqtl_hits <- bind_rows(sqtl_hits_list) %>%
  mutate(
    qtl_type = "sQTL",
    is_target_gene = gene_name == protein_gene,
    neglog10_p = -log10(pval_nominal)
  )

fwrite(sqtl_hits, file.path(out_dir, "gtex_v8_brain_sqtl_region_hits_detail.csv"))

##---------------------------------------------------------
## Combine eQTL and sQTL
##---------------------------------------------------------

qtl_hits <- bind_rows(eqtl_hits, sqtl_hits) %>%
  mutate(
    qtl_gene_class = case_when(
      is_target_gene ~ "protein-encoding gene",
      !is.na(gene_name) ~ "nearby gene",
      TRUE ~ "unknown"
    )
  )

fwrite(qtl_hits, file.path(out_dir, "gtex_v8_brain_eqtl_sqtl_region_hits_detail.csv"))

## locus-level summary

qtl_summary <- qtl_hits %>%
  group_by(
    protein,
    pqtl_class,
    protein_signal_group,
    locus_id,
    protein_gene,
    qtl_type
  ) %>%
  summarise(
    has_brain_qtl_any = TRUE,
    has_brain_qtl_target_gene = any(is_target_gene, na.rm = TRUE),
    n_brain_tissues = n_distinct(tissue),
    n_qtl_genes = n_distinct(gene_name),
    min_p = min(pval_nominal, na.rm = TRUE),
    max_neglog10_p = max(neglog10_p, na.rm = TRUE),
    top_gene = gene_name[which.min(pval_nominal)][1],
    top_tissue = tissue[which.min(pval_nominal)][1],
    top_variant = variant_id[which.min(pval_nominal)][1],
    all_genes = paste(sort(unique(na.omit(gene_name))), collapse = "; "),
    all_tissues = paste(sort(unique(na.omit(tissue))), collapse = "; "),
    .groups = "drop"
  )

qtl_grid <- expand.grid(
  locus_id = query_regions$locus_id,
  qtl_type = c("eQTL", "sQTL"),
  stringsAsFactors = FALSE
) %>%
  left_join(
    query_regions %>%
      select(
        protein,
        pqtl_class,
        protein_signal_group,
        locus_id,
        protein_gene
      ),
    by = "locus_id"
  )

qtl_summary_full <- qtl_grid %>%
  left_join(
    qtl_summary,
    by = c(
      "protein",
      "pqtl_class",
      "protein_signal_group",
      "locus_id",
      "protein_gene",
      "qtl_type"
    )
  ) %>%
  mutate(
    has_brain_qtl_any = ifelse(is.na(has_brain_qtl_any), FALSE, has_brain_qtl_any),
    has_brain_qtl_target_gene = ifelse(is.na(has_brain_qtl_target_gene), FALSE, has_brain_qtl_target_gene),
    n_brain_tissues = ifelse(is.na(n_brain_tissues), 0, n_brain_tissues),
    n_qtl_genes = ifelse(is.na(n_qtl_genes), 0, n_qtl_genes)
  )

fwrite(qtl_summary_full, file.path(out_dir, "gtex_v8_brain_eqtl_sqtl_region_hits_summary.csv"))

## protein-level summary

protein_summary <- qtl_summary_full %>%
  group_by(protein, protein_signal_group, qtl_type) %>%
  summarise(
    n_loci = n_distinct(locus_id),
    n_loci_with_brain_qtl = sum(has_brain_qtl_any, na.rm = TRUE),
    n_loci_with_target_gene_brain_qtl = sum(has_brain_qtl_target_gene, na.rm = TRUE),
    min_p = ifelse(all(is.na(min_p)), NA_real_, min(min_p, na.rm = TRUE)),
    top_gene = ifelse(all(is.na(min_p)), NA_character_, top_gene[which.min(min_p)][1]),
    top_tissue = ifelse(all(is.na(min_p)), NA_character_, top_tissue[which.min(min_p)][1]),
    .groups = "drop"
  ) %>%
  mutate(
    has_any_brain_qtl = n_loci_with_brain_qtl > 0,
    has_any_target_gene_brain_qtl = n_loci_with_target_gene_brain_qtl > 0
  )

fwrite(protein_summary, file.path(out_dir, "gtex_v8_brain_eqtl_sqtl_protein_summary.csv"))

##---------------------------------------------------------
## KLKB1/F11 hotspot input files
##---------------------------------------------------------

hotspot_dir <- "results/hotspot_KLKB1_F11"
dir.create(hotspot_dir, showWarnings = FALSE, recursive = TRUE)

##-------------------------------
## Define KLKB1/F11 trans-hotspot proteins
##-------------------------------

klkb1_hotspot_proteins <- c(
  "NPY",
  "POSTN",
  "SFTPD",
  "TAFA5",
  "VEGFD"
)

klkb1_trans_loci <- trans_loci %>%
  filter(protein %in% klkb1_hotspot_proteins) %>%
  arrange(lead_chr, lead_pos, protein)

fwrite(klkb1_trans_loci, file.path(hotspot_dir, "KLKB1_F11_hotspot_trans_loci.csv"))

klkb1_trans_loci %>%
  select(
    protein,
    protein_gene,
    lead_rsid,
    lead_chr,
    lead_pos,
    lead_beta,
    lead_p,
    region_file
  )

##-------------------------------
## GTEx brain regulatory annotation for KLKB1/F11 hotspot
##-------------------------------

qtl_hits <- fread("results/gtex_v8_brain_qtl_region_scan/gtex_v8_brain_eqtl_sqtl_region_hits_detail.csv")

hotspot_qtl <- qtl_hits %>%
  filter(locus_id %in% klkb1_trans_loci$locus_id) %>%
  mutate(neglog10_p = -log10(pval_nominal))

fwrite(hotspot_qtl, file.path(hotspot_dir, "KLKB1_F11_hotspot_GTEx_brain_QTL_detail.csv"))

hotspot_qtl_summary_gene <- hotspot_qtl %>%
  group_by(qtl_type, gene_name) %>%
  summarise(
    n_tissues = n_distinct(tissue),
    min_p = min(pval_nominal, na.rm = TRUE),
    max_neglog10_p = max(neglog10_p, na.rm = TRUE),
    top_tissue = tissue[which.min(pval_nominal)][1],
    top_variant = variant_id[which.min(pval_nominal)][1],
    .groups = "drop"
  ) %>%
  arrange(qtl_type, min_p)

fwrite(hotspot_qtl_summary_gene, file.path(hotspot_dir, "KLKB1_F11_hotspot_GTEx_brain_QTL_gene_summary.csv"))

hotspot_qtl_summary_gene %>% print(n = 50)

focus_hotspot_genes <- c(
  "KLKB1", "F11",
  "UFSP2", "FAM149A", "CYP4V2", "LRP2BP",
  "TLR3", "SORBS2", "FAT1", "SNX25",
  "ANKRD37"
)

hotspot_qtl_focus <- hotspot_qtl_summary_gene %>%
  filter(gene_name %in% focus_hotspot_genes)

fwrite(hotspot_qtl_focus, file.path(hotspot_dir, "KLKB1_F11_hotspot_GTEx_focus_gene_summary.csv"))

hotspot_qtl_focus

##-------------------------------
## Pathway/network
##-------------------------------

affected_proteins <- klkb1_trans_loci %>%
  distinct(protein, protein_gene) %>%
  arrange(protein)
fwrite(affected_proteins, file.path(hotspot_dir, "KLKB1_F11_hotspot_affected_proteins.csv"))

candidate_local_genes <- hotspot_qtl_summary_gene %>%
  filter(gene_name %in% focus_hotspot_genes) %>%
  distinct(gene_name) %>%
  bind_rows(
    data.frame(gene_name = c("KLKB1", "F11"))
  ) %>%
  distinct(gene_name) %>%
  arrange(gene_name)
fwrite(candidate_local_genes, file.path(hotspot_dir, "KLKB1_F11_hotspot_candidate_local_genes.csv"))

### Load input gene sets

affected_proteins <- fread(
  file.path(hotspot_dir, "KLKB1_F11_hotspot_affected_proteins.csv")
)

candidate_local_genes <- fread(
  file.path(hotspot_dir, "KLKB1_F11_hotspot_candidate_local_genes.csv")
)

affected_genes <- affected_proteins %>%
  filter(!is.na(protein_gene), protein_gene != "") %>%
  transmute(gene_symbol = protein_gene, gene_set = "trans-regulated protein") %>%
  distinct()

local_genes <- candidate_local_genes %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  transmute(gene_symbol = gene_name, gene_set = "local regulatory gene") %>%
  distinct()

module_genes <- bind_rows(local_genes, affected_genes) %>%
  distinct(gene_symbol, gene_set)

module_gene_symbols <- unique(module_genes$gene_symbol)
local_gene_symbols <- unique(local_genes$gene_symbol)
affected_gene_symbols <- unique(affected_genes$gene_symbol)

write_xlsx(
  list(
    local_genes = local_genes,
    affected_genes = affected_genes,
    module_genes = module_genes
  ),
  file.path(hotspot_dir, "KLKB1_F11_hotspot_gene_sets.xlsx")
)

### Map gene symbols to Entrez IDs

gene_map <- bitr(
  module_gene_symbols,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "ENSEMBL"),
  OrgDb = org.Hs.eg.db
)

local_map <- gene_map %>%
  filter(SYMBOL %in% local_gene_symbols)

affected_map <- gene_map %>%
  filter(SYMBOL %in% affected_gene_symbols)

module_entrez <- unique(gene_map$ENTREZID)
local_entrez <- unique(local_map$ENTREZID)
affected_entrez <- unique(affected_map$ENTREZID)

unmapped_genes <- setdiff(module_gene_symbols, gene_map$SYMBOL); unmapped_genes

### Reactome pathway enrichment

react_local <- enrichPathway(
  gene = local_entrez,
  organism = "human",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  readable = TRUE
)

react_affected <- enrichPathway(
  gene = affected_entrez,
  organism = "human",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  readable = TRUE
)

react_module <- enrichPathway(
  gene = module_entrez,
  organism = "human",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  readable = TRUE
)

react_local_df <- as.data.frame(react_local)
react_affected_df <- as.data.frame(react_affected)
react_module_df <- as.data.frame(react_module)

fwrite(react_local_df, file.path(hotspot_dir, "Reactome_local_regulatory_genes.csv"))
fwrite(react_affected_df, file.path(hotspot_dir, "Reactome_affected_proteins.csv"))
fwrite(react_module_df, file.path(hotspot_dir, "Reactome_combined_module.csv"))
