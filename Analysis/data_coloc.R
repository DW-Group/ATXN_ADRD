library(dplyr)
library(data.table)
library(readxl)
library(stringr)

rm(list = ls())

##-------------------------------
## Load pQTL summary statistics
##-------------------------------

pqtl_dat <- fread("results/pqtl_meta_all_variants.csv")

marker_dict <- read_excel("data/NULISA_pQTL_dictionary.xlsx") %>%
  select(
    targetName,
    gene,
    gene_chr,
    gene_start,
    gene_end
  ) %>%
  mutate(
    cis_start = pmax(gene_start - 1e6, 0),
    cis_end = gene_end + 1e6
  )

pqtl_dat <- pqtl_dat %>%
  left_join(marker_dict, by = c("Trait" = "targetName")) %>%
  mutate(
    cis_region = Chr == gene_chr & Pos >= cis_start & Pos <= cis_end,
    pval = as.numeric(Pval),
    se = (UCI_Effect - LCI_Effect) / (2 * 1.96)
  )

##-------------------------------
## Load independent trans lead pQTLs
##-------------------------------

trans_leads <- read_excel("results/pqtl_trans_leads_ldclumped_r2_0.1.xlsx") %>%
  select(
    Trait,
    lead_rsid = rsid,
    lead_chr = Chr,
    lead_pos = Pos
  ) %>%
  mutate(
    region_start = pmax(lead_pos - 1e6, 0),
    region_end = lead_pos + 1e6
  )

##-------------------------------
## Load outcome GWAS summary statistics
##-------------------------------

bellenguez <- fread("data/gwas_outcomes/bellenguez_ad.tsv.gz") %>%
  transmute(
    outcome = "bellenguez_ad",
    rsid = variant_id,
    Chr = chromosome,
    Pos = base_pair_location,
    beta_outcome = beta,
    se_outcome = standard_error,
    p_outcome = as.numeric(p_value),
    effect_allele_outcome = effect_allele,
    other_allele_outcome = other_allele,
    eaf_outcome = effect_allele_frequency,
    n_outcome = n_cases + n_controls
  )

wightman <- fread("data/gwas_outcomes/wightman_ad.txt") %>%
  transmute(
    outcome = "wightman_ad",
    rsid = rsid,
    Chr = chromosome,
    Pos = base_pair_location,
    beta_outcome = beta,
    se_outcome = standard_error,
    p_outcome = as.numeric(p_value),
    effect_allele_outcome = effect_allele,
    other_allele_outcome = other_allele,
    eaf_outcome = effect_allele_frequency,
    n_outcome = N
  )

wightmanukb <- fread("data/gwas_outcomes/wightmanukb_ad.txt") %>%
  filter(
    !is.na(rsid),
    rsid != ".",
    is.finite(z),
    is.finite(N)
  ) %>%
  transmute(
    outcome = "wightmanukb_ad",
    rsid = rsid,
    Chr = chr,
    Pos = PosGRCh37,
    beta_outcome = z / sqrt(N),
    se_outcome = 1 / sqrt(N),
    p_outcome = as.numeric(p),
    effect_allele_outcome = testedAllele,
    other_allele_outcome = otherAllele,
    eaf_outcome = NA_real_,
    n_outcome = N
  )

lee <- fread("data/gwas_outcomes/lee_cognitive_performance.txt") %>%
  transmute(
    outcome = "lee_cp",
    rsid = MarkerName,
    Chr = CHR,
    Pos = POS,
    beta_outcome = Beta,
    se_outcome = SE,
    p_outcome = as.numeric(Pval),
    effect_allele_outcome = A1,
    other_allele_outcome = A2,
    eaf_outcome = EAF,
    n_outcome = N
  )

outcome_dat <- bind_rows(
  bellenguez,
  wightman,
  wightmanukb,
  lee
) %>%
  filter(
    !is.na(rsid),
    !is.na(Chr),
    !is.na(Pos),
    !is.na(beta_outcome),
    !is.na(se_outcome),
    !is.na(p_outcome)
  )

##-------------------------------
## Format coloc input
##-------------------------------

format_coloc_region <- function(pqtl_sub, outcome_sub) {
  
  pqtl_fmt <- pqtl_sub %>%
    transmute(
      rsid,
      Chr,
      Pos,
      beta_pqtl = Effect,
      se_pqtl = se,
      p_pqtl = pval,
      effect_allele_pqtl = Alt,
      other_allele_pqtl = Ref,
      eaf_pqtl = if ("AAF" %in% names(pqtl_sub)) AAF else NA_real_,
      n_pqtl = if ("N" %in% names(pqtl_sub)) N else NA_real_
    )
  
  outcome_fmt <- outcome_sub %>%
    select(
      outcome,
      rsid,
      beta_outcome,
      se_outcome,
      p_outcome,
      effect_allele_outcome,
      other_allele_outcome,
      eaf_outcome,
      n_outcome
    )
  
  pqtl_fmt %>%
    inner_join(outcome_fmt, by = "rsid") %>%
    filter(
      !is.na(beta_pqtl),
      !is.na(se_pqtl),
      !is.na(p_pqtl),
      !is.na(beta_outcome),
      !is.na(se_outcome),
      !is.na(p_outcome)
    )
}

##-------------------------------
## Cis-region coloc inputs
##-------------------------------

cis_query_loci <- marker_dict %>%
  transmute(
    Trait = targetName,
    gene,
    Chr = gene_chr,
    region_start = cis_start,
    region_end = cis_end,
    pqtl_class = "cis"
  )

cis_file_list <- list()
k <- 1

for (i in seq_len(nrow(cis_query_loci))) {
  
  this_locus <- cis_query_loci[i, ]
  
  pqtl_sub <- pqtl_dat %>%
    filter(
      Trait == this_locus$Trait,
      Chr == this_locus$Chr,
      Pos >= this_locus$region_start,
      Pos <= this_locus$region_end
    )
  
  if (nrow(pqtl_sub) == 0) next
  
  for (outcome_name in unique(outcome_dat$outcome)) {
    
    outcome_sub <- outcome_dat %>%
      filter(
        outcome == outcome_name,
        Chr == this_locus$Chr,
        Pos >= this_locus$region_start,
        Pos <= this_locus$region_end
      )
    
    coloc_dat <- format_coloc_region(pqtl_sub, outcome_sub)
    
    if (nrow(coloc_dat) < 2) next
    
    out_file <- file.path(
      "data/colocalization/coloc_cis",
      paste0(this_locus$Trait, "_", outcome_name, "_cis.csv")
    )
    
    fwrite(coloc_dat, out_file)
    
    cis_file_list[[k]] <- tibble(
      Trait = this_locus$Trait,
      outcome = outcome_name,
      pqtl_class = "cis",
      Chr = this_locus$Chr,
      region_start = this_locus$region_start,
      region_end = this_locus$region_end,
      region_file = out_file,
      n_variants = nrow(coloc_dat)
    )
    
    k <- k + 1
  }
}

cis_file_index <- bind_rows(cis_file_list)

fwrite(
  cis_file_index,
  "data/colocalization/coloc_cis_file_index.csv"
)

##-------------------------------
## Trans-region coloc inputs
##-------------------------------

trans_file_list <- list()
k <- 1

for (i in seq_len(nrow(trans_leads))) {
  
  this_locus <- trans_leads[i, ]
  
  pqtl_sub <- pqtl_dat %>%
    filter(
      Trait == this_locus$Trait,
      Chr == this_locus$lead_chr,
      Pos >= this_locus$region_start,
      Pos <= this_locus$region_end
    )
  
  if (nrow(pqtl_sub) == 0) next
  
  for (outcome_name in unique(outcome_dat$outcome)) {
    
    outcome_sub <- outcome_dat %>%
      filter(
        outcome == outcome_name,
        Chr == this_locus$lead_chr,
        Pos >= this_locus$region_start,
        Pos <= this_locus$region_end
      )
    
    coloc_dat <- format_coloc_region(pqtl_sub, outcome_sub)
    
    if (nrow(coloc_dat) < 2) next
    
    out_file <- file.path(
      "data/colocalization/coloc_sig",
      paste0(
        this_locus$Trait,
        "_",
        this_locus$lead_rsid,
        "_",
        outcome_name,
        "_trans.csv"
      )
    )
    
    fwrite(coloc_dat, out_file)
    
    trans_file_list[[k]] <- tibble(
      Trait = this_locus$Trait,
      outcome = outcome_name,
      pqtl_class = "trans",
      lead_rsid = this_locus$lead_rsid,
      Chr = this_locus$lead_chr,
      region_start = this_locus$region_start,
      region_end = this_locus$region_end,
      region_file = out_file,
      n_variants = nrow(coloc_dat)
    )
    
    k <- k + 1
  }
}

trans_file_index <- bind_rows(trans_file_list)

fwrite(
  trans_file_index,
  "data/colocalization/coloc_sig_file_index.csv"
)

##-------------------------------
## Save query index
##-------------------------------

coloc_file_index <- bind_rows(
  cis_file_index,
  trans_file_index
)

fwrite(
  coloc_file_index,
  "data/colocalization/coloc_file_index.csv"
)