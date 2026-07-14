library(data.table)
library(dplyr)
library(readxl)
library(stringr)
library(TwoSampleMR)

rm(list = ls())

##-----------------------------------------
## Load data
##-----------------------------------------

lod_dict <- read_excel("data/LOD_dictionary.xlsx")

marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx") %>%
  left_join(lod_dict, by = "targetName")

cis_res_dat <- fread("results/pqtl_cis_res.csv")

##-----------------------------------------
## Select cis instruments
##-----------------------------------------

cis_res_sub <- cis_res_dat %>%
  mutate(
    SE = (UCI_Effect - LCI_Effect) / (2 * 1.96),
    F_stats = (Effect / SE)^2,
    INFO_list = str_extract(Info, "(?<=INFO_LIST=)[0-9\\.,]+"),
    INFO_min = sapply(
      strsplit(INFO_list, ","),
      function(x) min(as.numeric(x))
    )
  ) %>%
  filter(
    !(Trait %in% (
      marker_dict_final %>%
        filter(LOD_flag_all) %>%
        pull(targetName)
    )),
    INFO_min >= 0.8
  ) %>%
  group_by(Trait) %>%
  mutate(
    Bonf = 0.05 / n(),
    sig_Bonf = Pval <= Bonf
  ) %>%
  ungroup() %>%
  filter(
    sig_Bonf,
    F_stats > 10
  )

##-----------------------------------------
## Retain independent cis instruments
##-----------------------------------------

clumped_res_dat <- read_excel("results/mr/ld_clumped_instruments_r2_0.001.xlsx")

cis_res_final <- cis_res_sub %>%
  filter(
    paste0(Trait, "_", rsid) %in%
      (
        clumped_res_dat %>%
          mutate(prot_rsid = paste0(Trait, "_", RS_Number)) %>%
          pull(prot_rsid)
      )
  )

write.csv(
  cis_res_final,
  "results/mr/mr_instruments.csv",
  row.names = FALSE
)

##-----------------------------------------
## Format exposure data
##-----------------------------------------

cis_res_final_fmt <- cis_res_final %>%
  mutate(
    SNP = rsid,
    exposure = Trait,
    beta.exposure = Effect,
    se.exposure = SE,
    effect_allele.exposure = Alt,
    other_allele.exposure = Ref,
    pval.exposure = Pval,
    chr.exposure = Chr,
    position = Pos,
    eaf.exposure = if ("AAF" %in% names(.)) AAF else NA_real_,
    samplesize.exposure = if ("Num_Cases" %in% names(.)) Num_Cases else NA_real_,
    gene.exposure = if ("gene" %in% names(.)) gene else NA_character_
  ) %>%
  select(
    SNP,
    exposure,
    beta.exposure,
    se.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    eaf.exposure,
    pval.exposure,
    gene.exposure,
    samplesize.exposure,
    chr.exposure,
    position
  )

pqtl_pheno_list <- unique(cis_res_final_fmt$exposure)

instru_pqtl_list <- list()

for (pheno in pqtl_pheno_list) {
  instru_pqtl_list[[pheno]] <- cis_res_final_fmt %>%
    filter(exposure == pheno) %>%
    mutate(id.exposure = pheno)
}

##-----------------------------------------
## Load outcome GWAS summary statistics
##-----------------------------------------

wightman <- fread("data/gwas_outcomes/wightman_ad.txt") %>%
  mutate(outcome = "wightman_ad")

wightmanukb <- fread("data/gwas_outcomes/wightmanukb_ad.txt") %>%
  mutate(outcome = "wightmanukb_ad")

bellenguez <- fread("data/gwas_outcomes/bellenguez_ad.tsv.gz") %>%
  mutate(outcome = "bellenguez_ad")

lee <- fread("data/gwas_outcomes/lee_cognitive_performance.txt") %>%
  mutate(outcome = "lee_cp")

##-----------------------------------------
## Subset outcome data to MR instruments
##-----------------------------------------

wightman_sub <- wightman %>%
  filter(rsid %in% cis_res_final$rsid) %>%
  mutate(
    ncase_approx = (N - sqrt(N * (N - Neffective))) / 2,
    ncontrol_approx = N - ncase_approx
  )

wightmanukb_sub <- wightmanukb %>%
  filter(
    rsid %in% cis_res_final$rsid,
    is.finite(z),
    is.finite(N),
    !is.na(rsid),
    rsid != "."
  ) %>%
  mutate(
    beta_approx = z / sqrt(N),
    se_approx = 1 / sqrt(N),
    ncase_approx = N * (86531 / 762917),
    ncontrol_approx = N - ncase_approx
  )

bellenguez_sub <- bellenguez %>%
  filter(variant_id %in% cis_res_final$rsid) %>%
  mutate(
    n = n_cases + n_controls,
    p_value = as.numeric(p_value)
  ) %>%
  select(-base_pair_location) %>%
  left_join(
    cis_res_final %>% select(rsid, Pos),
    by = c("variant_id" = "rsid")
  )

lee_sub <- lee %>%
  filter(MarkerName %in% cis_res_final$rsid)

##-----------------------------------------
## Format outcome data
##-----------------------------------------

out_wightman <- format_data(
  as.data.frame(wightman_sub),
  type = "outcome",
  phenotype_col = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  ncase_col = "ncase_approx",
  ncontrol_col = "ncontrol_approx",
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)

out_wightmanukb <- format_data(
  as.data.frame(wightmanukb_sub),
  type = "outcome",
  phenotype_col = "outcome",
  snp_col = "rsid",
  beta_col = "beta_approx",
  se_col = "se_approx",
  effect_allele_col = "testedAllele",
  other_allele_col = "otherAllele",
  pval_col = "p",
  ncase_col = "ncase_approx",
  ncontrol_col = "ncontrol_approx",
  chr_col = "chr",
  pos_col = "PosGRCh37"
)

out_bellenguez <- format_data(
  as.data.frame(bellenguez_sub),
  type = "outcome",
  phenotype_col = "outcome",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  samplesize_col = "n",
  chr_col = "chromosome",
  pos_col = "Pos"
)

out_lee <- format_data(
  as.data.frame(lee_sub),
  type = "outcome",
  phenotype_col = "outcome",
  snp_col = "MarkerName",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "Pval",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "POS"
)

out_list <- list(
  "wightman_ad" = out_wightman,
  "wightmanukb_ad" = out_wightmanukb,
  "bellenguez_ad" = out_bellenguez,
  "lee_cp" = out_lee
)

out_combined <- Reduce(full_join, out_list)

##-----------------------------------------
## Create outcome list
##-----------------------------------------

outcome_list <- list()

for (pheno in names(instru_pqtl_list)) {
  outcome_list[[pheno]] <- out_combined %>%
    filter(SNP %in% instru_pqtl_list[[pheno]][["SNP"]])
}

##-----------------------------------------
## Save
##-----------------------------------------

save(
  instru_pqtl_list,
  outcome_list,
  file = "data/mr.RData"
)