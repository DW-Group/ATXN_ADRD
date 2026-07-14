library(dplyr)
library(data.table)
library(coloc)
library(readxl)
library(writexl)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx")

##-----------------------------------------
## coloc sig
##-----------------------------------------

### coloc.abf

sig_files <- list.files("data/colocalization/coloc_sig/", pattern = "coloc_sig_wightman_ad", full.names = FALSE)
sig_prot_rsid  <- gsub("coloc_sig_wightman_ad_|\\.tsv$", "", sig_files)

coloc_res_list <- list(); i <- 1

for (prot_rsid in sig_prot_rsid){
  
  message(prot_rsid)
  
  for (out in c("wightman_ad","wightmanukb_ad","bellenguez_ad","lee")){
    
    message(out)
    
    tmp_dat <- fread(paste0("data/colocalization/coloc_sig/coloc_sig_", out, "_", prot_rsid, ".tsv")) %>%
      filter(
        is.finite(beta_exp),
        is.finite(varbeta_exp),
        is.finite(beta_out),
        is.finite(varbeta_out),
        !is.na(snp),
        !is.na(position_exp),
        !is.na(position_out)
      )

    tmp_exp_list <- list()
    tmp_exp_list[["beta"]] <- tmp_dat$beta_exp
    tmp_exp_list[["varbeta"]] <- tmp_dat$varbeta_exp
    tmp_exp_list[["snp"]] <- tmp_dat$snp
    tmp_exp_list[["position"]] <- if ("position_exp" %in% names(tmp_dat)) tmp_dat$position_exp else tmp_dat$position
    tmp_exp_list[["type"]] <- "quant"
    tmp_exp_list[["sdY"]] <- 1
    if ("N_exp" %in% names(tmp_dat)) {
      tmp_exp_list[["N"]] <- median(tmp_dat$N_exp, na.rm = TRUE)
    }
    
    tmp_out_list <- list()
    tmp_out_list[["beta"]] <- tmp_dat$beta_out
    tmp_out_list[["varbeta"]] <- tmp_dat$varbeta_out
    tmp_out_list[["snp"]] <- tmp_dat$snp
    tmp_out_list[["position"]] <- if ("position_out" %in% names(tmp_dat)) tmp_dat$position_out else tmp_dat$position
    
    if (out == "lee"){
      tmp_out_list[["type"]] <- "quant"
      tmp_out_list[["sdY"]] <- 1
      if ("N_out" %in% names(tmp_dat)) {
        tmp_out_list[["N"]] <- median(tmp_dat$N_out, na.rm = TRUE)
      }
    } else{
      tmp_out_list[["type"]] <- "cc"
      if ("N_out" %in% names(tmp_dat)) {
        tmp_out_list[["N"]] <- median(tmp_dat$N_out, na.rm = TRUE)
      }
      if ("s_out" %in% names(tmp_dat)) {
        tmp_out_list[["s"]] <- median(tmp_dat$s_out, na.rm = TRUE)
      }
    }
    
    tmp_res <- coloc.abf(dataset1 = tmp_exp_list, dataset2 = tmp_out_list)
    coloc_res_list[[i]] <- c(prot_rsid, out, as.vector(tmp_res[["summary"]]))
    i <- i + 1
  }
}

coloc_res_dat <- as.data.frame(do.call("rbind", coloc_res_list))
colnames(coloc_res_dat) <- c("targetName","outcome","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")

coloc_res_dat <- coloc_res_dat %>%
  mutate(
    nsnps = as.numeric(nsnps),
    PP.H0.abf = as.numeric(PP.H0.abf),
    PP.H1.abf = as.numeric(PP.H1.abf),
    PP.H2.abf = as.numeric(PP.H2.abf),
    PP.H3.abf = as.numeric(PP.H3.abf),
    PP.H4.abf = as.numeric(PP.H4.abf),
    PP3_plus_PP4 = PP.H3.abf + PP.H4.abf,
    PP4_div_PP3_plus_PP4 = PP.H4.abf / PP3_plus_PP4,
    outcome_exposure = paste0(outcome, "_", targetName)
  )

write_xlsx(coloc_res_dat, "results/colocalization/coloc_sig_res.xlsx")

##-----------------------------------------
## coloc cis
##-----------------------------------------

### coloc.abf

rm(list = ls())

marker_dict_final <- read_excel("data/NULISA_pQTL_dictionary.xlsx")

coloc_res_list <- list(); i <- 1

for (prot in marker_dict_final %>% filter(chr != "X") %>% pull(targetName_unique)){
  
  message(prot)
  
  for (out in c("wightman_ad", "wightmanukb_ad", "bellenguez_ad", "lee")){
    
    message(out)
    
    tmp_dat <- fread(paste0("data/colocalization/coloc_cis/coloc_cis_", out, "_", prot, ".tsv")) %>%
      filter(
        is.finite(beta_exp),
        is.finite(varbeta_exp),
        is.finite(beta_out),
        is.finite(varbeta_out),
        !is.na(snp),
        !is.na(position_exp),
        !is.na(position_out)
      )
    
    if (nrow(tmp_dat) < 2) next
    
    tmp_exp_list <- list()
    tmp_exp_list[["beta"]] <- tmp_dat$beta_exp
    tmp_exp_list[["varbeta"]] <- tmp_dat$varbeta_exp
    tmp_exp_list[["snp"]] <- tmp_dat$snp
    tmp_exp_list[["position"]] <- if ("position_exp" %in% names(tmp_dat)) tmp_dat$position_exp else tmp_dat$position
    tmp_exp_list[["type"]] <- "quant"
    tmp_exp_list[["sdY"]] <- 1
    if ("N_exp" %in% names(tmp_dat)) {
      tmp_exp_list[["N"]] <- tmp_dat$N_exp[1]
    }
    
    tmp_out_list <- list()
    tmp_out_list[["beta"]] <- tmp_dat$beta_out
    tmp_out_list[["varbeta"]] <- tmp_dat$varbeta_out
    tmp_out_list[["snp"]] <- tmp_dat$snp
    tmp_out_list[["position"]] <- if ("position_out" %in% names(tmp_dat)) tmp_dat$position_out else tmp_dat$position
    
    if (out == "lee"){
      tmp_out_list[["type"]] <- "quant"
      tmp_out_list[["sdY"]] <- 1
      if ("N_out" %in% names(tmp_dat)) {
        tmp_out_list[["N"]] <- tmp_dat$N_out[1]
      }
    } else {
      tmp_out_list[["type"]] <- "cc"
      if ("N_out" %in% names(tmp_dat)) {
        tmp_out_list[["N"]] <- tmp_dat$N_out[1]
      }
      if ("s_out" %in% names(tmp_dat)) {
        tmp_out_list[["s"]] <- tmp_dat$s_out[1]
      }
    }
    
    tmp_res <- coloc.abf(dataset1 = tmp_exp_list, dataset2 = tmp_out_list)
    coloc_res_list[[i]] <- c(prot, out, as.vector(tmp_res[["summary"]]))
    i <- i + 1
  }
}

coloc_res_dat <- as.data.frame(do.call("rbind", coloc_res_list))
colnames(coloc_res_dat) <- c("targetName", "outcome", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

coloc_res_dat <- coloc_res_dat %>%
  mutate(
    nsnps = as.numeric(nsnps),
    PP.H0.abf = as.numeric(PP.H0.abf),
    PP.H1.abf = as.numeric(PP.H1.abf),
    PP.H2.abf = as.numeric(PP.H2.abf),
    PP.H3.abf = as.numeric(PP.H3.abf),
    PP.H4.abf = as.numeric(PP.H4.abf),
    PP3_plus_PP4 = PP.H3.abf + PP.H4.abf,
    PP4_div_PP3_plus_PP4 = PP.H4.abf / PP3_plus_PP4,
    outcome_exposure = paste0(outcome, "_", targetName)
  )

addmargins(table(coloc_res_dat$PP.H4.abf > 0.7, coloc_res_dat$PP4_div_PP3_plus_PP4 > 0.7))
table(coloc_res_dat$PP.H4.abf > 0.7)
table(coloc_res_dat$PP.H3.abf > 0.7)

write_xlsx(coloc_res_dat, "results/colocalization/coloc_cis_res.xlsx")
