module load regenie/3.0.3

# Phenotype info
PHENO_FILE="data/input_regenie_NULISA_dp_pheno.txt"
PHENO_NAME=$(head -1 "$PHENO_FILE" | cut -f3- | tr '\t' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")
COVAR="data/input_regenie_NULISA_dp_covar.txt"

# Output directory
OUTDIR="results/regenie_NULISA_dp"
STEP1_OUT="${OUTDIR}/step1_dp_${PHENO_NAME}"
STEP2_OUT="${OUTDIR}/step2_dp_${PHENO_NAME}"

echo "[`date`] Starting GWAS for phenotype: $PHENO_NAME"

# Step 1: Fit null model
regenie \
  --step 1 \
  --pgen "data/genotypes/for_regenie_step_1/DP_allchr_genotyped" \
  --phenoFile "$PHENO_FILE" \
  --covarFile "$COVAR" \
  --phenoCol "$PHENO_NAME" \
  --bsize 1000 \
  --threads 8 \
  --extract "data/genotypes/for_regenie_step_1/DP_1000G_pruned_snps_for_regenie_NULISA.prune.in" \
  --out "$STEP1_OUT"

# Step 2: Association test
PGEN="data/genotypes/DIRECT-PLUS_GSA_all_chr_1000G_eur_nonrel"
STEP2_CHR_OUT="${STEP2_OUT}_allchr"

echo "[`date`] Step 2: GWAS for $PHENO_NAME"

regenie \
  --step 2 \
  --pgen "$PGEN" \
  --phenoFile "$PHENO_FILE" \
  --covarFile "$COVAR" \
  --phenoCol "$PHENO_NAME" \
  --pred "${STEP1_OUT}_pred.list" \
  --bsize 400 \
  --threads 8 \
  --minMAC 5 \
  --htp "dp" \
  --out "$STEP2_CHR_OUT"

echo "[`date`] Completed GWAS for $PHENO_NAME"