module load regenie/3.0.3

# Phenotype info
PHENO_FILE="data/input_regenie_NULISA_nhs_final_pheno.txt"
PHENO_NAME=$(head -1 "$PHENO_FILE" | cut -f3- | tr '\t' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")
COVAR="data/input_regenie_NULISA_nhs_final_covar.txt"

# Output directory
OUTDIR="results/regenie_NULISA_nhs"

echo "[`date`] Starting GWAS for phenotype: $PHENO_NAME"

for chr in {1..23}; do

  STEP1_OUT="${OUTDIR}/step1_nhs_${PHENO_NAME}_chr${chr}"
  STEP2_OUT="${OUTDIR}/step2_nhs_${PHENO_NAME}_chr${chr}"

  # Step 1: Fit null model
  regenie \
    --step 1 \
    --pgen "data/genotypes/nhs/for_regenie_step1/NULISA_reimputation_chr${chr}_1000G_for_regenie_step1" \
    --phenoFile "$PHENO_FILE" \
    --covarFile "$COVAR" \
    --phenoCol "$PHENO_NAME" \
    --bsize 1000 \
    --threads 8 \
    --extract "data/genotypes/nhs/for_regenie_step1/NULISA_reimputation_chr${chr}_1000G_for_regenie_step1_prune.prune.in" \
    --out "$STEP1_OUT"

  # Step 2: Association test
  PGEN="data/genotypes/nhs/plink2/NULISA_reimputation_chr${chr}_1000G"

  echo "[`date`] Step 2: GWAS for $PHENO_NAME chr${chr}"

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
    --htp "nhs" \
    --out "$STEP2_OUT"

done

echo "[`date`] Completed GWAS for $PHENO_NAME"