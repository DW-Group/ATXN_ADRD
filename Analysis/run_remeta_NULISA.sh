module purge
module load gcc/8.2.0
module load cmake/3.23.0

export LD_LIBRARY_PATH="path/to/boost/lib:$LD_LIBRARY_PATH"

REMETA="path/to/remeta"
WORKDIR="results"
OUTDIR="results/remeta_NULISA"
PHENO_FILE="data/input_regenie_NULISA_nhs_final_pheno.txt"

mkdir -p "$OUTDIR"

PHENO_NAME=$(head -1 "$PHENO_FILE" | cut -f3- | tr '\t' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "PHENO_NAME=${PHENO_NAME}"

NHS_FILE="${WORKDIR}/regenie_NULISA_nhs_qc/step2_nhs_${PHENO_NAME}_allchr_${PHENO_NAME}.qc.regenie"
DP_FILE="${WORKDIR}/regenie_NULISA_dp_qc/step2_dp_${PHENO_NAME}_allchr_${PHENO_NAME}.qc.regenie"

HTP_FILES=()
COHORTS=()

if [ -f "${NHS_FILE}" ]; then
  HTP_FILES+=("${NHS_FILE}")
  COHORTS+=("nhs")
fi

if [ -f "${DP_FILE}" ]; then
  HTP_FILES+=("${DP_FILE}")
  COHORTS+=("dp")
fi

"${REMETA}" esma \
  --htp "${HTP_FILES[@]}" \
  --cohorts "${COHORTS[@]}" \
  --trait-name "${PHENO_NAME}" \
  --trait-type QT \
  --out "${OUTDIR}/remeta_${PHENO_NAME}"