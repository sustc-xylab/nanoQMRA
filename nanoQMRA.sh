#!/bin/bash
set -euo pipefail

show_help() {
cat << EOF
Usage: ${0##*/} -f reads.fa -t threads [-p Plascad_env]

version 1.1

Arguments:
  -h                display this help
  -f reads.fa       input nanopore reads fasta file
  -t threads        number of threads
  -p env_name       conda environment name for Plascad (default: Plascad)

Description:
  nanoQMRA wraps:
    1) ARGpore2
    2) Plascad
    3) custom scoring scripts

Required in current ARGpore2 directory:
  - argpore.sh
  - scripts/build_and_score_mother_table.py
  - scripts/make_risk_score_list.py
  - rules/
  - input/ARG.csv

Example:
  bash nanoQMRA.sh -f test.fa -t 60
EOF
}

SCRIPT=$(realpath "$0")
DIR=$(dirname "$SCRIPT")
NOWT=$(date +%Y-%m-%d.%H:%M:%S)

INPUT_FA=""
THREADS=""
PLASCAD_ENV="Plascad"

while getopts "f:t:p:h" opt; do
  case "$opt" in
    h)
      show_help
      exit 0
      ;;
    f)
      INPUT_FA="$OPTARG"
      ;;
    t)
      THREADS="$OPTARG"
      ;;
    p)
      PLASCAD_ENV="$OPTARG"
      ;;
    *)
      echo "[ERROR] Invalid option"
      show_help >&2
      exit 1
      ;;
  esac
done

if [ -z "${INPUT_FA}" ]; then
  echo "[ERROR] -f must be specified"
  exit 1
fi

if [ -z "${THREADS}" ]; then
  echo "[ERROR] -t must be specified"
  exit 1
fi

if [ ! -f "${INPUT_FA}" ]; then
  echo "[ERROR] input fasta not found: ${INPUT_FA}"
  exit 1
fi

if [ ! -f "${DIR}/argpore.sh" ]; then
  echo "[ERROR] argpore.sh not found in ${DIR}"
  exit 1
fi

if [ ! -f "${DIR}/scripts/build_and_score_mother_table.py" ]; then
  echo "[ERROR] scripts/build_and_score_mother_table.py not found"
  exit 1
fi

if [ ! -f "${DIR}/scripts/make_risk_score_list.py" ]; then
  echo "[ERROR] scripts/make_risk_score_list.py not found"
  exit 1
fi

if [ ! -d "${DIR}/rules" ]; then
  echo "[ERROR] rules directory not found"
  exit 1
fi

if [ ! -d "${DIR}/input" ]; then
  echo "[ERROR] input directory not found"
  exit 1
fi

if [ ! -f "${DIR}/input/ARG.csv" ]; then
  echo "[ERROR] input/ARG.csv not found"
  exit 1
fi

echo "
--------------------------------------------------
Start nanoQMRA @ $(date +"%Y-%m-%d %T")
Input fasta: ${INPUT_FA}
Threads: ${THREADS}
Plascad env: ${PLASCAD_ENV}
Working dir: ${DIR}
--------------------------------------------------
"

INPUT_NAME=$(basename "${INPUT_FA}")
INPUT_PREFIX="${INPUT_NAME%.*}"

mkdir -p "${DIR}/output"

############################################################
# Step 1. Run ARGpore2
############################################################
echo "[1/5] Running ARGpore2 ..."
bash "${DIR}/argpore.sh" -f "${INPUT_FA}" -t "${THREADS}" > "${DIR}/output/${INPUT_PREFIX}_ARGpore.log" 2>&1

# ARGpore2 output directory may contain timestamp and may be named like:
#   example.fa.uniq_ARGpore2_2026-04-22.09:44:02
#   example.fa_ARGpore2_2026-04-22.09:44:02
#   example.fa.uniq_ARGpore_*
# so search multiple possible patterns
ARGPORE_OUT=$(find "${DIR}" -maxdepth 1 -type d \( \
    -name "${INPUT_NAME}.uniq_ARGpore2_*" -o \
    -name "${INPUT_NAME}_ARGpore2_*" -o \
    -name "${INPUT_NAME}.uniq_ARGpore_*" -o \
    -name "${INPUT_NAME}_ARGpore_*" \
  \) | sort | tail -n 1)

if [ -z "${ARGPORE_OUT}" ]; then
  echo "[ERROR] ARGpore2 output directory not found for ${INPUT_NAME}"
  echo "[INFO] Tried patterns:"
  echo "       ${INPUT_NAME}.uniq_ARGpore2_*"
  echo "       ${INPUT_NAME}_ARGpore2_*"
  echo "       ${INPUT_NAME}.uniq_ARGpore_*"
  echo "       ${INPUT_NAME}_ARGpore_*"
  exit 1
fi

echo "[INFO] ARGpore2 output directory: ${ARGPORE_OUT}"

# Find taxa result file, e.g.:
#   example.fa.uniq_arg.w.taxa.tab
ARG_TAXA_FILE=$(find "${ARGPORE_OUT}" -type f \( \
    -name "${INPUT_NAME}.uniq_arg.w.taxa.tab" -o \
    -name "${INPUT_NAME}*.arg.w.taxa.tab" -o \
    -name "${INPUT_NAME}*.w.taxa.tab" -o \
    -name "*.arg.w.taxa.tab" -o \
    -name "*.w.taxa.tab" \
  \) | sort | head -n 1)

if [ -z "${ARG_TAXA_FILE}" ]; then
  echo "[ERROR] ARG taxa result file not found in ${ARGPORE_OUT}"
  exit 1
fi

echo "[INFO] ARG taxa file: ${ARG_TAXA_FILE}"
cp "${ARG_TAXA_FILE}" "${DIR}/input/"
############################################################
# Step 2. Run Plascad
############################################################
echo "[2/5] Running Plascad ..."

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda not found in PATH"
  exit 1
fi

CONDA_BASE=$(conda info --base)
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${PLASCAD_ENV}"

Plascad -i "${INPUT_FA}" > "${DIR}/output/${INPUT_PREFIX}_Plascad.log" 2>&1 || {
  echo "[ERROR] Plascad failed"
  conda deactivate || true
  exit 1
}

conda deactivate || true

# Find Plascad output summary, expected like test_plasmids_classification_sum.txt
PLASCAD_SUM=$(find "${DIR}" -maxdepth 1 -type f -name "${INPUT_PREFIX}*_plasmids_classification_sum.txt" | sort | tail -n 1)

if [ -z "${PLASCAD_SUM}" ]; then
  PLASCAD_SUM=$(find "${DIR}" -maxdepth 1 -type f -name "*_plasmids_classification_sum.txt" | sort | tail -n 1)
fi

if [ -z "${PLASCAD_SUM}" ]; then
  echo "[ERROR] Plascad summary file not found"
  exit 1
fi

echo "[INFO] Plascad summary file: ${PLASCAD_SUM}"
cp "${PLASCAD_SUM}" "${DIR}/input/"

############################################################
# Step 3. Build mother table
############################################################
echo "[3/5] Building scored mother table ..."

INPUT_ARG_TAXA=$(find "${DIR}/input" -maxdepth 1 -type f -name "${INPUT_NAME}*.arg.w.taxa.tab" | head -n 1)
if [ -z "${INPUT_ARG_TAXA}" ]; then
  INPUT_ARG_TAXA=$(find "${DIR}/input" -maxdepth 1 -type f -name "${INPUT_NAME}*.w.taxa.tab" | head -n 1)
fi

INPUT_PLASCAD_SUM=$(find "${DIR}/input" -maxdepth 1 -type f -name "${INPUT_PREFIX}*_plasmids_classification_sum.txt" | head -n 1)

if [ -z "${INPUT_ARG_TAXA}" ]; then
  echo "[ERROR] input ARG taxa file not found in input/"
  exit 1
fi

if [ -z "${INPUT_PLASCAD_SUM}" ]; then
  echo "[ERROR] input Plascad summary file not found in input/"
  exit 1
fi

MOTHER_OUT="${DIR}/output/${INPUT_PREFIX}_mother_table.scored.tsv"

python "${DIR}/scripts/build_and_score_mother_table.py" \
  --arg_taxa "${INPUT_ARG_TAXA}" \
  --plasmid_sum "${INPUT_PLASCAD_SUM}" \
  --rules "${DIR}/rules" \
  --arg_csv "${DIR}/input/ARG.csv" \
  -o "${MOTHER_OUT}"

############################################################
# Step 4. Build risk score list
############################################################
echo "[4/5] Building risk_score_list.tsv ..."

RISK_LIST_OUT="${DIR}/output/${INPUT_PREFIX}_risk_score_list.tsv"

python "${DIR}/scripts/make_risk_score_list.py" \
  -i "${MOTHER_OUT}" \
  -o "${RISK_LIST_OUT}"

############################################################
# Step 5. Calculate final sample-level risk_score
############################################################
echo "[5/6] Calculating final risk_score ..."

TOTAL_READS=$(grep -c "^>" "${INPUT_FA}" || true)

if [ -z "${TOTAL_READS}" ] || [ "${TOTAL_READS}" -eq 0 ]; then
  echo "[ERROR] total reads count is zero; cannot calculate final risk_score"
  exit 1
fi

SUM_SCORE=$(awk -F'\t' 'NR>1 && $4 != "" {sum += $4} END {print sum+0}' "${RISK_LIST_OUT}")
FINAL_SCORE=$(awk -v s="${SUM_SCORE}" -v n="${TOTAL_READS}" 'BEGIN {printf "%.6f", s/n}')

{
  echo ""
  echo -e "risk_score:\t${FINAL_SCORE}"
} >> "${RISK_LIST_OUT}"

echo ""
echo "[OK] total reads       : ${TOTAL_READS}"
echo "[OK] summed risk score : ${SUM_SCORE}"
echo "[OK] final risk_score  : ${FINAL_SCORE}"
echo "[OK] mother table      : ${MOTHER_OUT}"
echo "[OK] risk score list   : ${RISK_LIST_OUT}"

echo ""
echo "--------------------------------------------------"

############################################################
# Step 6. Organize intermediate files
############################################################
echo "[6/6] Organizing intermediate files ..."

INTER_DIR="${DIR}/intermediate/${INPUT_PREFIX}"
mkdir -p "${INTER_DIR}"

# Move ARGpore2 output directory
if [ -n "${ARGPORE_OUT:-}" ] && [ -d "${ARGPORE_OUT}" ]; then
  mv "${ARGPORE_OUT}" "${INTER_DIR}/"
fi

# Move Plascad intermediate directory if present
PLASCAD_MAP_DIR="${DIR}/${INPUT_PREFIX}_Conjugative_plasmids_map"
if [ -d "${PLASCAD_MAP_DIR}" ]; then
  mv "${PLASCAD_MAP_DIR}" "${INTER_DIR}/"
fi

# Move Plascad summary-related files if present
for f in \
  "${DIR}/${INPUT_PREFIX}_Conj_plasmids_loc_sum.txt" \
  "${DIR}/${INPUT_PREFIX}_mob_unconj_plasmids_loc_sum.txt" \
  "${DIR}/${INPUT_PREFIX}_plasmids_classification_sum.txt"
do
  if [ -f "$f" ]; then
    mv "$f" "${INTER_DIR}/"
  fi
done

echo "[OK] intermediate files moved to: ${INTER_DIR}"

echo "--------------------------------------------------"
echo "Done nanoQMRA @ $(date +"%Y-%m-%d %T")"
echo "--------------------------------------------------"