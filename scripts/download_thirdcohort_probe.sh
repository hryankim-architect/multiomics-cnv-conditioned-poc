#!/usr/bin/env bash
# v0.6 (3) download-probe: find a third BRCA cohort (not TCGA, not METABRIC) with
# genome-wide gene-level CNA + a HER2/subtype label, for external validation of the
# strength->transfer law. Tries several cBioPortal datahub studies via the git-LFS media
# endpoint (same mechanism that worked for METABRIC), validates each, and reports which
# is usable. Saves the first fully-valid candidate to data/thirdcohort/<study>/.
#
# Usage: bash scripts/download_thirdcohort_probe.sh   (run on the Mac; paste the report back)
set -uo pipefail

BASE="https://media.githubusercontent.com/media/cBioPortal/datahub/master/public"
OUT_DIR="data/thirdcohort"
mkdir -p "${OUT_DIR}"

# Candidate studies with plausibly genome-wide CNA + clinical subtype/HER2 (not TCGA/METABRIC).
CANDIDATES=(
  brca_smc_2018
  brca_mbcproject_wagle_2017
  brca_igr_2015
  brca_broad
  brca_sanger
  brca_bccrc
)

probe_one() {
  local study="$1"
  local d="${OUT_DIR}/${study}"
  mkdir -p "${d}"
  echo "----- ${study} -----"

  # 1) gene-level CNA
  curl -fsSL --retry 2 -o "${d}/data_cna.txt" "${BASE}/${study}/data_cna.txt" 2>/dev/null || true
  if [[ -s "${d}/data_cna.txt" ]] && head -1 "${d}/data_cna.txt" | grep -q "^Hugo_Symbol"; then
    local genes samples
    genes=$(($(wc -l < "${d}/data_cna.txt") - 1))
    samples=$(($(head -1 "${d}/data_cna.txt" | tr '\t' '\n' | wc -l) - 2))
    echo "  CNA: OK  (~${genes} genes x ${samples} samples, $(du -h "${d}/data_cna.txt" | cut -f1))"
    if [[ "${genes}" -lt 5000 ]]; then echo "  CNA: WARNING - only ${genes} genes (panel, not genome-wide?)"; fi
  else
    echo "  CNA: missing / not a valid table (skip)"
    rm -f "${d}/data_cna.txt"; rmdir "${d}" 2>/dev/null; return
  fi

  # 2) clinical with a HER2 / subtype label
  for cf in data_clinical_sample.txt data_clinical_patient.txt; do
    curl -fsSL --retry 2 -o "${d}/${cf}" "${BASE}/${study}/${cf}" 2>/dev/null || true
    [[ -s "${d}/${cf}" ]] || { rm -f "${d}/${cf}"; continue; }
    local hits
    hits=$(grep -v '^#' "${d}/${cf}" | head -1 | tr '\t' '\n' \
           | grep -inE "her2_status|claudin_subtype|pam50|subtype|er_status" | tr '\n' ' ')
    [[ -n "${hits}" ]] && echo "  ${cf}: label cols -> ${hits}"
  done
  echo "  ERBB2 present in CNA: $(cut -f1 "${d}/data_cna.txt" | grep -qx ERBB2 && echo yes || echo NO)"
}

for s in "${CANDIDATES[@]}"; do probe_one "${s}"; done

echo ""
echo "Report the per-candidate lines above. A usable 3rd cohort needs: CNA OK (>=5000 genes,"
echo "ERBB2 present) AND a HER2/subtype label column. We build the loader+eval for the first"
echo "that qualifies; if none do, v0.6 (3) is recorded as 'no clean 3rd cohort readily available'."
