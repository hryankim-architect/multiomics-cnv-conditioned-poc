#!/usr/bin/env bash
# Download TCGA-BRCA GISTIC2 gene-level copy number (UCSC Xena).
# Output: data/tcga_brca/Gistic2_CopyNumber_all_data_by_genes.gz
set -euo pipefail

URL="https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"
OUT_DIR="data/tcga_brca"
OUT="${OUT_DIR}/Gistic2_CopyNumber_all_data_by_genes.gz"

mkdir -p "${OUT_DIR}"

if [[ -s "${OUT}" ]]; then
  echo "Already present: ${OUT} ($(du -h "${OUT}" | cut -f1))"
else
  echo "Downloading TCGA-BRCA GISTIC2 gene-level CNV ..."
  curl -fL --retry 3 -o "${OUT}" "${URL}"
fi

if [[ ! -s "${OUT}" ]]; then
  echo "ERROR: download produced an empty file: ${OUT}" >&2
  exit 1
fi

echo "Saved: ${OUT} ($(du -h "${OUT}" | cut -f1))"
echo "SHA-256 (paste into data/manifest.yaml -> tcga_brca_cnv.sha256):"
shasum -a 256 "${OUT}" 2>/dev/null || sha256sum "${OUT}"
