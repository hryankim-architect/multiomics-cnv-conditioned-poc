#!/usr/bin/env bash
# Download METABRIC SNP6 gene-level CNA from cBioPortal (brca_metabric bundle)
# and extract just data_CNA.txt.
# Output: data/metabric/data_CNA.txt
set -euo pipefail

URL="https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz"
OUT_DIR="data/metabric"
TARBALL="${OUT_DIR}/brca_metabric.tar.gz"
OUT="${OUT_DIR}/data_CNA.txt"

mkdir -p "${OUT_DIR}"

if [[ -s "${OUT}" ]]; then
  echo "Already present: ${OUT} ($(du -h "${OUT}" | cut -f1))"
else
  echo "Downloading METABRIC bundle (only data_CNA.txt is kept) ..."
  curl -fL --retry 3 -o "${TARBALL}" "${URL}"
  echo "Extracting data_CNA.txt ..."
  # The CNA table lives at brca_metabric/data_CNA.txt inside the archive.
  tar -xzf "${TARBALL}" -C "${OUT_DIR}" --strip-components=1 "brca_metabric/data_CNA.txt"
  rm -f "${TARBALL}"
fi

if [[ ! -s "${OUT}" ]]; then
  echo "ERROR: data_CNA.txt missing or empty after extraction: ${OUT}" >&2
  exit 1
fi

echo "Saved: ${OUT} ($(du -h "${OUT}" | cut -f1))"
echo "SHA-256 (paste into data/manifest.yaml -> metabric_cna.sha256):"
shasum -a 256 "${OUT}" 2>/dev/null || sha256sum "${OUT}"
