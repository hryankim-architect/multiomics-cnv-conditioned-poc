#!/usr/bin/env bash
# Download METABRIC SNP6 gene-level CNA (data_CNA.txt) from the cBioPortal datahub.
#
# The S3 tarball (cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz) returned
# HTTP 403 on 2026-06-04, so this fetches the single CNA file via the datahub git-LFS
# media endpoint. If that fails, see the documented fallbacks below.
# Output: data/metabric/data_CNA.txt
set -euo pipefail

URL="https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/data_CNA.txt"
OUT_DIR="data/metabric"
OUT="${OUT_DIR}/data_CNA.txt"
mkdir -p "${OUT_DIR}"

if [[ -s "${OUT}" ]] && head -1 "${OUT}" | grep -q "^Hugo_Symbol"; then
  echo "Already present: ${OUT} ($(du -h "${OUT}" | cut -f1))"
else
  echo "Downloading METABRIC data_CNA.txt (cBioPortal datahub, git-LFS media)..."
  curl -fL --retry 3 -o "${OUT}" "${URL}" || true
fi

# Validate: a real CNA table starts with the Hugo_Symbol header, not an LFS
# pointer / HTML error page.
if [[ ! -s "${OUT}" ]] || ! head -1 "${OUT}" | grep -q "^Hugo_Symbol"; then
  echo "ERROR: ${OUT} is not a valid CNA table." >&2
  echo "  got: $(head -c 80 "${OUT}" 2>/dev/null)" >&2
  echo "Fallbacks (pick one, then re-run to validate + hash):" >&2
  echo "  1) cBioPortal UI: Studies -> 'Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)'" >&2
  echo "     -> Download -> extract brca_metabric/data_CNA.txt into ${OUT_DIR}/" >&2
  echo "  2) git clone the LFS file: git lfs pull from cBioPortal/datahub public/brca_metabric/" >&2
  echo "  3) retry the S3 tarball (was 403): https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz" >&2
  exit 1
fi

echo "Saved: ${OUT} ($(du -h "${OUT}" | cut -f1))"
echo "SHA-256 (paste into data/manifest.yaml -> metabric_cna.sha256):"
shasum -a 256 "${OUT}" 2>/dev/null || sha256sum "${OUT}"
