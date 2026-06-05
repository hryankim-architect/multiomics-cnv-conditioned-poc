# Data — multiomics-cnv-conditioned-poc

Public, aggregate molecular tables only — no patient-level identifiers. The unit
tests use **small synthetic fixtures** (`src/mocnv/synth.py`, from v0.1); the real
cohorts below are downloaded by the user and are gitignored.

## What is new here vs reused

Only the **CNV** modality is new. RNA and methylation are **reused from
`dmoi-brca-poc`** so sample IDs stay aligned across modalities — point this repo at
the dmoi data directory (or copy the relevant tables) rather than re-downloading.

## CNV downloads (run these)

| Cohort | File | Source | Script |
|---|---|---|---|
| TCGA-BRCA | GISTIC2 gene-level copy number | UCSC Xena | `scripts/download_tcga_cnv.sh` |
| METABRIC | SNP6 gene-level CNA (`data_CNA.txt`) | cBioPortal `brca_metabric` | `scripts/download_metabric_cna.sh` |

Each script downloads to `data/<cohort>/`, verifies the file is non-empty, and
prints a SHA-256 to paste back into `data/manifest.yaml` (the `sha256: TBD` fields).

## Honest note (cross-platform CNV)

TCGA GISTIC2 and METABRIC SNP6 CNA are **different platforms and pipelines**. They
are ingested at gene level and aligned to the shared RNA gene universe, but
cross-platform copy-number harmonization is not as clean as RNA quantile
normalization — the v0.3 cross-cohort check is explicitly weaker and caveated. See
`docs/architecture.md` and `docs/what-is-out-of-scope.md`.
