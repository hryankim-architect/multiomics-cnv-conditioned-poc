"""mocnv — hypothesis-conditioned multi-omics with copy-number as a third modality.

A clean-room capability portrait: the DMOI pole-conditioned architecture (RNA +
DNA-methylation) extended with a third **copy-number (CNV)** branch and pole mask.
The deliverable is an honest **modality ablation** — does adding CNV over RNA+meth
help, per task axis, and does per-pole attribution key on the expected amplicons
(ERBB2 for HER2; MYC / CCND1 for proliferation)? A null or axis-specific result is
acceptable and on-brand (the dmoi v0.13 posture).

This package ships the harness and the engineering. Real cohorts (TCGA CNV via
GISTIC2; METABRIC SNP6 CNA) are downloaded by the user; the unit tests run on
small synthetic fixtures only.
"""
from __future__ import annotations

__version__ = "0.0.0"

# Modality order used everywhere (encoders, masks, fusion, ablation labels).
MODALITIES = ("rna", "meth", "cnv")
