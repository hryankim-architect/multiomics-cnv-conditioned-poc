"""CNV pole masks — amplicon-locus gene sets for hypothesis-conditioning.

Copy number is amplicon-structured, so a CNV pole "mask" is a small set of genes
at the locus a hypothesis says matters. Unlike RNA pathway masks (dozens-to-hundreds
of genes), these are tight, locus-anchored sets:

  HER2 pole          -> ERBB2 amplicon at 17q12
  Proliferation pole -> MYC amplicon at 8q24  +  CCND1 amplicon at 11q13

Gene symbols are facts; the sets below are the canonical co-amplified genes in the
smallest regions of recurrent amplification at each locus. The keystone anchors
(ERBB2 / MYC / CCND1) are what the v0.2 attribution check verifies the CNV branch
keys on — co-amplified neighbours are included so the mask is not a single feature.
"""
from __future__ import annotations

import numpy as np

# 17q12 ERBB2 amplicon (HER2-driven subtypes).
POLE_HER2_CNV: tuple[str, ...] = (
    "ERBB2", "GRB7", "STARD3", "MIEN1", "PGAP3", "TCAP", "PNMT", "MED1",
)

# 8q24 MYC amplicon + 11q13 CCND1 amplicon (proliferation-driven).
POLE_PROLIFERATION_CNV: tuple[str, ...] = (
    "MYC", "PVT1", "CASC8", "POU5F1B",                      # 8q24
    "CCND1", "ORAOV1", "FGF3", "FGF4", "FGF19", "ANO1", "FADD", "CTTN",  # 11q13
)

CNV_POLES: dict[str, tuple[str, ...]] = {
    "HER2": POLE_HER2_CNV,
    "Proliferation": POLE_PROLIFERATION_CNV,
}

# The single anchor gene per locus — the amplification driver the attribution
# check (v0.2) expects the CNV branch to weight most.
ANCHOR_GENES: dict[str, str] = {
    "HER2": "ERBB2",
    "Proliferation_8q24": "MYC",
    "Proliferation_11q13": "CCND1",
}


def cnv_pole_mask(pole: str, gene_universe: list[str]) -> np.ndarray:
    """Boolean mask (len == len(gene_universe)) selecting the pole's amplicon genes.

    Genes in the pole set that are absent from ``gene_universe`` are simply not
    selected (the data layer decides coverage); see ``mask_coverage`` to report it.
    """
    if pole not in CNV_POLES:
        raise KeyError(f"unknown CNV pole {pole!r}; have {sorted(CNV_POLES)}")
    wanted = set(CNV_POLES[pole])
    return np.array([g in wanted for g in gene_universe], dtype=bool)


def mask_coverage(pole: str, gene_universe: list[str]) -> dict[str, object]:
    """How many of the pole's genes are present in the universe (honest coverage)."""
    present = [g for g in CNV_POLES[pole] if g in set(gene_universe)]
    missing = [g for g in CNV_POLES[pole] if g not in set(gene_universe)]
    anchors_present = [a for k, a in ANCHOR_GENES.items() if k.startswith(pole) and a in present]
    return {
        "pole": pole,
        "n_in_set": len(CNV_POLES[pole]),
        "n_present": len(present),
        "present": present,
        "missing": missing,
        "anchors_present": anchors_present,
    }
