#!/usr/bin/env python3
"""v0.6 (2) — is the strength->transfer law robust to how 'strength' is defined?

v0.5 (1) used within-cohort single-gene AUROC as 'strength' (pooled Spearman vs
cross-cohort transfer = +0.836). v0.6 (2) re-tests with a different, model-free
metric: **raw GISTIC2 amplification amplitude** (mean positive - mean negative,
un-z-scored). If the law holds with amplitude too, it is not an artifact of the AUROC
definition. sklearn/numpy (no torch). Reproduce:
  python scripts/build_metabric_cohort_v2.py && python scripts/amplicon_strength_v0.6.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
DMOI = Path(os.environ.get("MOCNV_DMOI", str(REPO.parent / "dmoi-brca-poc")))

from mocnv import audit  # noqa: E402
from mocnv import cnv as cnvmod  # noqa: E402
from mocnv.priors import CNV_POLES  # noqa: E402
from mocnv.strength import (  # noqa: E402
    per_gene_amplitude,
    per_gene_strength_transfer,
    rank_spearman,
)

GISTIC2 = REPO / "data" / "tcga_brca" / "Gistic2_CopyNumber_all_data_by_genes.gz"
METABRIC_CNA = REPO / "data" / "metabric" / "data_CNA.txt"
AUDIT = REPO / "audit"
AMPLICON = sorted({g for genes in CNV_POLES.values() for g in genes})
AXES = [("HER2", "cohort_v4.tsv", "HER2"), ("LumB", "cohort_v2.tsv", "LumB")]


def _align_raw(matrix, want, *, tcga):
    """Amplicon-aligned RAW CNV (no z-score) + keep-index."""
    samp = cnvmod.barcode_to_sample(matrix.sample_ids) if tcga else matrix.sample_ids
    idx = {s: i for i, s in enumerate(samp)}
    rows, keep = [], []
    for k, sid in enumerate(want):
        s = cnvmod.barcode_to_sample([sid])[0] if tcga else sid
        if s in idx:
            rows.append(idx[s])
            keep.append(k)
    return cnvmod.align_to_genes(matrix.values[rows], matrix.gene_names, AMPLICON, fill_value=0.0), keep


def run_axis(cohort_file: str, positive: str) -> dict | None:
    tco_path, mco_path = DMOI / "data/tcga_brca" / cohort_file, DMOI / "data/metabric" / cohort_file
    if not tco_path.exists() or not mco_path.exists():
        sys.stderr.write(f"skip {positive}: missing cohort\n")
        return None
    tco = pd.read_csv(tco_path, sep="\t")
    tco = tco[tco["has_rna"] & tco["has_meth"]]
    raw_t, keep = _align_raw(cnvmod.load_gistic2(GISTIC2), tco["sample_id"].astype(str).tolist(), tcga=True)
    y_t = (tco["group"] == positive).astype(int).to_numpy()[keep]

    mco = pd.read_csv(mco_path, sep="\t")
    mco = mco[mco["has_rna"]]
    ymap = dict(zip(mco["sample_id"].astype(str), (mco["group"] == positive).astype(int), strict=False))
    cna = cnvmod.load_cbioportal_cna(METABRIC_CNA, sample_ids=set(ymap))
    raw_m, keepm = _align_raw(cna, list(cna.sample_ids), tcga=False)
    y_m = np.array([ymap[cna.sample_ids[i]] for i in keepm])

    z_t = cnvmod.harmonize_gene_level(raw_t, method="zscore")
    z_m = cnvmod.harmonize_gene_level(raw_m, method="zscore")
    within, transfer = per_gene_strength_transfer(z_t, y_t, z_m, y_m)
    amplitude = per_gene_amplitude(raw_t, y_t)
    return {"axis": positive, "within": within, "transfer": transfer, "amplitude": amplitude,
            "rho_within": rank_spearman(within, transfer),
            "rho_amp": rank_spearman(amplitude, transfer)}


def main() -> int:
    for p in (GISTIC2, METABRIC_CNA):
        if not p.exists():
            sys.stderr.write(f"missing input: {p}\n")
            return 1
    AUDIT.mkdir(exist_ok=True)
    results = [r for _, cf, pos in AXES if (r := run_axis(cf, pos)) is not None]
    if not results:
        sys.stderr.write("no axes ran\n")
        return 1

    aw, at, aa = [], [], []
    for r in results:
        aw.extend(r["within"])
        at.extend(r["transfer"])
        aa.extend(r["amplitude"])
        print(f"  {r['axis']:5} rho(within-AUROC, transfer)={r['rho_within']:+.3f} "
              f"| rho(raw amplitude, transfer)={r['rho_amp']:+.3f}")
    pooled_w = rank_spearman(np.array(aw), np.array(at))
    pooled_a = rank_spearman(np.array(aa), np.array(at))
    print(f"\n  POOLED rho(within-AUROC, transfer)={pooled_w:+.3f} "
          f"| rho(raw amplitude, transfer)={pooled_a:+.3f}")

    audit.emit("amplicon_strength_v0.6", "robustness",
               {"pooled_within": round(pooled_w, 3), "pooled_amplitude": round(pooled_a, 3)})
    md = AUDIT / "amplicon_strength_v0.6.md"
    body = "\n".join(
        f"| {r['axis']} | {r['rho_within']:+.3f} | {r['rho_amp']:+.3f} |" for r in results
    )
    md.write_text(
        "# Strength->transfer law robustness (raw amplitude) — v0.6 (2)\n\n"
        "v0.5 (1) measured strength as within-cohort single-gene AUROC. Here strength is the "
        "**raw GISTIC2 amplification amplitude** (mean positive - mean negative, un-z-scored), "
        "a model-free metric. If the Spearman correlation with cross-cohort transfer stays "
        "strong, the law is not an artifact of the AUROC definition.\n\n"
        f"**Pooled (40 gene-points): rho(within-AUROC, transfer) = {pooled_w:+.3f}; "
        f"rho(raw amplitude, transfer) = {pooled_a:+.3f}.**\n\n"
        "| axis | rho(within-AUROC, transfer) | rho(raw amplitude, transfer) |\n"
        "|---|---|---|\n" + body + "\n\n"
        "Honest reading: a strong amplitude correlation that agrees with the AUROC one confirms "
        "the strength->transfer law is metric-robust — genes with larger copy-number amplitude "
        "in the positive class transfer better across platforms, however 'strength' is measured.\n",
    )
    print(f"\nwrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
