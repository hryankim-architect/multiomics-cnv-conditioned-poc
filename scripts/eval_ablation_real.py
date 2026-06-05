#!/usr/bin/env python3
"""Real-cohort modality ablation: does CNV add signal over RNA+meth, per axis?

Wires the TCGA RNA + methylation feature pipeline from `dmoi-brca-poc` with the
GISTIC2 CNV modality (restricted to the amplicon pole masks), then runs the
{rna,meth} vs {rna,meth,cnv} ablation 5-fold per task axis, with per-pole
Integrated Gradients on the CNV branch.

The honest hypothesis: CNV helps the amplicon-driven axis (HER2, via the ERBB2
amplicon) and adds little to a non-amplicon axis (LumA-vs-LumB). A null delta is
a valid result.

REQUIRES `dmoi-brca-poc` checked out alongside (for the RNA/meth loaders + data):
  - RNA   : <dmoi>/data/tcga_brca/HiSeqV2.gz
  - meth  : <dmoi>/data/tcga_brca/HumanMethylation450.gz
  - cohort: <dmoi>/data/tcga_brca/cohort_v4.tsv (HER2 axis), cohort_v2.tsv (LumA/LumB)
CNV (this repo): data/tcga_brca/Gistic2_CopyNumber_all_data_by_genes.gz
Override the dmoi location with the MOCNV_DMOI env var if it lives elsewhere.

Reproduce (M-series Mac, torch; minutes after the one-time methylation load):
  python scripts/eval_ablation_real.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from sklearn.model_selection import StratifiedKFold

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
DMOI = Path(os.environ.get("MOCNV_DMOI", str(REPO.parent / "dmoi-brca-poc")))
sys.path.insert(0, str(DMOI / "src"))

import torch  # noqa: E402
from dmoi_brca.features import load_features  # noqa: E402

from mocnv import audit, eval_ablation  # noqa: E402
from mocnv import cnv as cnvmod  # noqa: E402
from mocnv.attribution import integrated_gradients, rank_genes_by_attribution  # noqa: E402
from mocnv.priors import CNV_POLES, POLE_HER2_CNV  # noqa: E402

TCGA = DMOI / "data" / "tcga_brca"
CNV_GZ = REPO / "data" / "tcga_brca" / "Gistic2_CopyNumber_all_data_by_genes.gz"
AUDIT = REPO / "audit"

# (axis label, cohort file under <dmoi>/data/tcga_brca, positive-class group)
AXES = [
    ("HER2-vs-Luminal", "cohort_v4.tsv", "HER2"),
    ("LumA-vs-LumB", "cohort_v2.tsv", "LumB"),
]
AMPLICON_GENES = sorted({g for genes in CNV_POLES.values() for g in genes})
LATENT_DIM = 128
N_EPOCHS = 120
SEED = 42


def _cnv_aligned(sample_ids: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Align the amplicon-masked, z-scored CNV matrix to the cohort's sample IDs.

    Returns (cnv [n_keep x n_amplicon], keep_idx into sample_ids, amplicon gene names).
    """
    m = cnvmod.load_gistic2(CNV_GZ)
    cnv_samp = cnvmod.barcode_to_sample(m.sample_ids)
    idx = {s: i for i, s in enumerate(cnv_samp)}
    rows, keep = [], []
    for k, sid in enumerate(sample_ids):
        s = cnvmod.barcode_to_sample([sid])[0]
        if s in idx:
            rows.append(idx[s])
            keep.append(k)
    sub = m.values[rows]
    aligned = cnvmod.align_to_genes(sub, m.gene_names, AMPLICON_GENES, fill_value=0.0)
    aligned = cnvmod.harmonize_gene_level(aligned, method="zscore")
    return aligned, np.array(keep, dtype=int), AMPLICON_GENES


def _five_fold(arrays: dict[str, np.ndarray], y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    base, full = [], []
    for tr, te in skf.split(arrays["rna"], y):
        base.append(eval_ablation.train_and_auroc(
            arrays, y, eval_ablation.BASELINE_SET, train_idx=tr, val_idx=te,
            latent_dim=LATENT_DIM, n_epochs=N_EPOCHS, seed=SEED))
        full.append(eval_ablation.train_and_auroc(
            arrays, y, eval_ablation.FULL_SET, train_idx=tr, val_idx=te,
            latent_dim=LATENT_DIM, n_epochs=N_EPOCHS, seed=SEED))
    return np.array(base), np.array(full)


def main() -> int:
    AUDIT.mkdir(exist_ok=True)
    rows = []
    for axis, cohort_file, pos in AXES:
        cohort = TCGA / cohort_file
        if not cohort.exists():
            print(f"skip {axis}: {cohort} not found")
            continue
        print(f"\n=== {axis} (cohort={cohort_file}, positive={pos}) ===")
        feats = load_features(
            cohort, TCGA / "HiSeqV2.gz", TCGA / "HumanMethylation450.gz",
            meth_topk=10_000, dual_modality_only=True, positive_label=pos,
        )
        cnv, keep, cnv_genes = _cnv_aligned(feats.sample_ids)
        print(f"  RNA/meth n={len(feats.sample_ids)}; with CNV n={len(keep)}; CNV amplicon genes={len(cnv_genes)}")
        if len(keep) < 30:
            print(f"  too few samples with CNV ({len(keep)}); skipping axis")
            continue
        arrays = {"rna": feats.rna[keep], "meth": feats.meth[keep], "cnv": cnv}
        y = feats.y[keep]

        base, full = _five_fold(arrays, y)
        delta = full - base
        print(f"  rna+meth AUROC {base.mean():.3f}±{base.std():.3f} | "
              f"+cnv {full.mean():.3f}±{full.std():.3f} | delta {delta.mean():+.3f}")

        model = eval_ablation.fit_model(
            arrays, y, eval_ablation.FULL_SET, train_idx=np.arange(len(y)),
            latent_dim=LATENT_DIM, n_epochs=N_EPOCHS, seed=SEED)
        inputs = {m: torch.tensor(arrays[m], dtype=torch.float32) for m in eval_ablation.FULL_SET}
        attr = integrated_gradients(model, inputs, "cnv", steps=32)
        top5 = rank_genes_by_attribution(attr, cnv_genes)[:5]
        print(f"  CNV IG top-5: {[g for g, _ in top5]}")

        rows.append({
            "axis": axis, "n": int(len(y)),
            "base_mean": float(base.mean()), "base_std": float(base.std()),
            "full_mean": float(full.mean()), "full_std": float(full.std()),
            "delta": float(delta.mean()), "ig_top5": [g for g, _ in top5],
            "ig_has_erbb2": any(g in POLE_HER2_CNV for g, _ in top5),
        })

    if not rows:
        print("no axes ran; check data paths")
        return 1

    audit.emit("ablation_real_v0.2", "tcga", {r["axis"]: r["delta"] for r in rows})
    _write_audit_doc(rows)
    return 0


def _write_audit_doc(rows: list[dict]) -> None:
    table = "\n".join(
        f"| {r['axis']} | {r['n']} | {r['base_mean']:.3f} ± {r['base_std']:.3f} | "
        f"{r['full_mean']:.3f} ± {r['full_std']:.3f} | {r['delta']:+.3f} | "
        f"{', '.join(r['ig_top5'])} |"
        for r in rows
    )
    md = AUDIT / "ablation_v0.2.md"
    md.write_text(
        "# Modality ablation (TCGA real cohort) — v0.2\n\n"
        "Does adding the CNV branch (amplicon pole-masked) help over RNA+meth, per "
        "task axis? Same `MultiOmicsModel`, 5-fold StratifiedKFold (seed 42), "
        "v0.6 sizing. CNV is gene-level GISTIC2 restricted to the HER2 (ERBB2 17q12) "
        "and proliferation (MYC 8q24 / CCND1 11q13) amplicon masks, z-scored.\n\n"
        "| Axis | n | rna+meth AUROC | rna+meth+cnv AUROC | CNV delta | CNV IG top-5 |\n"
        "|---|---|---|---|---|---|\n" + table + "\n\n"
        "Honest reading: CNV is expected to help the amplicon-driven HER2 axis "
        "(ERBB2 should lead the CNV attribution) and add little to LumA-vs-LumB. "
        "Cross-cohort CNV (METABRIC SNP6) is a recorded harder problem (v0.3); this "
        "is the TCGA within-cohort ablation. A null delta on an axis is a valid result.\n",
    )
    print(f"\nwrote {md}")


if __name__ == "__main__":
    raise SystemExit(main())
