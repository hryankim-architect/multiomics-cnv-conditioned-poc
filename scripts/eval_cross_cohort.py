#!/usr/bin/env python3
"""v0.3 cross-cohort CNV transfer: does the TCGA-trained CNV branch hold on METABRIC?

Trains the 3-modality model on TCGA (RNA + meth + GISTIC2 CNV, HER2 axis) and
scores it on METABRIC (RNA microarray quantile-normalized to TCGA + meth silenced
+ **SNP6** CNA, amplicon-masked). The ablation — METABRIC AUROC with vs without the
CNV branch — asks whether CNV *transfers across platforms* (GISTIC2 -> SNP6).

Honest expectation (stated in the README/scope): cross-platform CNV is harder than
RNA, so the cross-cohort CNV lift should be **weaker** than the TCGA within-cohort
result (v0.2: +0.125 on HER2). A null/negative cross-cohort delta is a valid,
on-brand result — it sharpens the recorded limit, it does not hide it.

REQUIRES (his Mac): `dmoi-brca-poc` alongside (RNA/meth loaders + data + QN/align),
TCGA GISTIC2 (this repo), and METABRIC `data_CNA.txt` (run scripts/download_metabric_cna.sh).

Reproduce:  python scripts/eval_cross_cohort.py
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
sys.path.insert(0, str(DMOI / "src"))

import torch  # noqa: E402
from dmoi_brca.external import (  # noqa: E402
    align_to_train_genes,
    make_silenced_meth,
    quantile_normalize_to_train,
)
from dmoi_brca.features import load_features  # noqa: E402

from mocnv import audit, eval_ablation  # noqa: E402
from mocnv import cnv as cnvmod  # noqa: E402
from mocnv.attribution import integrated_gradients, rank_genes_by_attribution  # noqa: E402
from mocnv.priors import CNV_POLES, POLE_HER2_CNV  # noqa: E402

TCGA = DMOI / "data" / "tcga_brca"
METABRIC = DMOI / "data" / "metabric"
TCGA_CNV_GZ = REPO / "data" / "tcga_brca" / "Gistic2_CopyNumber_all_data_by_genes.gz"
METABRIC_CNA = REPO / "data" / "metabric" / "data_CNA.txt"
AUDIT = REPO / "audit"

POSITIVE = "HER2"
AMPLICON_GENES = sorted({g for genes in CNV_POLES.values() for g in genes})
LATENT_DIM = 128
N_EPOCHS = 120
SEED = 42


def _amplicon_cnv(matrix: cnvmod.CNVMatrix, want_samples: list[str], *, tcga: bool) -> tuple[np.ndarray, list[str]]:
    """Amplicon-masked, z-scored CNV aligned to `want_samples` (order preserved)."""
    samp = cnvmod.barcode_to_sample(matrix.sample_ids) if tcga else matrix.sample_ids
    idx = {s: i for i, s in enumerate(samp)}
    rows, keep = [], []
    for k, sid in enumerate(want_samples):
        s = cnvmod.barcode_to_sample([sid])[0] if tcga else sid
        if s in idx:
            rows.append(idx[s])
            keep.append(k)
    sub = matrix.values[rows]
    aligned = cnvmod.align_to_genes(sub, matrix.gene_names, AMPLICON_GENES, fill_value=0.0)
    return cnvmod.harmonize_gene_level(aligned, method="zscore"), keep


def main() -> int:
    for p in (TCGA / "cohort_v4.tsv", METABRIC / "cohort_v4.tsv", METABRIC / "mrna_microarray.txt",
              TCGA_CNV_GZ, METABRIC_CNA):
        if not p.exists():
            sys.stderr.write(f"missing input: {p}\n"
                             "  (METABRIC CNA? run scripts/download_metabric_cna.sh)\n")
            return 1
    AUDIT.mkdir(exist_ok=True)

    # --- TCGA train (RNA + meth + GISTIC2 CNV), HER2 axis ---
    print("--- TCGA train (cohort_v4, HER2) ---")
    feats = load_features(TCGA / "cohort_v4.tsv", TCGA / "HiSeqV2.gz", TCGA / "HumanMethylation450.gz",
                          meth_topk=10_000, dual_modality_only=True, positive_label=POSITIVE)
    tcga_cnv, tk = _amplicon_cnv(cnvmod.load_gistic2(TCGA_CNV_GZ), feats.sample_ids, tcga=True)
    tcga = {"rna": feats.rna[tk], "meth": feats.meth[tk], "cnv": tcga_cnv}
    y_tcga = feats.y[tk]
    print(f"  TCGA n={len(y_tcga)} (HER2={int(y_tcga.sum())}); rna={tcga['rna'].shape[1]} meth={tcga['meth'].shape[1]} cnv={len(AMPLICON_GENES)}")

    # --- METABRIC score (RNA QN->TCGA + meth silenced + SNP6 CNA amplicon-masked) ---
    print("--- METABRIC (cohort_v4) ---")
    cohort = pd.read_csv(METABRIC / "cohort_v4.tsv", sep="\t")
    cohort = cohort[cohort["has_rna"]].copy()
    want = cohort["sample_id"].astype(str).tolist()
    label = dict(zip(cohort["sample_id"].astype(str), (cohort["group"] == POSITIVE).astype(int), strict=False))

    rna_m = cnvmod.load_cbioportal_cna(METABRIC / "mrna_microarray.txt", sample_ids=set(want))
    cna_m = cnvmod.load_cbioportal_cna(METABRIC_CNA, sample_ids=set(want))
    common = [s for s in rna_m.sample_ids if s in set(cna_m.sample_ids)]
    if len(common) < 30:
        sys.stderr.write(f"too few METABRIC samples with RNA+CNV: {len(common)}\n")
        return 1

    rna_idx = {s: i for i, s in enumerate(rna_m.sample_ids)}
    rna_sub = rna_m.values[[rna_idx[s] for s in common]]
    rna_aligned = align_to_train_genes(rna_sub, rna_m.gene_names, feats.rna_features, fill_value=0.0)
    rna_qn = quantile_normalize_to_train(rna_aligned, feats.rna)
    meth_sil = make_silenced_meth(len(common), feats.meth.shape[1])
    cnv_m, mk = _amplicon_cnv(cna_m, common, tcga=False)
    common = [common[i] for i in mk]                       # keep only samples present in CNV
    rna_qn, meth_sil = rna_qn[mk], meth_sil[mk]
    y_metab = np.array([label[s] for s in common], dtype=np.int64)
    metab = {"rna": rna_qn.astype(np.float32), "meth": meth_sil, "cnv": cnv_m}
    print(f"  METABRIC n={len(common)} (HER2={int(y_metab.sum())}); RNA QN->TCGA, meth silenced, CNV=SNP6 amplicon")

    # --- cross-cohort ablation: train TCGA, score METABRIC ---
    all_t = np.arange(len(y_tcga))
    all_m = np.arange(len(y_metab))
    base_model = eval_ablation.fit_model(tcga, y_tcga, eval_ablation.BASELINE_SET, train_idx=all_t,
                                         latent_dim=LATENT_DIM, n_epochs=N_EPOCHS, seed=SEED)
    full_model = eval_ablation.fit_model(tcga, y_tcga, eval_ablation.FULL_SET, train_idx=all_t,
                                         latent_dim=LATENT_DIM, n_epochs=N_EPOCHS, seed=SEED)
    auroc_base = eval_ablation.auroc_of(base_model, metab, y_metab, eval_ablation.BASELINE_SET, all_m)
    auroc_full = eval_ablation.auroc_of(full_model, metab, y_metab, eval_ablation.FULL_SET, all_m)
    delta = auroc_full - auroc_base
    print(f"\n  METABRIC (cross-cohort) rna+meth AUROC {auroc_base:.3f} | +cnv {auroc_full:.3f} | delta {delta:+.3f}")

    inputs = {m: torch.tensor(metab[m], dtype=torch.float32) for m in eval_ablation.FULL_SET}
    attr = integrated_gradients(full_model, inputs, "cnv", steps=32)
    top5 = rank_genes_by_attribution(attr, AMPLICON_GENES)[:5]
    erbb2_leads = any(g in POLE_HER2_CNV for g, _ in top5)
    print(f"  METABRIC CNV IG top-5: {[g for g, _ in top5]}  (ERBB2-amplicon present: {erbb2_leads})")

    audit.emit("cross_cohort_v0.3", "tcga->metabric",
               {"auroc_base": auroc_base, "auroc_full": auroc_full, "delta": delta})
    _write_doc(len(y_tcga), int(y_tcga.sum()), len(y_metab), int(y_metab.sum()),
               auroc_base, auroc_full, delta, [g for g, _ in top5])
    return 0


def _write_doc(n_t, her2_t, n_m, her2_m, base, full, delta, ig_top5) -> None:
    md = AUDIT / "cross_cohort_v0.3.md"
    md.write_text(
        "# Cross-cohort CNV transfer (TCGA -> METABRIC) — v0.3\n\n"
        f"Train TCGA cohort_v4 (n={n_t}, HER2={her2_t}; RNA+meth+GISTIC2 CNV), score "
        f"METABRIC cohort_v4 (n={n_m}, HER2={her2_m}; RNA QN->TCGA + meth silenced + "
        "**SNP6** CNA, amplicon-masked + z-scored). Does the CNV branch transfer across "
        "platforms?\n\n"
        "| Setting | METABRIC AUROC |\n|---|---|\n"
        f"| rna + meth(silenced) | {base:.3f} |\n"
        f"| + CNV (SNP6, amplicon) | {full:.3f} |\n"
        f"| **cross-cohort CNV delta** | **{delta:+.3f}** |\n\n"
        f"METABRIC CNV IG top-5: {', '.join(ig_top5)}\n\n"
        "Honest reading: TCGA GISTIC2 and METABRIC SNP6 are different platforms; gene-level "
        "alignment + per-gene z-scoring is the cross-platform bridge, not a validated "
        "normalization. The within-cohort TCGA result was +0.125 (v0.2); a smaller or null "
        "cross-cohort delta here is expected and reported, not hidden — it is the recorded "
        "limit of cross-platform CNV, sharpened into a measurement.\n",
    )
    print(f"\nwrote {md}")


if __name__ == "__main__":
    raise SystemExit(main())
