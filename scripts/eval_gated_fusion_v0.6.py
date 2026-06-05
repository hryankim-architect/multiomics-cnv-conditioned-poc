#!/usr/bin/env python3
"""v0.6 (1) — does modality-dropout fix the v0.5 gate collapse?

v0.5 (2): the input-conditioned softmax gate collapsed to CNV (~0.99) on both axes
(it learned a TCGA modality preference that does not transfer), so the gated model was
effectively CNV-only and *worsened* LumB. v0.6 (1) regularizes training with
**modality-dropout**: per epoch, with prob p, one modality's input is zeroed, forcing
each modality to stay usable and the gate to stay adaptive instead of collapsing.

Per axis (cross-cohort TCGA -> METABRIC, meth-free, class-weighted), four models:
RNA-only, RNA+CNV concat, RNA+CNV gated, RNA+CNV gated+dropout. If dropout works, the
gated+dropout CNV gate moves off ~0.99 and LumB's delta rises toward 0 (the gate stops
discarding RNA) while HER2 keeps its gain. If not, the collapse is more fundamental.

torch (run on his Mac). Reproduce:
  python scripts/build_metabric_cohort_v2.py && python scripts/eval_gated_fusion_v0.6.py
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
from dmoi_brca.external import align_to_train_genes, quantile_normalize_to_train  # noqa: E402
from dmoi_brca.features import load_features  # noqa: E402

from mocnv import audit, eval_ablation  # noqa: E402
from mocnv import cnv as cnvmod  # noqa: E402
from mocnv.priors import CNV_POLES  # noqa: E402

TCGA = DMOI / "data" / "tcga_brca"
METABRIC = DMOI / "data" / "metabric"
TCGA_CNV_GZ = REPO / "data" / "tcga_brca" / "Gistic2_CopyNumber_all_data_by_genes.gz"
METABRIC_CNA = REPO / "data" / "metabric" / "data_CNA.txt"
AUDIT = REPO / "audit"
AXES = [("HER2-vs-Luminal", "cohort_v4.tsv", "HER2"), ("LumA-vs-LumB", "cohort_v2.tsv", "LumB")]
AMPLICON_GENES = sorted({g for genes in CNV_POLES.values() for g in genes})
RNA_ONLY = ("rna",)
RNA_CNV = ("rna", "cnv")
LATENT_DIM = 128
N_EPOCHS = 120
SEED = 42
DROPOUT_P = 0.5


def _amplicon_cnv(matrix, want, *, tcga):
    samp = cnvmod.barcode_to_sample(matrix.sample_ids) if tcga else matrix.sample_ids
    idx = {s: i for i, s in enumerate(samp)}
    rows, keep = [], []
    for k, sid in enumerate(want):
        s = cnvmod.barcode_to_sample([sid])[0] if tcga else sid
        if s in idx:
            rows.append(idx[s])
            keep.append(k)
    aligned = cnvmod.align_to_genes(matrix.values[rows], matrix.gene_names, AMPLICON_GENES, fill_value=0.0)
    return cnvmod.harmonize_gene_level(aligned, method="zscore"), keep


def run_axis(cohort_file: str, positive: str) -> dict | None:
    tcga_cohort, metab_cohort = TCGA / cohort_file, METABRIC / cohort_file
    if not tcga_cohort.exists() or not metab_cohort.exists():
        sys.stderr.write(f"skip {positive}: missing cohort (build_metabric_cohort_v2.py?)\n")
        return None
    print(f"\n=== {positive} ({cohort_file}) ===")
    feats = load_features(tcga_cohort, TCGA / "HiSeqV2.gz", TCGA / "HumanMethylation450.gz",
                          meth_topk=10_000, dual_modality_only=True, positive_label=positive)
    tcga_cnv, tk = _amplicon_cnv(cnvmod.load_gistic2(TCGA_CNV_GZ), feats.sample_ids, tcga=True)
    tcga = {"rna": feats.rna[tk], "cnv": tcga_cnv}
    y_tcga = feats.y[tk]

    cohort = pd.read_csv(metab_cohort, sep="\t")
    cohort = cohort[cohort["has_rna"]].copy()
    want = cohort["sample_id"].astype(str).tolist()
    label = dict(zip(cohort["sample_id"].astype(str), (cohort["group"] == positive).astype(int), strict=False))
    rna_m = cnvmod.load_cbioportal_cna(METABRIC / "mrna_microarray.txt", sample_ids=set(want))
    cna_m = cnvmod.load_cbioportal_cna(METABRIC_CNA, sample_ids=set(want))
    common = [s for s in rna_m.sample_ids if s in set(cna_m.sample_ids)]
    rna_idx = {s: i for i, s in enumerate(rna_m.sample_ids)}
    rna_aligned = align_to_train_genes(rna_m.values[[rna_idx[s] for s in common]],
                                       rna_m.gene_names, feats.rna_features, fill_value=0.0)
    rna_qn = np.nan_to_num(quantile_normalize_to_train(rna_aligned, feats.rna), nan=0.0)
    cnv_m, mk = _amplicon_cnv(cna_m, common, tcga=False)
    common = [common[i] for i in mk]
    y_metab = np.array([label[s] for s in common], dtype=np.int64)
    metab = {"rna": rna_qn[mk].astype(np.float32), "cnv": cnv_m}
    print(f"  TCGA n={len(y_tcga)} | METABRIC n={len(y_metab)} ({positive}={int(y_metab.sum())})")

    all_t, all_m = np.arange(len(y_tcga)), np.arange(len(y_metab))
    metab_in = {m: torch.tensor(metab[m], dtype=torch.float32) for m in RNA_CNV}

    def _fit(mods, *, gated=False, md=0.0):
        m = eval_ablation.fit_model(tcga, y_tcga, mods, train_idx=all_t, latent_dim=LATENT_DIM,
                                    n_epochs=N_EPOCHS, seed=SEED, pos_weight=True, gated=gated,
                                    modality_dropout=md)
        return eval_ablation.auroc_of(m, metab, y_metab, mods, all_m), m

    auroc_rna, _ = _fit(RNA_ONLY)
    auroc_concat, _ = _fit(RNA_CNV)
    auroc_gated, m_gated = _fit(RNA_CNV, gated=True)
    auroc_drop, m_drop = _fit(RNA_CNV, gated=True, md=DROPOUT_P)
    gate_g = (m_gated.gate_weights(metab_in) or {}).get("cnv", float("nan"))
    gate_d = (m_drop.gate_weights(metab_in) or {}).get("cnv", float("nan"))
    print(f"  RNA {auroc_rna:.3f} | concat {auroc_concat:.3f} ({auroc_concat - auroc_rna:+.3f}) "
          f"| gated {auroc_gated:.3f} ({auroc_gated - auroc_rna:+.3f}, gate {gate_g:.3f}) "
          f"| gated+drop {auroc_drop:.3f} ({auroc_drop - auroc_rna:+.3f}, gate {gate_d:.3f})")
    return {"axis": positive, "rna": auroc_rna, "concat": auroc_concat, "gated": auroc_gated,
            "drop": auroc_drop, "gate_g": gate_g, "gate_d": gate_d}


def main() -> int:
    for p in (TCGA_CNV_GZ, METABRIC_CNA, METABRIC / "mrna_microarray.txt"):
        if not p.exists():
            sys.stderr.write(f"missing input: {p}\n")
            return 1
    AUDIT.mkdir(exist_ok=True)
    rows = [r for _, cf, pos in AXES if (r := run_axis(cf, pos)) is not None]
    if not rows:
        sys.stderr.write("no axes ran\n")
        return 1
    audit.emit("gated_fusion_v0.6", "tcga->metabric",
               {r["axis"]: {"d_concat": round(r["concat"] - r["rna"], 3),
                            "d_gated": round(r["gated"] - r["rna"], 3),
                            "d_drop": round(r["drop"] - r["rna"], 3),
                            "gate_drop": round(r["gate_d"], 3)} for r in rows})
    _write_doc(rows)
    return 0


def _write_doc(rows: list[dict]) -> None:
    table = "\n".join(  # deltas vs RNA-only; parentheses = mean CNV gate
        f"| {r['axis']} | {r['rna']:.3f} | {r['concat'] - r['rna']:+.3f} | "
        f"{r['gated'] - r['rna']:+.3f} ({r['gate_g']:.2f}) | {r['drop'] - r['rna']:+.3f} ({r['gate_d']:.2f}) |"
        for r in rows
    )
    md = AUDIT / "gated_fusion_v0.6.md"
    md.write_text(
        f"# Modality-dropout vs the gate collapse — v0.6 (1)\n\n"
        f"v0.5 (2) showed the softmax gate collapsed to CNV (~0.99) and worsened LumB. v0.6 (1) "
        f"adds **modality-dropout** training (p={DROPOUT_P}: per epoch, with this probability one "
        f"modality's input is zeroed) to keep each modality usable and the gate adaptive. "
        f"Cross-cohort, deltas vs RNA-only; values in parentheses are the mean CNV gate.\n\n"
        "| Axis | RNA-only | concat delta | gated delta (CNV gate) | gated+dropout delta (CNV gate) |\n"
        "|---|---|---|---|---|\n" + table + "\n\n"
        "Honest reading: if **gated+dropout** moves the CNV gate off ~0.99 and lifts LumB's delta "
        "toward >= 0 (the gate stops discarding the cross-platform-robust RNA) while HER2 keeps its "
        "gain, modality-dropout fixes the v0.5 collapse — the dilution was a regularization problem. "
        "If the gate still collapses and LumB stays negative, the failure is deeper than "
        "regularization (the gate has no held-out signal to learn cross-cohort modality trust), and "
        "plain concat remains the honest default. Either outcome is reported.\n",
    )
    print(f"\nwrote {md}")


if __name__ == "__main__":
    raise SystemExit(main())
