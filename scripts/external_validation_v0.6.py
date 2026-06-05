#!/usr/bin/env python3
"""v0.6 (3) — external validation of the strength->transfer law on a THIRD cohort.

v0.5 (1) found per-gene amplicon strength predicts cross-cohort transfer (TCGA->METABRIC,
pooled Spearman +0.836). Does it hold on an independent third cohort? This trains the
per-gene single-gene classifiers on TCGA (HER2 axis) and tests transfer on the **MBC
Project** (`brca_mbcproject_wagle_2017`, cBioPortal) — genome-wide gene-level CNA, HER2
label from `HER2_STATUS_REPORTED`. Get the data first:
  bash scripts/download_thirdcohort_probe.sh   (writes data/thirdcohort/brca_mbcproject_wagle_2017/)

Honest caveats (recorded): the MBC Project is metastatic (not primary) and HER2 is
patient-reported (not centrally assayed) -> a noisier test than METABRIC. If the law
still holds, that is strong external evidence; if weaker, these caveats explain it.

sklearn (no torch). Reproduce: python scripts/external_validation_v0.6.py
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
from mocnv.cohort import read_cbioportal_clinical  # noqa: E402
from mocnv.priors import CNV_POLES, POLE_HER2_CNV  # noqa: E402
from mocnv.strength import per_gene_strength_transfer, rank_spearman  # noqa: E402

GISTIC2 = REPO / "data" / "tcga_brca" / "Gistic2_CopyNumber_all_data_by_genes.gz"
THIRD = REPO / "data" / "thirdcohort" / "brca_mbcproject_wagle_2017"
AUDIT = REPO / "audit"
AMPLICON = sorted({g for genes in CNV_POLES.values() for g in genes})


def _align_z(matrix, want, *, tcga):
    samp = cnvmod.barcode_to_sample(matrix.sample_ids) if tcga else matrix.sample_ids
    idx = {s: i for i, s in enumerate(samp)}
    rows, keep = [], []
    for k, sid in enumerate(want):
        s = cnvmod.barcode_to_sample([sid])[0] if tcga else sid
        if s in idx:
            rows.append(idx[s])
            keep.append(k)
    aligned = cnvmod.align_to_genes(matrix.values[rows], matrix.gene_names, AMPLICON, fill_value=0.0)
    return cnvmod.harmonize_gene_level(aligned, method="zscore"), keep


def _third_cohort_her2():
    """MBC Project amplicon CNV (z-scored) + HER2 label (YES=1/NO=0) via sample->patient join."""
    pat = read_cbioportal_clinical(THIRD / "data_clinical_patient.txt")
    samp = read_cbioportal_clinical(THIRD / "data_clinical_sample.txt")
    her2 = dict(zip(pat["PATIENT_ID"].astype(str), pat["HER2_STATUS_REPORTED"].astype(str), strict=False))
    s2p = dict(zip(samp["SAMPLE_ID"].astype(str), samp["PATIENT_ID"].astype(str), strict=False))
    cna = cnvmod.load_cbioportal_cna(THIRD / "data_cna.txt")
    labels = {}
    for s in cna.sample_ids:
        status = her2.get(s2p.get(s, ""), "")
        if status in ("YES", "NO"):
            labels[s] = 1 if status == "YES" else 0
    want = [s for s in cna.sample_ids if s in labels]
    z, keep = _align_z(cna, want, tcga=False)
    y = np.array([labels[want[i]] for i in keep])
    return z, y


def main() -> int:
    if not GISTIC2.exists():
        sys.stderr.write(f"missing {GISTIC2}\n")
        return 1
    if not (THIRD / "data_cna.txt").exists():
        sys.stderr.write(f"missing {THIRD}/data_cna.txt\n  run scripts/download_thirdcohort_probe.sh first\n")
        return 1
    AUDIT.mkdir(exist_ok=True)

    # TCGA HER2 axis, amplicon CNV z-scored
    tco = pd.read_csv(DMOI / "data/tcga_brca/cohort_v4.tsv", sep="\t")
    tco = tco[tco["has_rna"] & tco["has_meth"]]
    z_t, keep = _align_z(cnvmod.load_gistic2(GISTIC2), tco["sample_id"].astype(str).tolist(), tcga=True)
    y_t = (tco["group"] == "HER2").astype(int).to_numpy()[keep]

    z_3, y_3 = _third_cohort_her2()
    print(f"TCGA n={len(y_t)} (HER2={int(y_t.sum())}) | MBCProject n={len(y_3)} (HER2+={int(y_3.sum())})")

    within, transfer = per_gene_strength_transfer(z_t, y_t, z_3, y_3)
    rho = rank_spearman(within, transfer)
    order = np.argsort(-within)
    print(f"\n  {'gene':9} {'pole':6} {'within(TCGA)':>12} {'transfer(MBC)':>14}")
    for j in order:
        g = AMPLICON[j]
        pole = "HER2" if g in POLE_HER2_CNV else "prolif"
        print(f"  {g:9} {pole:6} {within[j]:12.3f} {transfer[j]:14.3f}")
    print(f"\n  Spearman(strength, transfer to MBCProject) rho = {rho:+.3f}")
    print("  (METABRIC reference, HER2 axis, v0.5: rho = +0.877)")

    audit.emit("external_validation_v0.6", "tcga->mbcproject",
               {"rho_strength_transfer": round(rho, 3), "n_mbc": len(y_3), "her2_pos": int(y_3.sum())})
    _write_doc(within, transfer, rho, len(y_t), int(y_t.sum()), len(y_3), int(y_3.sum()))
    return 0


def _write_doc(within, transfer, rho, n_t, pos_t, n_3, pos_3) -> None:
    order = np.argsort(-within)
    rows = "\n".join(
        f"| {AMPLICON[j]} | {'HER2' if AMPLICON[j] in POLE_HER2_CNV else 'prolif'} | "
        f"{within[j]:.3f} | {transfer[j]:.3f} |"
        for j in order
    )
    md = AUDIT / "external_validation_v0.6.md"
    md.write_text(
        "# External validation of the strength->transfer law (third cohort) — v0.6 (3)\n\n"
        f"Train per-gene single-gene classifiers on TCGA (HER2 axis, n={n_t}, HER2={pos_t}); test "
        f"transfer on the **MBC Project** (`brca_mbcproject_wagle_2017`, n={n_3}, HER2+={pos_3}; "
        "genome-wide gene-level CNA, HER2 from `HER2_STATUS_REPORTED`). Does the v0.5 (1) "
        "strength->transfer law hold on an independent third cohort?\n\n"
        f"**Spearman(within-cohort strength, transfer to MBCProject) = {rho:+.3f}** "
        "(METABRIC reference, HER2 axis: +0.877).\n\n"
        "| gene | pole | within-cohort strength (TCGA AUROC) | transfer (MBCProject AUROC) |\n"
        "|---|---|---|---|\n" + rows + "\n\n"
        "Honest caveats: the MBC Project is **metastatic** (not primary) and HER2 is "
        "**patient-reported** (`HER2_STATUS_REPORTED`, not centrally assayed) — a noisier test "
        "than METABRIC. A positive rho here is independent-cohort evidence that the strength->"
        "transfer law generalizes; a weaker rho than METABRIC's +0.877 is consistent with the "
        "metastatic + self-reported-label noise, not necessarily a failure of the law. Reported, "
        "not hidden.\n",
    )
    print(f"\nwrote {md}")


if __name__ == "__main__":
    raise SystemExit(main())
