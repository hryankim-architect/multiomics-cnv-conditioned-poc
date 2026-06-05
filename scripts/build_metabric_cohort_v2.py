#!/usr/bin/env python3
"""Build METABRIC `cohort_v2.tsv` (LumA-vs-LumB) for the v0.4 proliferation axis.

METABRIC `cohort_v4.tsv` carries Luminal as a single lump; the v0.4 cross-cohort
check needs it split into LumA / LumB (the proliferation / MYC-8q24 axis, mirroring
the TCGA `cohort_v2.tsv` that v0.2 used). The split comes from the PAM50 +
Claudin-low call (`CLAUDIN_SUBTYPE`) in the cBioPortal clinical_patient.txt.

Reads from the dmoi-brca-poc data layer (RNA/meth/clinical are reused from there):
  <dmoi>/data/metabric/clinical_patient.txt   (CLAUDIN_SUBTYPE per PATIENT_ID)
  <dmoi>/data/metabric/cohort_v4.tsv          (sample_id, group, has_rna, has_meth)
Writes:
  <dmoi>/data/metabric/cohort_v2.tsv          (sample_id, group in {LumA,LumB}, has_rna, has_meth)

Override the dmoi location with the MOCNV_DMOI env var. Reproduce:
  python scripts/build_metabric_cohort_v2.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
DMOI = Path(os.environ.get("MOCNV_DMOI", str(REPO.parent / "dmoi-brca-poc")))

from mocnv.cohort import build_metabric_cohort_v2  # noqa: E402

METABRIC = DMOI / "data" / "metabric"
CLINICAL = METABRIC / "clinical_patient.txt"
COHORT_V4 = METABRIC / "cohort_v4.tsv"
OUT = METABRIC / "cohort_v2.tsv"


def main() -> int:
    for p in (CLINICAL, COHORT_V4):
        if not p.exists():
            sys.stderr.write(f"missing input: {p}\n  (expected in the dmoi-brca-poc data layer)\n")
            return 1

    df = build_metabric_cohort_v2(CLINICAL, COHORT_V4)
    if df.empty:
        sys.stderr.write("no LumA/LumB samples found; check CLAUDIN_SUBTYPE column\n")
        return 1

    df.to_csv(OUT, sep="\t", index=False)
    counts = df["group"].value_counts().to_dict()
    rna = int(df["has_rna"].astype(str).isin({"True", "1", "true"}).sum())
    print(f"wrote {OUT}")
    print(f"  n={len(df)}  groups={counts}  has_rna={rna}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
