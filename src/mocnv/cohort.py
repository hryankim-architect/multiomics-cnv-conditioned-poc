"""Cohort table: (sample_id, group, has_rna, has_meth, has_cnv).

Mirrors the dmoi cohort schema, extended with the CNV modality. v0.1 provides the
table builder + the three-modality intersection. The real TCGA/METABRIC
clinical-to-group assignment (PAM50 / HER2 status) is reused from dmoi when the
model is wired (v0.2); keeping it out here lets the data layer be tested on
sample-ID sets alone.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class CohortSummary:
    n_total: int
    n_per_group: dict[str, int]
    n_dual: int      # rna + meth
    n_triple: int    # rna + meth + cnv


def build_cohort_table(
    labels: dict[str, str | None],
    *,
    rna_ids: set[str] | list[str],
    meth_ids: set[str] | list[str],
    cnv_ids: set[str] | list[str],
) -> pd.DataFrame:
    """Build (sample_id, group, has_rna, has_meth, has_cnv) from labels + ID sets.

    ``labels`` maps sample_id -> group (None drops the sample). The has_* flags
    record which modalities measured each sample.
    """
    rna_ids, meth_ids, cnv_ids = set(rna_ids), set(meth_ids), set(cnv_ids)
    rows = [
        {
            "sample_id": sid,
            "group": group,
            "has_rna": sid in rna_ids,
            "has_meth": sid in meth_ids,
            "has_cnv": sid in cnv_ids,
        }
        for sid, group in labels.items()
        if group is not None
    ]
    return pd.DataFrame(rows, columns=["sample_id", "group", "has_rna", "has_meth", "has_cnv"])


def read_cbioportal_clinical(path: str | Path) -> pd.DataFrame:
    """Read a cBioPortal clinical_*.txt, skipping leading ``#``-prefixed metadata.

    cBioPortal clinical files prefix several metadata header lines with ``#``; the
    first non-``#`` line is the real column header (``PATIENT_ID``, ...). We count
    and skip exactly those leading comment lines (not ``comment='#'``, which would
    also strip ``#`` inside data cells).
    """
    skip = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=skip)


def build_metabric_cohort_v2(
    clinical_patient: str | Path,
    cohort_v4: str | Path,
    *,
    subtype_col: str = "CLAUDIN_SUBTYPE",
    id_col: str = "PATIENT_ID",
    keep: tuple[str, ...] = ("LumA", "LumB"),
) -> pd.DataFrame:
    """Split the METABRIC cohort_v4 Luminal lump into LumA/LumB via clinical PAM50.

    METABRIC ``cohort_v4.tsv`` carries Luminal as a single group; the PAM50 +
    Claudin-low call (``CLAUDIN_SUBTYPE`` in ``clinical_patient.txt``) refines it.
    We join each cohort_v4 ``sample_id`` to its patient-level subtype and keep the
    requested classes, inheriting the cohort's ``has_rna``/``has_meth`` flags so
    the modality availability matches the rest of the pipeline.

    Returns columns ``(sample_id, group, has_rna, has_meth)``. Samples whose
    subtype is not in ``keep`` (HER2, Basal, Normal, blank, ...) are dropped.
    """
    clin = read_cbioportal_clinical(clinical_patient)
    subtype = dict(zip(clin[id_col].astype(str), clin[subtype_col].astype(str), strict=False))
    v4 = pd.read_csv(cohort_v4, sep="\t")
    keep_set = set(keep)
    rows = [
        {"sample_id": r["sample_id"], "group": g,
         "has_rna": r["has_rna"], "has_meth": r["has_meth"]}
        for _, r in v4.iterrows()
        if (g := subtype.get(str(r["sample_id"]), "")) in keep_set
    ]
    return pd.DataFrame(rows, columns=["sample_id", "group", "has_rna", "has_meth"])


def summarize(df: pd.DataFrame) -> CohortSummary:
    triple = df[df.has_rna & df.has_meth & df.has_cnv]
    dual = df[df.has_rna & df.has_meth]
    return CohortSummary(
        n_total=len(df),
        n_per_group={str(k): int(v) for k, v in df.group.value_counts().items()},
        n_dual=len(dual),
        n_triple=len(triple),
    )
