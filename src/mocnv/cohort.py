"""Cohort table: (sample_id, group, has_rna, has_meth, has_cnv).

Mirrors the dmoi cohort schema, extended with the CNV modality. v0.1 provides the
table builder + the three-modality intersection. The real TCGA/METABRIC
clinical-to-group assignment (PAM50 / HER2 status) is reused from dmoi when the
model is wired (v0.2); keeping it out here lets the data layer be tested on
sample-ID sets alone.
"""
from __future__ import annotations

from dataclasses import dataclass

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


def summarize(df: pd.DataFrame) -> CohortSummary:
    triple = df[df.has_rna & df.has_meth & df.has_cnv]
    dual = df[df.has_rna & df.has_meth]
    return CohortSummary(
        n_total=len(df),
        n_per_group={str(k): int(v) for k, v in df.group.value_counts().items()},
        n_dual=len(dual),
        n_triple=len(triple),
    )
