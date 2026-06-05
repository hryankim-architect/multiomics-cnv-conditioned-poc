"""Tests for the cohort table builder (3-modality intersection)."""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv import cohort  # noqa: E402


def test_build_cohort_table_drops_none_and_flags_modalities():
    labels = {"s1": "HER2", "s2": "Luminal", "s3": "HER2", "s4": None}
    df = cohort.build_cohort_table(
        labels, rna_ids={"s1", "s2", "s3"}, meth_ids={"s1", "s2"}, cnv_ids={"s1", "s3"}
    )
    assert len(df) == 3  # s4 (None) dropped
    assert list(df.columns) == ["sample_id", "group", "has_rna", "has_meth", "has_cnv"]
    row1 = df[df.sample_id == "s1"].iloc[0]
    assert bool(row1.has_rna) and bool(row1.has_meth) and bool(row1.has_cnv)


def test_summarize_counts_dual_and_triple():
    labels = {"s1": "HER2", "s2": "Luminal", "s3": "HER2"}
    df = cohort.build_cohort_table(
        labels, rna_ids={"s1", "s2", "s3"}, meth_ids={"s1", "s2"}, cnv_ids={"s1", "s3"}
    )
    s = cohort.summarize(df)
    assert s.n_total == 3
    assert s.n_per_group["HER2"] == 2
    assert s.n_dual == 2     # s1, s2 have rna + meth
    assert s.n_triple == 1   # only s1 has all three
