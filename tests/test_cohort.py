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


# --- v0.4: METABRIC LumA/LumB cohort builder (PAM50 split of the Luminal lump) ---

_CLINICAL = (
    "#Patient Identifier\tPam50 + Claudin-low subtype\n"  # cBioPortal metadata header lines
    "#patient id\tsubtype\n"
    "#STRING\tSTRING\n"
    "#1\t1\n"
    "PATIENT_ID\tCLAUDIN_SUBTYPE\n"                        # the real column header (no '#')
    "MB-0001\tLumA\n"
    "MB-0002\tLumB\n"
    "MB-0003\tHer2\n"
    "MB-0004\tLumA\n"
    "MB-0005\t\n"                                          # blank subtype
)
_COHORT_V4 = (
    "sample_id\tgroup\thas_rna\thas_meth\n"
    "MB-0001\tLuminal\tTrue\tFalse\n"
    "MB-0002\tLuminal\tTrue\tTrue\n"
    "MB-0003\tHER2\tTrue\tFalse\n"
    "MB-0004\tLuminal\tFalse\tTrue\n"
    "MB-9999\tLuminal\tTrue\tTrue\n"                       # not in clinical -> dropped
)


def test_read_cbioportal_clinical_skips_comment_header(tmp_path):
    p = tmp_path / "clinical_patient.txt"
    p.write_text(_CLINICAL)
    df = cohort.read_cbioportal_clinical(p)
    assert list(df.columns) == ["PATIENT_ID", "CLAUDIN_SUBTYPE"]
    assert len(df) == 5  # five data rows, four '#' lines skipped


def test_build_metabric_cohort_v2_splits_luminal(tmp_path):
    clin = tmp_path / "clinical_patient.txt"
    clin.write_text(_CLINICAL)
    v4 = tmp_path / "cohort_v4.tsv"
    v4.write_text(_COHORT_V4)

    df = cohort.build_metabric_cohort_v2(clin, v4)
    assert list(df.columns) == ["sample_id", "group", "has_rna", "has_meth"]
    assert set(df["sample_id"]) == {"MB-0001", "MB-0002", "MB-0004"}  # LumA/LumB only
    assert df["group"].value_counts().to_dict() == {"LumA": 2, "LumB": 1}
    assert "MB-0003" not in set(df["sample_id"])   # Her2 dropped
    assert "MB-9999" not in set(df["sample_id"])   # absent-from-clinical dropped
    # has_rna inherited from cohort_v4 (MB-0004 was False)
    assert not bool(df[df.sample_id == "MB-0004"].iloc[0]["has_rna"])
