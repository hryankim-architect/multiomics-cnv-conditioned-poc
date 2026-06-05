"""Tests for the CNV data layer: GISTIC2 ingest, alignment, harmonization."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv import cnv  # noqa: E402

GISTIC2 = REPO / "data" / "tcga_brca" / "Gistic2_CopyNumber_all_data_by_genes.gz"


def test_barcode_to_sample_truncates_to_four_fields():
    assert cnv.barcode_to_sample(["TCGA-3C-AAAU-01A-11"]) == ["TCGA-3C-AAAU-01"]
    assert cnv.barcode_to_sample(["TCGA-3C-AAAU-01"]) == ["TCGA-3C-AAAU-01"]


def test_align_to_genes_reorders_and_fills():
    vals = np.array([[1.0, 2.0, 3.0]], dtype=np.float32)
    out = cnv.align_to_genes(vals, ["A", "B", "C"], ["C", "A", "MISSING"], fill_value=0.0)
    assert out.tolist() == [[3.0, 1.0, 0.0]]


def test_harmonize_zscore_standardizes_per_gene():
    rng = np.random.default_rng(0)
    vals = rng.normal(5.0, 3.0, size=(100, 4)).astype(np.float32)
    h = cnv.harmonize_gene_level(vals, method="zscore")
    assert np.allclose(h.mean(axis=0), 0.0, atol=1e-4)
    assert np.allclose(h.std(axis=0), 1.0, atol=1e-3)


def test_harmonize_none_is_passthrough():
    vals = np.array([[1.0, 2.0]], dtype=np.float32)
    assert np.array_equal(cnv.harmonize_gene_level(vals, method="none"), vals)


def test_harmonize_unknown_method_raises():
    with pytest.raises(ValueError, match="harmonize"):
        cnv.harmonize_gene_level(np.zeros((2, 2), dtype=np.float32), method="quantile")


@pytest.mark.skipif(not GISTIC2.exists(), reason="GISTIC2 download not present")
def test_load_gistic2_real_file_shape_and_anchors():
    m = cnv.load_gistic2(GISTIC2)
    assert m.n_samples > 500 and m.n_genes > 10000
    assert m.sample_ids[0].startswith("TCGA-")
    # the CNV pole anchors must exist in the real gene universe
    for anchor in ("ERBB2", "MYC", "CCND1"):
        assert anchor in m.gene_names
