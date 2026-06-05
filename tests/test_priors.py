"""Tests for CNV pole masks (amplicon-locus gene sets)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv import priors  # noqa: E402


def test_cnv_poles_contain_anchor_genes():
    assert "ERBB2" in priors.POLE_HER2_CNV
    assert "MYC" in priors.POLE_PROLIFERATION_CNV
    assert "CCND1" in priors.POLE_PROLIFERATION_CNV


def test_cnv_pole_mask_selects_only_pole_genes():
    universe = ["FOO", "ERBB2", "BAR", "GRB7", "MYC"]
    m = priors.cnv_pole_mask("HER2", universe)
    assert m.dtype == bool
    assert m.tolist() == [False, True, False, True, False]


def test_cnv_pole_mask_unknown_pole_raises():
    with pytest.raises(KeyError, match="pole"):
        priors.cnv_pole_mask("Basal", ["ERBB2"])


def test_mask_coverage_reports_present_and_anchors():
    universe = [*priors.POLE_HER2_CNV, "X", "Y"]
    cov = priors.mask_coverage("HER2", universe)
    assert cov["n_present"] == cov["n_in_set"]
    assert cov["missing"] == []
    assert "ERBB2" in cov["anchors_present"]
