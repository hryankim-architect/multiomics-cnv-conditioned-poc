"""Tests for the synthetic 3-modality fixtures + the planted ablation hypothesis.

These pin the repo's central claim on synthetic data: CNV separates an
amplicon-driven (HER2) axis and is flat on a non-amplicon (Luminal) axis, where
RNA carries the signal instead.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv import synth  # noqa: E402
from mocnv.priors import POLE_HER2_CNV  # noqa: E402


def _cnv_amplicon_auroc(s: synth.SyntheticMultiomics) -> float:
    cols = [s.cnv_genes.index(g) for g in POLE_HER2_CNV]
    score = s.cnv[:, cols].mean(axis=1)
    return float(roc_auc_score(s.y, score))


def _best_rna_auroc(s: synth.SyntheticMultiomics) -> float:
    return max(float(roc_auc_score(s.y, s.rna[:, j])) for j in range(10))


def test_her2_axis_cnv_is_informative():
    s = synth.generate("HER2", seed=0, n=300)
    assert _cnv_amplicon_auroc(s) > 0.75          # CNV amplicon separates the HER2 axis


def test_luminal_axis_cnv_is_flat_but_rna_carries_signal():
    s = synth.generate("Luminal", seed=0, n=300)
    assert 0.38 <= _cnv_amplicon_auroc(s) <= 0.62  # CNV ~ chance on the Luminal axis
    assert _best_rna_auroc(s) > 0.70               # RNA carries it


def test_generate_is_deterministic():
    a = synth.generate("HER2", seed=1, n=120)
    b = synth.generate("HER2", seed=1, n=120)
    assert np.array_equal(a.cnv, b.cnv)
    assert np.array_equal(a.y, b.y)
    assert a.cnv_genes == b.cnv_genes


def test_shapes_and_modalities_consistent():
    s = synth.generate("HER2", seed=2, n=80)
    assert s.cnv.shape == (80, len(s.cnv_genes))
    assert s.rna.shape == (80, len(s.rna_genes))
    assert s.meth.shape == (80, len(s.meth_probes))
    assert s.y.shape == (80,)
    # all HER2 amplicon genes are in the synthetic CNV universe
    assert all(g in s.cnv_genes for g in POLE_HER2_CNV)


def test_invalid_axis_raises():
    with pytest.raises(ValueError, match="axis"):
        synth.generate("Basal")
