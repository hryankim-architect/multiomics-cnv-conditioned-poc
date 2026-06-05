"""Tests for the calibration metrics (v0.5 (3); pure numpy, no torch)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv.calibration import (  # noqa: E402
    brier_score,
    expected_calibration_error,
    reliability_bins,
)


def test_brier_perfect_and_worst():
    y = np.array([1, 0, 1, 0])
    assert brier_score(y, y.astype(float)) == 0.0       # perfect predictions
    assert brier_score(y, 1.0 - y) == 1.0               # confidently wrong


def test_ece_perfectly_calibrated_is_low():
    rng = np.random.default_rng(0)
    p = rng.random(2000)
    y = (rng.random(2000) < p).astype(int)              # y drawn with prob p -> calibrated
    assert expected_calibration_error(y, p, n_bins=10) < 0.08


def test_ece_overconfident_is_high():
    y = np.zeros(100, dtype=int)                         # all negative
    p = np.full(100, 0.9)                                # predicted 0.9 -> badly miscalibrated
    assert expected_calibration_error(y, p) > 0.7


def test_reliability_bins_cover_all_points():
    y = np.array([0, 1, 1, 0])
    p = np.array([0.1, 0.9, 0.8, 0.2])
    bins = reliability_bins(y, p, n_bins=2)
    assert sum(count for *_, count in bins) == 4
