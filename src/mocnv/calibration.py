"""Probability calibration metrics for the cross-cohort transfer (v0.5 (3)).

AUROC measures ranking; calibration measures whether the predicted probabilities are
trustworthy. A TCGA-trained model can rank METABRIC well (good AUROC) yet be badly
calibrated cross-platform. Brier score + expected calibration error (ECE) make that
explicit. Pure numpy (no torch), so the metrics are unit-tested directly.
"""
from __future__ import annotations

import numpy as np


def brier_score(y: np.ndarray, p: np.ndarray) -> float:
    """Mean squared error of probabilistic predictions (lower is better)."""
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    return float(np.mean((p - y) ** 2))


def reliability_bins(
    y: np.ndarray, p: np.ndarray, *, n_bins: int = 10
) -> list[tuple[float, float, int]]:
    """Equal-width reliability bins: (mean predicted prob, observed frequency, count)."""
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    out: list[tuple[float, float, int]] = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        mask = (p >= lo) & (p < hi) if hi < 1.0 else (p >= lo) & (p <= hi)
        n = int(mask.sum())
        if n:
            out.append((float(p[mask].mean()), float(y[mask].mean()), n))
    return out


def expected_calibration_error(y: np.ndarray, p: np.ndarray, *, n_bins: int = 10) -> float:
    """ECE: count-weighted mean gap between confidence and accuracy across bins."""
    n_total = len(p)
    if n_total == 0:
        return float("nan")
    return float(sum(
        (count / n_total) * abs(conf - acc)
        for conf, acc, count in reliability_bins(y, p, n_bins=n_bins)
    ))
