"""Tests for the per-gene strength/transfer helpers (v0.5 (1); sklearn, no torch)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv.strength import per_gene_strength_transfer, rank_spearman  # noqa: E402


def test_rank_spearman_monotone():
    assert rank_spearman([1, 2, 3, 4], [10, 20, 30, 40]) == 1.0
    assert rank_spearman([1, 2, 3, 4], [40, 30, 20, 10]) == -1.0


def test_rank_spearman_constant_is_nan():
    assert np.isnan(rank_spearman([1, 1, 1], [1, 2, 3]))


def test_per_gene_strong_transfers_noise_does_not():
    n = 300

    def cohort(seed):
        r = np.random.default_rng(seed)
        y = (r.random(n) < 0.4).astype(int)
        strong = y + r.normal(0, 0.3, n)     # cleanly separates the classes
        noise = r.normal(0, 1, n)            # unrelated to y
        return np.column_stack([strong, noise]).astype(np.float32), y

    cnv_tr, y_tr = cohort(1)
    cnv_te, y_te = cohort(2)
    within, transfer = per_gene_strength_transfer(cnv_tr, y_tr, cnv_te, y_te)

    assert within[0] > 0.85 and transfer[0] > 0.85    # strong gene: high within + transfers
    assert abs(within[1] - 0.5) < 0.15                # noise gene: ~chance within
    assert abs(transfer[1] - 0.5) < 0.15              # noise gene: ~chance transfer
    assert within[0] > within[1] and transfer[0] > transfer[1]
