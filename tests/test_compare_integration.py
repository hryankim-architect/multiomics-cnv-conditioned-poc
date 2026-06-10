from __future__ import annotations

import numpy as np

from mocnv.compare_integration import (
    cnv_prior_universe,
    eval_selector,
    jaccard_index,
    prior_cnv_indices,
    topvar_indices,
    topvar_within,
)

FAKE = {"HER2": ("ERBB2", "GRB7"), "Prolif": ("MYC", "CCND1")}


def test_cnv_prior_universe_unions():
    assert cnv_prior_universe(FAKE) == {"ERBB2", "GRB7", "MYC", "CCND1"}


def test_prior_cnv_indices_label_free():
    genes = ["TP53", "ERBB2", "FOO", "MYC"]
    assert prior_cnv_indices(genes, FAKE) == [1, 3]


def test_real_priors_present():
    # The shipped CNV_POLES universe is non-empty and includes the keystone anchors.
    u = cnv_prior_universe()
    assert {"ERBB2", "MYC", "CCND1"} <= u


def test_topvar_helpers():
    X = np.array([[0.0, 0.0, 0.0], [0.0, 5.0, 1.0], [0.0, -5.0, -1.0]])
    assert topvar_indices(X, 2) == [1, 2]
    assert topvar_within(X, [0, 2], 1) == [2]
    assert jaccard_index([1, 2, 3], [2, 3, 4]) == 0.5


def test_eval_selector_binary_has_auroc():
    rng = np.random.default_rng(0)
    n = 60
    y = np.array(["a"] * (n // 2) + ["b"] * (n // 2))
    sig = np.where(y == "b", 3.0, -3.0)[:, None] + rng.normal(0, 0.3, (n, 3))
    out = eval_selector(sig, y, n_splits=5, seed=0)
    assert out["n_features"] == 3
    assert out["auroc"] is not None and out["auroc"] > 0.9
    assert out["lr_weighted_f1"] > 0.9


def test_eval_selector_multiclass_auroc_none():
    rng = np.random.default_rng(1)
    y = np.array(["a", "b", "c"] * 20)
    x = rng.normal(0, 1, (60, 3))
    out = eval_selector(x, y, n_splits=3, seed=1)
    assert out["auroc"] is None  # AUROC only computed for binary tasks
