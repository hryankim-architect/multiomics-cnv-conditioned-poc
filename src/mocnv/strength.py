"""Per-amplicon-gene strength vs cross-cohort transfer (v0.5 (1)).

The v0.4 axis-level finding ("CNV helps where the amplicon defines the axis") made
amplicon strength a *binary* axis property. v0.5 (1) turns it into a *continuous*
per-gene variable: for each amplicon gene, its within-cohort discriminative strength
(single-gene AUROC on the train cohort) vs its cross-cohort transfer (a 1-feature
classifier trained on the train cohort and scored on the held-out cohort). The
hypothesis: strength predicts transfer — strongly-amplified, discriminative genes
carry across platforms, weak ones do not.

Single-gene classifiers are one feature, so this is sklearn (no torch).
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


def rank_spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation = Pearson on ranks (kept scipy-free)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def per_gene_strength_transfer(
    cnv_train: np.ndarray, y_train: np.ndarray,
    cnv_test: np.ndarray, y_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """For each gene column, return (within-cohort strength, cross-cohort transfer) AUROC.

    ``within`` is the single-gene AUROC on the train cohort (the raw z-scored CNV used
    directly as the score: higher copy number -> positive class). ``transfer`` trains a
    1-feature logistic on the train cohort and scores the test cohort, so the learned
    direction must hold across platforms for the transfer AUROC to stay above chance.
    """
    n_genes = cnv_train.shape[1]
    within = np.empty(n_genes, dtype=float)
    transfer = np.empty(n_genes, dtype=float)
    for j in range(n_genes):
        within[j] = roc_auc_score(y_train, cnv_train[:, j])
        clf = LogisticRegression(max_iter=1000).fit(cnv_train[:, [j]], y_train)
        transfer[j] = roc_auc_score(y_test, clf.predict_proba(cnv_test[:, [j]])[:, 1])
    return within, transfer
