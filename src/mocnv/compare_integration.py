"""Label-free CNV feature-selection comparison: amplicon prior vs baselines.

The sibling `dmoi-brca-poc` showed that a *label-free* biological prior (Hallmark +
HM450-cis) beats statistical (top-variance) feature selection for PAM50 subtyping, on
the same footing as MOFA+/MoGCN. This ports that comparison to the CNV modality.

CNV is gene-level (GISTIC2) and amplicon-structured, so the prior here is the tight
amplicon-locus gene set (`priors.CNV_POLES`: ERBB2 17q12 + MYC/CCND1) — knowledge-based
and label-free, exactly comparable to a top-variance CNV selector. Per this repo's
thesis, the amplicon prior is expected to help *where the amplicon defines the axis*
(HER2) and not elsewhere (LumA-vs-LumB) — an axis-specific, honestly-reported result.

Selectors take only knowledge/variance, never `y`; the downstream classifier is the
only supervised step.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from mocnv.priors import CNV_POLES


def cnv_prior_universe(poles: Mapping[str, Sequence[str]] | None = None) -> set[str]:
    """Union of amplicon genes across CNV poles (default: all of CNV_POLES)."""
    src = poles if poles is not None else CNV_POLES
    universe: set[str] = set()
    for genes in src.values():
        universe.update(genes)
    return universe


def prior_cnv_indices(
    gene_names: Sequence[str],
    poles: Mapping[str, Sequence[str]] | None = None,
) -> list[int]:
    """Column indices of CNV genes in the amplicon-prior universe (label-free)."""
    universe = cnv_prior_universe(poles)
    return [i for i, g in enumerate(gene_names) if g in universe]


def topvar_indices(X: np.ndarray, k: int) -> list[int]:
    """Indices of the k highest-variance columns of X (label-free statistical baseline)."""
    if k <= 0:
        return []
    k = min(k, X.shape[1])
    var = np.nanvar(X, axis=0)
    return np.argsort(-var)[:k].tolist()


def topvar_within(X: np.ndarray, candidate_idx: Sequence[int], k: int) -> list[int]:
    """Top-k highest-variance columns among ``candidate_idx``, as ORIGINAL indices."""
    cand = list(candidate_idx)
    if k <= 0 or not cand:
        return []
    var = np.nanvar(X[:, cand], axis=0)
    order = np.argsort(-var)[: min(k, len(cand))]
    return [cand[i] for i in order]


def jaccard_index(a: Sequence[int], b: Sequence[int]) -> float:
    """Jaccard overlap |A∩B| / |A∪B| (0.0 if both empty)."""
    sa, sb = set(a), set(b)
    union = sa | sb
    return len(sa & sb) / len(union) if union else 0.0


def eval_selector(
    X: np.ndarray,
    y: np.ndarray,
    *,
    n_splits: int = 5,
    seed: int = 42,
) -> dict[str, float | None]:
    """Downstream evaluation of a selected-feature matrix (label-free selection).

    LogisticRegression + linear SVC weighted-F1 under stratified k-fold, plus AUROC
    when the task is binary (DMOI/mocnv headline metric), and Calinski-Harabasz /
    Davies-Bouldin of the standardized feature space vs labels. ``y`` is used only by
    the supervised downstream step / cluster metrics, never by the selectors.
    """
    from sklearn.linear_model import LogisticRegression  # noqa: PLC0415
    from sklearn.metrics import (  # noqa: PLC0415
        calinski_harabasz_score,
        davies_bouldin_score,
    )
    from sklearn.model_selection import StratifiedKFold, cross_val_score  # noqa: PLC0415
    from sklearn.pipeline import make_pipeline  # noqa: PLC0415
    from sklearn.preprocessing import StandardScaler  # noqa: PLC0415
    from sklearn.svm import SVC  # noqa: PLC0415

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    lr = make_pipeline(StandardScaler(),
                       LogisticRegression(max_iter=2000, class_weight="balanced"))
    sv = make_pipeline(StandardScaler(), SVC(kernel="linear", class_weight="balanced"))
    f1_lr = float(cross_val_score(lr, X, y, cv=cv, scoring="f1_weighted").mean())
    f1_svc = float(cross_val_score(sv, X, y, cv=cv, scoring="f1_weighted").mean())

    auroc: float | None = None
    classes = np.unique(y)
    if classes.size == 2:
        y01 = (y == classes[1]).astype(int)
        auroc = float(cross_val_score(lr, X, y01, cv=cv, scoring="roc_auc").mean())

    Xz = StandardScaler().fit_transform(X)
    return {
        "n_features": int(X.shape[1]),
        "auroc": auroc,
        "lr_weighted_f1": f1_lr,
        "svc_weighted_f1": f1_svc,
        "calinski_harabasz": float(calinski_harabasz_score(Xz, y)),
        "davies_bouldin": float(davies_bouldin_score(Xz, y)),
    }
