"""Modality-ablation harness: does CNV add signal over RNA+meth, per axis?

Trains the same ``MultiOmicsModel`` on two modality sets — {rna, meth} and
{rna, meth, cnv} — on identical splits and reports the held-out AUROC delta. The
honest expectation: a positive delta on amplicon-driven axes (HER2) and ~0 on
others (Luminal). A null delta is a valid result.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from torch import nn

from mocnv.model import MultiOmicsModel

BASELINE_SET = ("rna", "meth")
FULL_SET = ("rna", "meth", "cnv")
RNA_ONLY = ("rna",)
CNV_ONLY = ("cnv",)


@dataclass
class AblationResult:
    axis: str
    auroc_baseline: float        # rna + meth
    auroc_full: float            # rna + meth + cnv
    delta: float                 # full - baseline (CNV's contribution)


def _tensor(a: np.ndarray) -> torch.Tensor:
    return torch.tensor(a, dtype=torch.float32)


def fit_model(
    modality_arrays: dict[str, np.ndarray],
    y: np.ndarray,
    modalities: tuple[str, ...],
    *,
    train_idx: np.ndarray,
    latent_dim: int = 32,
    n_epochs: int = 60,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    seed: int = 0,
    pos_weight: bool = False,
    gated: bool = False,
    modality_dropout: float = 0.0,
) -> MultiOmicsModel:
    """Build + train a MultiOmicsModel on the given modality set; return it.

    ``pos_weight`` (default False, so within-cohort v0.2 behavior is unchanged):
    when True, weight the positive class in BCE by n_neg/n_pos. This matters for
    the cross-cohort transfer, where the train cohort (TCGA HER2 ~13%) is imbalanced
    and an unweighted boundary transfers poorly across platforms.

    ``gated`` (default False): use the input-conditioned softmax gate over modalities
    (v0.5 (2)), so the model can down-weight an unhelpful modality instead of letting
    a plain concat fusion dilute a strong one.

    ``modality_dropout`` (default 0.0; v0.6 (1)): per epoch, with this probability,
    zero one randomly chosen modality's input. It regularizes the model (and the gate)
    against over-relying on a single modality — the v0.5 gate collapsed to CNV because
    it learned a train-cohort preference; dropout forces each modality to stay usable.
    """
    torch.manual_seed(seed)
    dims = {m: modality_arrays[m].shape[1] for m in modalities}
    model = MultiOmicsModel(dims, latent_dim=latent_dim, gated=gated)
    tr = {m: _tensor(modality_arrays[m][train_idx]) for m in modalities}
    y_tr = _tensor(y[train_idx])

    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    if pos_weight:
        n_pos = float(y_tr.sum())
        pw = torch.tensor([(len(y_tr) - n_pos) / max(n_pos, 1.0)], dtype=torch.float32)
        loss_fn = nn.BCEWithLogitsLoss(pos_weight=pw)
    else:
        loss_fn = nn.BCEWithLogitsLoss()
    model.train()
    md_gen = torch.Generator().manual_seed(seed + 1)
    n_mod = len(modalities)
    for _ in range(n_epochs):
        tr_ep = tr
        if modality_dropout > 0.0 and n_mod > 1 and torch.rand(1, generator=md_gen).item() < modality_dropout:
            drop = modalities[int(torch.randint(n_mod, (1,), generator=md_gen).item())]
            tr_ep = {m: (torch.zeros_like(v) if m == drop else v) for m, v in tr.items()}
        opt.zero_grad()
        loss_fn(model(tr_ep), y_tr).backward()
        opt.step()
    return model


def auroc_of(
    model: MultiOmicsModel,
    modality_arrays: dict[str, np.ndarray],
    y: np.ndarray,
    modalities: tuple[str, ...],
    val_idx: np.ndarray,
) -> float:
    va = {m: _tensor(modality_arrays[m][val_idx]) for m in modalities}
    model.eval()
    with torch.no_grad():
        proba = torch.sigmoid(model(va)).cpu().numpy()
    return float(roc_auc_score(y[val_idx], proba))


def train_and_auroc(
    modality_arrays: dict[str, np.ndarray],
    y: np.ndarray,
    modalities: tuple[str, ...],
    *,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    **fit_kw,
) -> float:
    """Train on the modality set and return held-out AUROC."""
    model = fit_model(modality_arrays, y, modalities, train_idx=train_idx, **fit_kw)
    return auroc_of(model, modality_arrays, y, modalities, val_idx)


def run_ablation(
    modality_arrays: dict[str, np.ndarray],
    y: np.ndarray,
    *,
    axis: str = "",
    val_fraction: float = 0.3,
    seed: int = 0,
    **fit_kw,
) -> AblationResult:
    """Single train/val split ablation: {rna,meth} vs {rna,meth,cnv}."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(y))
    n_val = max(2, int(round(len(y) * val_fraction)))
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    base = train_and_auroc(
        modality_arrays, y, BASELINE_SET, train_idx=train_idx, val_idx=val_idx, seed=seed, **fit_kw
    )
    full = train_and_auroc(
        modality_arrays, y, FULL_SET, train_idx=train_idx, val_idx=val_idx, seed=seed, **fit_kw
    )
    return AblationResult(axis=axis, auroc_baseline=base, auroc_full=full, delta=full - base)
