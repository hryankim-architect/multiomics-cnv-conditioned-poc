"""Integrated Gradients for one modality branch (lightweight, no Captum dependency).

Attributes the model's output to the features of a chosen modality by integrating
the input gradient along a straight path from a zero baseline to the input
(Sundararajan et al. 2017). For the CNV branch on an amplicon-driven axis, the
top-attributed genes should be the expected amplicon (ERBB2 for HER2) — the
evidence that any CNV gain is real, not noise.
"""
from __future__ import annotations

import numpy as np
import torch


def integrated_gradients(
    model: torch.nn.Module,
    inputs: dict[str, torch.Tensor],
    target_modality: str,
    *,
    baseline: torch.Tensor | None = None,
    steps: int = 32,
) -> np.ndarray:
    """Per-feature IG attribution for ``target_modality``.

    Returns a (n_samples, n_features) array of attributions. Other modalities are
    held at their input values; only the target modality is integrated.
    """
    model.eval()
    x = inputs[target_modality].detach()
    base = torch.zeros_like(x) if baseline is None else baseline.detach()
    grad_sum = torch.zeros_like(x)

    for step in range(1, steps + 1):
        alpha = step / steps
        interp = (base + alpha * (x - base)).clone().requires_grad_(True)
        fed = {m: (interp if m == target_modality else inputs[m].detach()) for m in inputs}
        out = model(fed).sum()
        (grad,) = torch.autograd.grad(out, interp)
        grad_sum = grad_sum + grad.detach()

    avg_grad = grad_sum / steps
    attr = (x - base) * avg_grad
    return attr.cpu().numpy()


def rank_genes_by_attribution(attr: np.ndarray, gene_names: list[str]) -> list[tuple[str, float]]:
    """Rank genes by mean absolute attribution across samples (descending)."""
    mean_abs = np.abs(attr).mean(axis=0)
    order = np.argsort(mean_abs)[::-1]
    return [(gene_names[i], float(mean_abs[i])) for i in order]
