"""Per-modality MLP encoders (RNA / meth / CNV).

Mirrors the DMOI encoder design: LayerNorm + ReLU + Dropout between hidden layers,
a bare Linear final layer projecting to a shared latent dim. A single generic
``ModalityEncoder`` covers all three modalities; ``make_encoder`` supplies sensible
per-modality default hidden widths (CNV is the narrowest — a small, amplicon-masked
gene set rather than the full transcriptome).
"""
from __future__ import annotations

from collections.abc import Sequence

import torch
from torch import nn


def build_mlp(in_dim: int, hidden_dims: Sequence[int], out_dim: int, dropout: float) -> nn.Sequential:
    """MLP with LayerNorm+ReLU+Dropout between hidden layers; bare Linear final layer."""
    dims = [in_dim, *list(hidden_dims), out_dim]
    layers: list[nn.Module] = []
    for i in range(len(dims) - 1):
        layers.append(nn.Linear(dims[i], dims[i + 1]))
        if i < len(dims) - 2:
            layers.append(nn.LayerNorm(dims[i + 1]))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.Dropout(dropout))
    return nn.Sequential(*layers)


class ModalityEncoder(nn.Module):
    """Generic MLP encoder: (batch, in_dim) -> (batch, out_dim) latent."""

    def __init__(
        self,
        in_dim: int,
        hidden_dims: Sequence[int] = (256,),
        out_dim: int = 128,
        dropout: float = 0.3,
    ) -> None:
        super().__init__()
        self.in_dim = in_dim
        self.out_dim = out_dim
        self.net = build_mlp(in_dim, hidden_dims, out_dim, dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


# Per-modality default hidden widths (scaled to the modality's input size).
_DEFAULT_HIDDEN = {
    "rna": (1024, 256),
    "meth": (512,),
    "cnv": (256,),
}


def make_encoder(modality: str, in_dim: int, *, latent_dim: int = 128, dropout: float = 0.3) -> ModalityEncoder:
    """Build a ModalityEncoder with modality-appropriate default hidden widths.

    Falls back to a single 256-wide hidden layer for unknown modalities.
    """
    hidden = _DEFAULT_HIDDEN.get(modality, (256,))
    return ModalityEncoder(in_dim=in_dim, hidden_dims=hidden, out_dim=latent_dim, dropout=dropout)


def count_parameters(module: nn.Module) -> int:
    return sum(p.numel() for p in module.parameters() if p.requires_grad)
