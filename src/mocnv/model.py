"""Three-modality (or two, for the ablation) fusion classifier.

Per-modality encoders project each input to a shared latent; the latents are
concatenated and fused, then a head emits one binary logit. The active modality
set is configurable, so the **modality ablation** is just two instances of the
same model — ``{rna, meth}`` vs ``{rna, meth, cnv}`` — with everything else fixed.

This reuses the DMOI encoder/fusion building blocks (MLP encoders, concat-then-MLP
fusion). It is deliberately the focused ablation backbone; the per-pole
disagreement signal is a v0.3 refinement, not part of the v0.2 modality claim.
"""
from __future__ import annotations

from collections.abc import Sequence

import torch
from torch import nn

from mocnv.encoder import build_mlp, make_encoder


class MultiOmicsModel(nn.Module):
    """Encode each active modality, concat-fuse, classify.

    Args:
        modality_dims: ordered dict {modality_name -> input_dim}. Its keys define
            which modalities are active — drop ``cnv`` for the ablation baseline.
        latent_dim:    shared per-modality latent width.
        fuse_hidden / fuse_out / head_hidden / dropout: fusion + head sizing
            (v0.6 DMOI defaults).
    """

    def __init__(
        self,
        modality_dims: dict[str, int],
        *,
        latent_dim: int = 128,
        fuse_hidden: Sequence[int] = (128,),
        fuse_out: int = 64,
        head_hidden: int = 32,
        dropout: float = 0.3,
        gated: bool = False,
    ) -> None:
        super().__init__()
        if not modality_dims:
            raise ValueError("modality_dims must contain at least one modality")
        self.modalities: tuple[str, ...] = tuple(modality_dims)
        self.latent_dim = latent_dim
        self.gated = gated

        self.encoders = nn.ModuleDict({
            m: make_encoder(m, d, latent_dim=latent_dim, dropout=dropout)
            for m, d in modality_dims.items()
        })

        fuse_in = len(self.modalities) * latent_dim
        # Optional input-conditioned softmax gate over modalities: the model learns,
        # per sample, how much to weight each modality, so an unhelpful modality can be
        # down-weighted toward zero (default off, so v0.2-v0.4 concat fusion is unchanged).
        # Only meaningful with >=2 modalities.
        self.gate_net = nn.Linear(fuse_in, len(self.modalities)) if gated and len(self.modalities) > 1 else None
        self.fuse = build_mlp(fuse_in, fuse_hidden, fuse_out, dropout)
        self.head = nn.Sequential(
            nn.LayerNorm(fuse_out),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(fuse_out, head_hidden),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(head_hidden, 1),
        )

    def _fuse_input(self, latents: list[torch.Tensor]) -> tuple[torch.Tensor, torch.Tensor | None]:
        """Concatenate modality latents, applying the softmax gate if enabled.

        Returns (fused_input, gates) where gates is (batch, n_modalities) or None.
        """
        z = torch.cat(latents, dim=-1)
        if self.gate_net is None:
            return z, None
        gates = torch.softmax(self.gate_net(z), dim=-1)          # (batch, M)
        stacked = torch.stack(latents, dim=1)                    # (batch, M, latent)
        gated = (stacked * gates.unsqueeze(-1)).reshape(z.shape[0], -1)
        return gated, gates

    def forward(self, inputs: dict[str, torch.Tensor]) -> torch.Tensor:
        """inputs: {modality -> (batch, in_dim)}. Returns (batch,) logits."""
        missing = [m for m in self.modalities if m not in inputs]
        if missing:
            raise KeyError(f"missing modality inputs: {missing}")
        latents = [self.encoders[m](inputs[m]) for m in self.modalities]
        z, _ = self._fuse_input(latents)
        return self.head(self.fuse(z)).squeeze(-1)

    def gate_weights(self, inputs: dict[str, torch.Tensor]) -> dict[str, float] | None:
        """Mean per-modality gate over the batch (None unless gated with >=2 modalities)."""
        if self.gate_net is None:
            return None
        self.eval()
        with torch.no_grad():
            latents = [self.encoders[m](inputs[m]) for m in self.modalities]
            _, gates = self._fuse_input(latents)
        mean = gates.mean(dim=0)
        return {m: float(mean[i]) for i, m in enumerate(self.modalities)}
