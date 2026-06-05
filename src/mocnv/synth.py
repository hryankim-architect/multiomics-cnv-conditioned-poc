"""Synthetic 3-modality fixtures with planted amplicon structure (for tests).

The fixtures encode the repo's central hypothesis so tests can assert it without
real data or a trained model:

  axis="HER2"     -> CNV is INFORMATIVE: positives carry an elevated ERBB2-amplicon
                     copy-number signal; RNA/meth carry only a weak signal. A
                     CNV-only readout separates the classes.
  axis="Luminal"  -> CNV is FLAT (no amplicon difference): the signal lives in RNA
                     (+ a little meth). A CNV-only readout is at chance; CNV adds
                     nothing over RNA+meth.

That contrast is exactly what the v0.2 modality ablation should recover on real
data — CNV helps amplicon-driven axes and not others.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mocnv.priors import POLE_HER2_CNV, POLE_PROLIFERATION_CNV

VALID_AXES = ("HER2", "Luminal")


@dataclass
class SyntheticMultiomics:
    rna: np.ndarray          # n x n_rna
    meth: np.ndarray         # n x n_meth
    cnv: np.ndarray          # n x n_cnv
    y: np.ndarray            # n, in {0, 1}
    rna_genes: list[str]
    meth_probes: list[str]
    cnv_genes: list[str]
    axis: str

    @property
    def n(self) -> int:
        return self.y.shape[0]


def generate(
    axis: str = "HER2",
    *,
    seed: int = 0,
    n: int = 200,
    positive_fraction: float = 0.4,
    cnv_effect: float = 1.5,
    rna_effect: float = 1.0,
    noise: float = 1.0,
) -> SyntheticMultiomics:
    """Generate a planted-amplicon 3-modality cohort for the given task axis."""
    if axis not in VALID_AXES:
        raise ValueError(f"axis must be one of {VALID_AXES}, got {axis!r}")
    rng = np.random.default_rng(seed)

    cnv_genes = [*POLE_HER2_CNV, *POLE_PROLIFERATION_CNV] + [f"FILLER{i}" for i in range(10)]
    rna_genes = [f"RNA{i}" for i in range(40)]
    meth_probes = [f"cg{i:06d}" for i in range(30)]

    y = (rng.random(n) < positive_fraction).astype(np.int64)
    pos = y == 1

    cnv = rng.normal(0.0, noise, size=(n, len(cnv_genes))).astype(np.float32)
    rna = rng.normal(0.0, noise, size=(n, len(rna_genes))).astype(np.float32)
    meth = rng.normal(0.0, noise, size=(n, len(meth_probes))).astype(np.float32)

    her2_cols = [cnv_genes.index(g) for g in POLE_HER2_CNV]

    if axis == "HER2":
        # CNV informative: ERBB2 amplicon elevated in positives; RNA only weakly.
        cnv[np.ix_(pos, her2_cols)] += cnv_effect
        rna[pos, :5] += 0.4 * rna_effect
    else:  # Luminal
        # CNV flat (no amplicon difference); signal is in RNA (+ a little meth).
        rna[pos, :10] += rna_effect
        meth[pos, :5] += 0.5 * rna_effect

    return SyntheticMultiomics(
        rna=rna, meth=meth, cnv=cnv, y=y,
        rna_genes=rna_genes, meth_probes=meth_probes, cnv_genes=cnv_genes, axis=axis,
    )
