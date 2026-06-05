"""Model + ablation + attribution tests (torch-guarded; run on a machine with torch).

These pin the v0.2 headline on synthetic data: a 3-modality model trains, CNV
helps the amplicon-driven (HER2) axis more than the non-amplicon (Luminal) axis,
and the CNV-branch attribution keys on the HER2 amplicon. Thresholds are traced
from the planted synthetic structure; the local run is authoritative (tune if a
bound is brittle).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

torch = pytest.importorskip("torch")

from mocnv import attribution, eval_ablation, synth  # noqa: E402
from mocnv.model import MultiOmicsModel  # noqa: E402
from mocnv.priors import POLE_HER2_CNV  # noqa: E402


def _arrays(s: synth.SyntheticMultiomics) -> dict[str, np.ndarray]:
    return {"rna": s.rna, "meth": s.meth, "cnv": s.cnv}


def _inputs(s: synth.SyntheticMultiomics) -> dict[str, torch.Tensor]:
    return {m: torch.tensor(v, dtype=torch.float32) for m, v in _arrays(s).items()}


def test_forward_shape_three_modalities():
    s = synth.generate("HER2", seed=0, n=50)
    dims = {m: a.shape[1] for m, a in _arrays(s).items()}
    out = MultiOmicsModel(dims, latent_dim=16)(_inputs(s))
    assert out.shape == (50,)


def test_modality_subset_drops_cnv():
    s = synth.generate("HER2", seed=0, n=20)
    m = MultiOmicsModel({"rna": s.rna.shape[1], "meth": s.meth.shape[1]}, latent_dim=8)
    assert m.modalities == ("rna", "meth")
    out = m({"rna": torch.tensor(s.rna), "meth": torch.tensor(s.meth)})
    assert out.shape == (20,)


def test_cnv_helps_her2_more_than_luminal():
    s_h = synth.generate("HER2", seed=0, n=300)
    s_l = synth.generate("Luminal", seed=0, n=300)
    r_h = eval_ablation.run_ablation(_arrays(s_h), s_h.y, axis="HER2", latent_dim=16, n_epochs=80, seed=0)
    r_l = eval_ablation.run_ablation(_arrays(s_l), s_l.y, axis="Luminal", latent_dim=16, n_epochs=80, seed=0)
    # CNV contributes more on the amplicon-driven axis than the non-amplicon one.
    assert r_h.delta > r_l.delta
    assert r_h.auroc_full > 0.75       # CNV makes the HER2 axis work
    assert r_l.auroc_baseline > 0.70   # RNA carries the Luminal axis without CNV


def test_ig_keys_on_her2_amplicon():
    s = synth.generate("HER2", seed=0, n=300)
    arrays = _arrays(s)
    model = eval_ablation.fit_model(
        arrays, s.y, ("rna", "meth", "cnv"), train_idx=np.arange(300), latent_dim=16, n_epochs=120, seed=0
    )
    inputs = {m: torch.tensor(arrays[m], dtype=torch.float32) for m in ("rna", "meth", "cnv")}
    attr = attribution.integrated_gradients(model, inputs, "cnv", steps=16)
    top5 = {g for g, _ in attribution.rank_genes_by_attribution(attr, s.cnv_genes)[:5]}
    assert top5 & set(POLE_HER2_CNV)   # a HER2 amplicon gene is among the top CNV attributions


def test_pos_weight_and_single_modality_sets_train():
    # v0.3/v0.4 path: class-weighted training + the RNA-only / CNV-only modality sets.
    s = synth.generate("HER2", seed=0, n=200)
    arrays = _arrays(s)
    idx = np.arange(200)
    for mods in (eval_ablation.RNA_ONLY, eval_ablation.CNV_ONLY, eval_ablation.FULL_SET):
        model = eval_ablation.fit_model(
            arrays, s.y, mods, train_idx=idx, latent_dim=16, n_epochs=60, seed=0, pos_weight=True
        )
        auroc = eval_ablation.auroc_of(model, arrays, s.y, mods, idx)
        assert 0.0 <= auroc <= 1.0     # single-modality subsets are valid models
    full = eval_ablation.fit_model(
        arrays, s.y, eval_ablation.FULL_SET, train_idx=idx, latent_dim=16, n_epochs=120, seed=0, pos_weight=True
    )
    # class-weighted full model still recovers the planted HER2 structure
    assert eval_ablation.auroc_of(full, arrays, s.y, eval_ablation.FULL_SET, idx) > 0.75
