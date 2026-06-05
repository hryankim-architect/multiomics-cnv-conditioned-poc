#!/usr/bin/env python3
"""Demonstrate the modality ablation on synthetic data (no real cohorts, no GPU).

Trains {rna, meth} vs {rna, meth, cnv} on a planted HER2 axis (CNV-informative)
and a Luminal axis (CNV-flat, RNA carries it), printing the held-out AUROC delta
and writing an audit doc. This is the headline shape the real-cohort ablation
(TCGA GISTIC2 + dmoi RNA/meth) should reproduce.

Reproduce:  python scripts/run_ablation_synth.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv import audit, eval_ablation, synth  # noqa: E402


def main() -> int:
    rows = []
    print(f"{'axis':10s}{'rna+meth':>12s}{'+cnv':>10s}{'delta':>10s}")
    for axis in ("HER2", "Luminal"):
        s = synth.generate(axis, seed=0, n=400)
        r = eval_ablation.run_ablation(
            {"rna": s.rna, "meth": s.meth, "cnv": s.cnv}, s.y,
            axis=axis, latent_dim=16, n_epochs=100, seed=0,
        )
        rows.append(r)
        print(f"{axis:10s}{r.auroc_baseline:>12.3f}{r.auroc_full:>10.3f}{r.delta:>+10.3f}")

    audit.emit(
        "ablation_synth_v0.2", "demo",
        {r.axis: {"baseline": r.auroc_baseline, "full": r.auroc_full, "delta": r.delta} for r in rows},
    )

    out = REPO / "audit"
    out.mkdir(exist_ok=True)
    md = out / "ablation_synth_v0.2.md"
    table = "\n".join(
        f"| {r.axis} | {r.auroc_baseline:.3f} | {r.auroc_full:.3f} | {r.delta:+.3f} |" for r in rows
    )
    md.write_text(
        "# Modality ablation (synthetic) — v0.2\n\n"
        "Planted structure: the HER2 axis is CNV-informative (ERBB2 amplicon); the "
        "Luminal axis is CNV-flat (RNA carries the signal). The ablation trains the "
        "same model with and without the CNV branch on identical splits.\n\n"
        "| Axis | rna+meth AUROC | rna+meth+cnv AUROC | CNV delta |\n|---|---|---|---|\n"
        + table
        + "\n\nExpected: a positive CNV delta on HER2 and ~0 on Luminal. The "
        "real-cohort ablation (TCGA GISTIC2 + dmoi RNA/meth, 5-fold) lands next, "
        "with per-pole Integrated Gradients confirming the CNV branch keys on the "
        "ERBB2 amplicon.\n"
    )
    print(f"\nwrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
