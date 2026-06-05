# Modality ablation (synthetic) — v0.2

Planted structure: the HER2 axis is CNV-informative (ERBB2 amplicon); the Luminal axis is CNV-flat (RNA carries the signal). The ablation trains the same model with and without the CNV branch on identical splits.

| Axis | rna+meth AUROC | rna+meth+cnv AUROC | CNV delta |
|---|---|---|---|
| HER2 | 0.643 | 0.998 | +0.355 |
| Luminal | 0.935 | 0.960 | +0.024 |

Expected: a positive CNV delta on HER2 and ~0 on Luminal. The real-cohort ablation (TCGA GISTIC2 + dmoi RNA/meth, 5-fold) lands next, with per-pole Integrated Gradients confirming the CNV branch keys on the ERBB2 amplicon.
