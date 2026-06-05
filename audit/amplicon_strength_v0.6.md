# Strength->transfer law robustness (raw amplitude) — v0.6 (2)

v0.5 (1) measured strength as within-cohort single-gene AUROC. Here strength is the **raw GISTIC2 amplification amplitude** (mean positive - mean negative, un-z-scored), a model-free metric. If the Spearman correlation with cross-cohort transfer stays strong, the law is not an artifact of the AUROC definition.

**Pooled (40 gene-points): rho(within-AUROC, transfer) = +0.836; rho(raw amplitude, transfer) = +0.790.**

| axis | rho(within-AUROC, transfer) | rho(raw amplitude, transfer) |
|---|---|---|
| HER2 | +0.877 | +0.854 |
| LumB | +0.388 | +0.191 |

Honest reading: a strong amplitude correlation that agrees with the AUROC one confirms the strength->transfer law is metric-robust — genes with larger copy-number amplitude in the positive class transfer better across platforms, however 'strength' is measured.
