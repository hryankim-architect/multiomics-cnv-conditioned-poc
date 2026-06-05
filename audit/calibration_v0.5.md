# Cross-cohort calibration — v0.5 (3)

Are the cross-cohort probabilities trustworthy, not just well-ranked? TCGA-trained (meth-free, class-weighted) models scored on METABRIC: **Brier** (lower better) and **ECE** (expected calibration error, lower better) next to AUROC, per modality set.

### HER2

| setting | AUROC | Brier | ECE |
|---|---|---|---|
| RNA-only | 0.684 | 0.241 | 0.327 |
| CNV-only | 0.762 | 0.177 | 0.166 |
| RNA+CNV | 0.786 | 0.176 | 0.152 |

### LumB

| setting | AUROC | Brier | ECE |
|---|---|---|---|
| RNA-only | 0.922 | 0.249 | 0.093 |
| CNV-only | 0.686 | 0.249 | 0.164 |
| RNA+CNV | 0.723 | 0.246 | 0.166 |

Honest reading: calibration tracks the v0.4 axis-specific value story. On **HER2**, where CNV helps the ranking, it also **improves calibration**: RNA-only is badly miscalibrated cross-platform (ECE 0.327 despite AUROC 0.684 — its probabilities are not trustworthy), and adding CNV cuts ECE to 0.152 (and Brier 0.241 -> 0.176). So CNV's HER2 value is, if anything, clearer in calibration than in AUROC. On **LumB**, where CNV hurts the ranking, it also **worsens calibration**: RNA-only is the best-calibrated (ECE 0.093) and adding CNV raises ECE to 0.166. Conclusion: a modality that genuinely carries the axis (CNV on HER2, RNA on LumB) is both better-ranking and better-calibrated cross-platform; adding the wrong modality degrades both. (No post-hoc recalibration is applied; these are the raw cross-cohort probabilities.)
