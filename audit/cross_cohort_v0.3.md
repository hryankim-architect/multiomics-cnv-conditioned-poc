# Cross-cohort CNV transfer (TCGA -> METABRIC) — v0.3

Train TCGA cohort_v4 (n=421, HER2=56; RNA+meth+GISTIC2 CNV, class-weighted), score METABRIC cohort_v4 (n=1399, HER2=224; RNA QN->TCGA + meth silenced + **SNP6** CNA, amplicon-masked + z-scored). Which modality transfers across platforms?

| Setting (cross-cohort, METABRIC) | AUROC |
|---|---|
| RNA only | 0.684 |
| CNV only (SNP6, amplicon) | 0.762 |
| RNA + meth(silenced) — baseline | 0.770 |
| + CNV (full) | 0.752 |
| **CNV delta (full − base)** | **-0.018** |

METABRIC CNV IG top-5: STARD3, PGAP3, MIEN1, CTTN, ERBB2

Honest reading: TCGA GISTIC2 and METABRIC SNP6 are different platforms; gene-level alignment + per-gene z-scoring is the cross-platform bridge, not a validated normalization. With the TCGA HER2 class (~13%) trained class-weighted (`pos_weight`), the RNA+meth baseline transfers cleanly — an *unweighted* baseline scored sub-chance (an inverted artifact, not a finding) and inflated the CNV delta to a misleading +0.4; the per-modality table above is reported precisely so the headline does not hinge on a fragile baseline.

Against the calibrated baseline the **CNV delta is ~null**: the within-cohort +0.125 (v0.2, same-platform GISTIC2) does **not** survive the GISTIC2->SNP6 jump. The CNV-only row shows whether the SNP6 amplicon carries standalone cross-platform signal at all; where it does, it is redundant with RNA cross-cohort, so the full model gains nothing incremental. Strikingly, CNV-only here transfers at least as well as RNA-only — copy-number amplification is a discrete, platform-robust event, whereas microarray->RNA-seq expression transfer is noisier even after quantile normalization; the null full-model delta is therefore RNA/CNV redundancy on the HER2 axis (amplification drives over-expression), not a CNV transfer failure. The CNV IG top-5 still concentrates on the ERBB2 17q12 amplicon (STARD3/PGAP3/MIEN1/ERBB2) — the branch attends to the biologically correct locus, but cross-platform that signal is not additive over RNA. A null cross-cohort CNV delta is the recorded limit of cross-platform CNV, reported not hidden — exactly what the scope doc anticipated.
