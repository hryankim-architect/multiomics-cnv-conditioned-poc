# Modality ablation (TCGA real cohort) — v0.2

Does adding the CNV branch (amplicon pole-masked) help over RNA+meth, per task axis? Same `MultiOmicsModel`, 5-fold StratifiedKFold (seed 42), v0.6 sizing. CNV is gene-level GISTIC2 restricted to the HER2 (ERBB2 17q12) and proliferation (MYC 8q24 / CCND1 11q13) amplicon masks, z-scored.

| Axis | n | rna+meth AUROC | rna+meth+cnv AUROC | CNV delta | CNV IG top-5 |
|---|---|---|---|---|---|
| HER2-vs-Luminal | 421 | 0.809 ± 0.118 | 0.934 ± 0.066 | +0.125 | PGAP3, ERBB2, GRB7, MIEN1, STARD3 |
| LumA-vs-LumB | 403 | 0.855 ± 0.030 | 0.848 ± 0.054 | -0.007 | PVT1, CASC8, MYC, POU5F1B, FGF19 |

Honest reading: CNV is expected to help the amplicon-driven HER2 axis (ERBB2 should lead the CNV attribution) and add little to LumA-vs-LumB. Cross-cohort CNV (METABRIC SNP6) is a recorded harder problem (v0.3); this is the TCGA within-cohort ablation. A null delta on an axis is a valid result.
