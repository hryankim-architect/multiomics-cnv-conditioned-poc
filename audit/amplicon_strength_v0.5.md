# Amplicon strength predicts CNV transfer — v0.5 (1)

Per amplicon gene: **within-cohort strength** (single-gene AUROC on TCGA) vs **cross-cohort transfer** (1-feature logistic trained on TCGA, scored on METABRIC). Makes the v0.4 axis-level result (CNV helps where the amplicon defines the axis) a continuous per-gene law.

**Pooled across both axes (40 gene-points): Spearman rho = +0.836.**

### HER2 (TCGA n=421, 56 pos -> METABRIC n=1399, 224 pos)

Spearman(strength, transfer) **rho = +0.877**

| gene | pole | within-cohort strength (AUROC) | cross-cohort transfer (AUROC) |
|---|---|---|---|
| GRB7 | HER2 | 0.950 | 0.738 |
| ERBB2 | HER2 | 0.950 | 0.741 |
| PGAP3 | HER2 | 0.950 | 0.742 |
| MIEN1 | HER2 | 0.950 | 0.741 |
| PNMT | HER2 | 0.948 | 0.738 |
| STARD3 | HER2 | 0.948 | 0.732 |
| TCAP | HER2 | 0.948 | 0.738 |
| MED1 | HER2 | 0.838 | 0.669 |
| CASC8 | prolif | 0.639 | 0.607 |
| POU5F1B | prolif | 0.637 | 0.606 |
| MYC | prolif | 0.636 | 0.602 |
| PVT1 | prolif | 0.627 | 0.596 |
| CCND1 | prolif | 0.559 | 0.487 |
| ORAOV1 | prolif | 0.558 | 0.488 |
| FGF19 | prolif | 0.550 | 0.492 |
| FADD | prolif | 0.545 | 0.498 |
| ANO1 | prolif | 0.542 | 0.494 |
| CTTN | prolif | 0.537 | 0.497 |
| FGF4 | prolif | 0.537 | 0.495 |
| FGF3 | prolif | 0.534 | 0.494 |

### LumB (TCGA n=403, 127 pos -> METABRIC n=1175, 475 pos)

Spearman(strength, transfer) **rho = +0.388**

| gene | pole | within-cohort strength (AUROC) | cross-cohort transfer (AUROC) |
|---|---|---|---|
| PVT1 | prolif | 0.735 | 0.658 |
| CASC8 | prolif | 0.732 | 0.659 |
| POU5F1B | prolif | 0.731 | 0.658 |
| MYC | prolif | 0.728 | 0.657 |
| CCND1 | prolif | 0.594 | 0.637 |
| ORAOV1 | prolif | 0.593 | 0.635 |
| MED1 | HER2 | 0.591 | 0.521 |
| GRB7 | HER2 | 0.587 | 0.528 |
| ERBB2 | HER2 | 0.587 | 0.532 |
| PGAP3 | HER2 | 0.587 | 0.531 |
| MIEN1 | HER2 | 0.587 | 0.530 |
| PNMT | HER2 | 0.586 | 0.531 |
| TCAP | HER2 | 0.586 | 0.531 |
| STARD3 | HER2 | 0.586 | 0.532 |
| FGF19 | prolif | 0.574 | 0.632 |
| FADD | prolif | 0.570 | 0.631 |
| CTTN | prolif | 0.567 | 0.616 |
| ANO1 | prolif | 0.566 | 0.630 |
| FGF4 | prolif | 0.566 | 0.627 |
| FGF3 | prolif | 0.564 | 0.628 |

Honest reading: stronger within-cohort amplicon genes transfer better across platforms — the relationship is tight where the amplicon is focal and high-amplitude (HER2 / ERBB2 17q12) and looser where it is diffuse (LumB / proliferation). This also resolves the v0.4 IG puzzle: the ERBB2/17q12 block is individually the strongest for HER2, but it is near-perfectly co-amplified (collinear), so a multi-gene model splits weight across it and each member's *marginal* attribution is diluted — the v0.4 IG spread to proliferation was a collinearity artifact, not biology. The univariate view here is collinearity-free and puts ERBB2/17q12 on top. CNV's cross-cohort value is, gene for gene, a function of amplicon strength.
