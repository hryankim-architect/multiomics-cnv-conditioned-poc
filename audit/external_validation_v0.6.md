# External validation of the strength->transfer law (third cohort) — v0.6 (3)

Train per-gene single-gene classifiers on TCGA (HER2 axis, n=421, HER2=56); test transfer on the **MBC Project** (`brca_mbcproject_wagle_2017`, n=218, HER2+=76; genome-wide gene-level CNA, HER2 from `HER2_STATUS_REPORTED`). Does the v0.5 (1) strength->transfer law hold on an independent third cohort?

**Spearman(within-cohort strength, transfer to MBCProject) = +0.889** (METABRIC reference, HER2 axis: +0.877).

| gene | pole | within-cohort strength (TCGA AUROC) | transfer (MBCProject AUROC) |
|---|---|---|---|
| GRB7 | HER2 | 0.950 | 0.775 |
| ERBB2 | HER2 | 0.950 | 0.777 |
| PGAP3 | HER2 | 0.950 | 0.773 |
| MIEN1 | HER2 | 0.950 | 0.783 |
| PNMT | HER2 | 0.948 | 0.772 |
| STARD3 | HER2 | 0.948 | 0.776 |
| TCAP | HER2 | 0.948 | 0.772 |
| MED1 | HER2 | 0.838 | 0.727 |
| CASC8 | prolif | 0.639 | 0.662 |
| POU5F1B | prolif | 0.637 | 0.662 |
| MYC | prolif | 0.636 | 0.660 |
| PVT1 | prolif | 0.627 | 0.636 |
| CCND1 | prolif | 0.559 | 0.505 |
| ORAOV1 | prolif | 0.558 | 0.498 |
| FGF19 | prolif | 0.550 | 0.507 |
| FADD | prolif | 0.545 | 0.507 |
| ANO1 | prolif | 0.542 | 0.507 |
| CTTN | prolif | 0.537 | 0.505 |
| FGF4 | prolif | 0.537 | 0.503 |
| FGF3 | prolif | 0.534 | 0.515 |

Honest caveats: the MBC Project is **metastatic** (not primary) and HER2 is **patient-reported** (`HER2_STATUS_REPORTED`, not centrally assayed) — a noisier test than METABRIC. A positive rho here is independent-cohort evidence that the strength->transfer law generalizes; a weaker rho than METABRIC's +0.877 is consistent with the metastatic + self-reported-label noise, not necessarily a failure of the law. Reported, not hidden.
