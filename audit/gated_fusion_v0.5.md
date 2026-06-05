# Gated fusion vs the LumB dilution — v0.5 (2)

v0.4 (plain concat) diluted a strong RNA signal when the weak CNV branch was added on LumB (delta -0.199). An input-conditioned softmax gate over modalities lets the model down-weight CNV. Cross-cohort (TCGA -> METABRIC, meth-free, class-weighted); `CNV gate` is the mean gate the gated model puts on the CNV modality at METABRIC.

| Axis | RNA-only | +CNV concat | concat delta | +CNV gated | gated delta | CNV gate |
|---|---|---|---|---|---|---|
| HER2 | 0.684 | 0.786 | +0.101 | 0.752 | +0.068 | 0.994 |
| LumB | 0.922 | 0.723 | -0.199 | 0.673 | -0.249 | 0.991 |

Honest reading (negative result): gating did **not** fix the LumB dilution — it made it worse (-0.249 vs concat -0.199) and also cost HER2 a little (+0.068 vs +0.101). The reason is in the gate column: the gate **collapsed to CNV on both axes** (CNV gate 0.99), so the gated model is effectively CNV-only (gated HER2 0.752 ≈ CNV-only 0.762; gated LumB 0.673 ≈ CNV-only 0.686). The gate is trained on TCGA, where the 20-gene amplicon CNV is low-dimensional and clean while RNA is 20530-dim and prone to overfit on n~400 — so the gate learns to *trust CNV* there. But that modality preference does not transfer: on LumB, RNA is the cross-platform-robust modality, yet the gate has already committed to CNV. So plain concat (which keeps both modalities) beats this gate on both axes, and the LumB dilution is **not** a simple concat artifact — even an adaptive gate fails because *which modality to trust does not itself transfer cross-cohort*. A regularized or cross-cohort-aware gate is future work; this naive gate is reported as a null/negative, not hidden.
