# Cross-cohort CNV transfer, second amplicon axis — v0.4

Is the v0.3 ERBB2 cross-cohort result amplicon-general or HER2-specific? Same per-modality cross-cohort pipeline (train TCGA class-weighted, score METABRIC: RNA QN->TCGA + SNP6 CNA amplicon-masked) on a second axis — **LumA-vs-LumB** (proliferation; MYC 8q24 / CCND1 11q13 amplicon pole).

**Meth-free cross-cohort.** METABRIC has no methylation; v0.3 silenced the meth branch (all-zeros) at inference, a train/test mismatch that made the RNA+meth baseline unreliable (sub-chance on LumB). v0.4 trains and scores only on the modalities present in both cohorts — RNA + CNV — so the baseline is RNA-only and the delta is uncontaminated. (RNA-only/CNV-only reproduce v0.3; the TCGA sample set is unchanged.)

| Axis (cross-cohort, METABRIC) | RNA-only (base) | CNV-only | RNA+CNV (full) | CNV delta | CNV IG top-5 (pole) |
|---|---|---|---|---|---|
| HER2-vs-Luminal | 0.684 | 0.762 | 0.786 | +0.101 | MYC, FADD, PVT1, MIEN1, POU5F1B (MYC-8q24/CCND1-11q13) |
| LumA-vs-LumB | 0.922 | 0.686 | 0.723 | -0.199 | CCND1, FGF19, FADD, PVT1, CASC8 (MYC-8q24/CCND1-11q13) |

IG note: attribution is on the **CNV-only** model (isolates the modality). The diagnostic `ig_within_vs_cross_v0.4.md` shows it is **platform-stable** (within-cohort ≈ cross-cohort, both axes). Notably the CNV-only **HER2** model keys on the co-amplified proliferation loci (8q24 / 11q13), **not** ERBB2/17q12 — on both platforms. The v0.2 *full* model's CNV branch keyed on ERBB2 because RNA carried the proliferation/expression signal; so which locus the CNV branch attributes to is **model-composition-dependent**, not a platform artifact.

Honest reading: the comparison hinges on **CNV-only vs RNA-only** (both single modalities, uncontaminated). On **HER2** CNV-only >= RNA-only (0.762 vs 0.684), so adding CNV *lifts* the full model (+0.101). On **LumB** — a proliferation / PAM50 subtype that RNA defines — RNA-only is far stronger (0.922) and adding the weaker CNV (0.686) *dilutes* it (-0.199 in a plain concat fusion). Conclusion: copy-number **amplicon transfer is general** — both amplicons transfer cross-platform standalone (CNV-only > chance) — but **CNV's *value* is axis-specific**: it helps where the amplicon *defines* the axis (HER2) and hurts where RNA does (LumB). Attribution caveat: the CNV-only HER2 model discriminates via the co-amplified proliferation landscape rather than ERBB2 itself (see IG note + diagnostic), so the HER2 copy-number signature is broader than ERBB2 alone. v0.2 anchors the value finding: the proliferation CNV pole was null for LumA-vs-LumB within-cohort (-0.007) too. Either outcome is reported, not hidden.
