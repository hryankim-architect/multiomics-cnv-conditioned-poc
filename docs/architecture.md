# Architecture — three-modality, CNV-conditioned

## Base (inherited from DMOI)

DMOI is a **hypothesis-conditioned** multi-omics classifier. For each pole (a
hypothesis, e.g. "this is HER2-driven"), each modality is projected through a
**pole mask** — the genes/loci that hypothesis says matter — encoded, and fused;
the model exposes a **disagreement signal** between poles as an interpretability
and calibration lever. RNA + DNA-methylation are the two original modalities.

```
modality x ── pole mask(pole p) ──▶ encoder_x ──▶ h(x, p)
                                                     │  for each pole p
                                fusion(h(rna,p), h(meth,p), …) ──▶ s(p)
                                disagreement(s(p1), s(p2)) ─┐
                                                            ▼
                                                    classifier head ──▶ call
```

## The CNV extension

A third branch is added with the **same shape contract** as the RNA/meth encoders:

- **`CNVEncoder`** — mirrors `RNAEncoder` / `MethEncoder` (the same small MLP
  family); input is gene-level copy number.
- **CNV pole masks** — the hypothesis-conditioning that makes CNV interpretable.
  Copy number is *amplicon-structured*, so a pole mask is a small set of loci:
  - HER2 pole → **ERBB2 amplicon (17q12)**
  - proliferation pole → **MYC (8q24)**, **CCND1 (11q13)**
  These are arguably *cleaner* for CNV than for RNA, since amplicon-driven subtypes
  are defined by copy number.
- **Three-perspective fusion** — the existing `PoleFuser` is extended from two
  inputs to three per pole; v0.6 hyperparameters are unchanged. Only modality
  count, masks, and the data layer differ.

## Data layer + harmonization

| Modality | TCGA-BRCA | METABRIC | Processing |
|---|---|---|---|
| RNA | HiSeqV2 | microarray | reused from dmoi |
| Methylation | HM450 | — | TCGA-only |
| **CNV (new)** | GISTIC2 gene-level (UCSC Xena) | SNP6 CNA (cBioPortal) | gene-level, see below |

**Harmonization decision (the real risk).** GISTIC2 emits discrete-ish gene-level
copy-number scores; METABRIC SNP6 CNA is a different platform/pipeline. We ingest
both at **gene level** and align to the shared RNA gene universe. Unlike RNA,
quantile normalization across platforms is **not obviously valid** for copy-number
calls — so the cross-cohort step (v0.3) is treated as a weaker, explicitly-caveated
check, not a headline. The decision and its limits are recorded in
`audit/cross_cohort_v0.3.md`.

## The deliverable — a modality ablation

For each task axis (e.g. HER2-vs-Luminal):

1. Train the model on **RNA + meth** and on **RNA + meth + CNV**, v0.6 config, 5-fold CV.
2. Report the **AUROC / balanced-accuracy delta** (mean ± std) — the modality effect.
3. **Per-pole Integrated Gradients** (Captum) on the CNV branch: does it key on the
   expected amplicon (ERBB2 for HER2)? Attribution that lands on the right locus is
   the evidence that any gain is real, not noise.
4. State the verdict honestly per axis — **CNV is expected to help amplicon-driven
   axes and add little elsewhere**; a null is a valid, on-brand result.

## What stays fixed

The model, encoders' hyperparameters, fusion, disagreement, and training loop are
the DMOI v0.6 design. This repo is a *modality-count + data-layer* extension, so the
ablation isolates the effect of the third modality rather than confounding it with
architecture changes.
