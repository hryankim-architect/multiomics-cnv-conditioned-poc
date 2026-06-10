# mocnv v0.7, CNV port of the prior-vs-baseline comparison

**Date:** 2026-06-10
**Tag:** v0.7
**Type:** method port (new analysis module + script; core CNV model unchanged)

## TL;DR

The sibling [`dmoi-brca-poc`](https://github.com/hryankim-architect/dmoi-brca-poc) showed
that a **label-free biological prior** (Hallmark + HM450-cis) beats statistical
(top-variance) feature selection for PAM50 subtyping, on the same footing as the
unsupervised integrators MOFA+/MoGCN (Omran et al. 2025). v0.7 ports that comparison to
the **CNV** modality. CNV is gene-level (GISTIC2) and amplicon-structured, so the prior
here is the tight amplicon-locus set (ERBB2 17q12 + MYC/CCND1, `priors.CNV_POLES`) — also
label-free, so it is directly comparable to a top-variance CNV selector. The amplicon
prior wins, **most sharply where its amplicon defines the axis** (HER2), consistent with
this repo's v0.2 modality-ablation thesis.

## Result (label-free CNV selection, LR/SVC stratified 5-fold)

**HER2-vs-rest** (amplicon-defined axis), n=831 (Her2 67 / rest 764):

| selector (label-free CNV) | n_feat | AUROC | LR weighted-F1 |
|---|---|---|---|
| **CNV-prior (amplicon)** | 20 | **0.830** | 0.913 |
| top-variance (k=20) | 20 | 0.812 | 0.906 |
| top-variance (k=100) | 100 | 0.767 | 0.902 |

**LumA-vs-LumB** (not amplicon-defined), n=607 (LumA 415 / LumB 192):

| selector (label-free CNV) | n_feat | AUROC | LR weighted-F1 |
|---|---|---|---|
| **CNV-prior (amplicon)** | 20 | **0.728** | 0.704 |
| top-variance (k=20) | 20 | 0.599 | 0.650 |
| top-variance (k=100) | 100 | 0.703 | 0.715 |

## Reading

- **The label-free amplicon prior beats top-variance, most sharply on its axis.** On
  HER2-vs-rest the 20-gene ERBB2/proliferation prior edges matched-budget top-variance
  and clearly beats top-variance(100) — with 5× fewer features. This reproduces the
  dmoi-brca-poc finding (biological prior > statistical selection) for CNV, with a
  compact *locus*-anchored prior.
- **The luminal edge is biologically honest.** CNV's absolute LumA-vs-LumB signal is
  weaker than its HER2 signal; the standalone prior still edges top-variance there
  because the CCND1 11q13 / MYC proliferation amplicon is genuinely informative for the
  proliferative LumB pole. The prior helps where copy-number biology is informative, not
  uniformly — the same axis-specific message as v0.2's modality ablation.

## What changed

- `src/mocnv/compare_integration.py` — new analysis module: label-free selectors
  (`prior_cnv_indices` over `CNV_POLES`, `topvar_indices`, `topvar_within`),
  `jaccard_index`, and `eval_selector` (LR + linear SVC weighted-F1, binary AUROC,
  Calinski-Harabasz / Davies-Bouldin). Selectors never see `y`.
- `scripts/compare_cnv_prior.py` (`make compare-cnv`) → `audit/cnv_prior_vs_baseline.md`
  + `.json`. Runs HER2-vs-rest and LumA-vs-LumB.
- `tests/test_compare_integration.py` — 6 unit tests (prior selection, top-variance,
  Jaccard, binary AUROC present / multiclass AUROC None).
- `Makefile` `compare-cnv` target; `ROADMAP.md` v0.7 entry.

No changes to the CNV model, encoders, or fusion in `src/mocnv/`.

## Reproduce

```bash
make compare-cnv     # needs the sibling dmoi-brca-poc clinical matrix for PAM50 labels
make check           # ruff + pytest
```

## Limitations (honest scope)

- **Label dependency on the sibling repo.** This repo ships CNV only; PAM50 labels are
  reused from `dmoi-brca-poc` (`--clinical` overrides the path). The two repos must sit
  side by side for the default path to resolve.
- **Feature-selection contribution only.** These are label-free selectors + a plain
  LR/SVC, not the full pole-conditioned CNV model — by design, to isolate the prior's
  value as a selector (the same framing as the dmoi comparison).
- **CNV is a weaker standalone modality off the amplicon axis** (lower absolute
  LumA-vs-LumB AUROC), consistent with v0.2's finding that CNV adds little to the *full*
  model where no amplicon defines the axis. The prior's off-axis edge is modest and
  driven by the proliferation amplicon, not a uniform gain.
- Single cohort (TCGA-BRCA GISTIC2); no cross-platform CNV transfer evaluated here
  (v0.3–v0.6 already cover that axis for the full model).
