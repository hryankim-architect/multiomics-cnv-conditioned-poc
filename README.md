# `multiomics-cnv-conditioned-poc`

> **One principle, applied here.** Pick the smallest representation that could carry the signal; measure it against an honest baseline; report the verdict faithfully — whether the added modality wins, ties, or loses. *Knowing whether a capability is real, rather than a flattering benchmark, is the point.*
>
> In this repo: **representation** a third, copy-number (CNV) pole-conditioned branch added to the DMOI RNA + methylation architecture → **baseline** the same model on RNA + methylation only → **verdict** a per-axis **modality ablation** — CNV is expected to help *amplicon-driven* axes (HER2 via the ERBB2 amplicon) and add little elsewhere. A null or axis-specific result is the honest deliverable.

This repo extends the [`dmoi-brca-poc`](https://github.com/hryankim-architect/dmoi-brca-poc) hypothesis-conditioned multi-omics architecture from **two modalities to three**, adding copy-number variation (CNV) as a pole-relevant modality alongside RNA-seq and DNA-methylation. It exists to answer one question honestly: **does a third modality add signal, where, and does its attribution key on the expected biology?**

It is a spin-off, not a `dmoi-brca-poc` point release — a third modality is a larger claim than "the architecture is reusable," with its own data layer and honest cross-cohort limits (see below), so it gets its own clean-room repo.

## The capability, in one diagram

```
 RNA-seq ──▶ RNA encoder ───┐
 HM450   ──▶ meth encoder ──┼─▶ pole-conditioned fusion ─▶ disagreement ─▶ head ─▶ call
 CNV     ──▶ CNV encoder ───┘        (per-pole gene/locus masks)
   (NEW: GISTIC2 gene-level copy number, pole-masked to amplicon loci)
```

Each modality is projected through **pole-specific masks** (the genes/loci a hypothesis says matter). For a HER2 axis the CNV pole mask centers on the **ERBB2 amplicon**; for proliferation, **MYC / CCND1**. The model fuses the three perspectives and exposes a disagreement signal, exactly as in DMOI — only the modality count, the masks, and the data layer change. v0.6 hyperparameters are kept.

## Headline deliverable — a modality ablation (not an accuracy race)

The result is **RNA+meth vs RNA+meth+CNV**, per task axis, with:
- AUROC / balanced-accuracy delta per axis (5-fold CV band),
- per-pole Integrated-Gradients attribution showing whether the CNV branch keys on the **expected amplicons**,
- an honest verdict — including "CNV adds nothing here" where that is what the data shows.

## Honest scope (stated up front)

- **Cross-cohort CNV is harder than RNA.** TCGA GISTIC2 and METABRIC SNP6 CNA are different platforms/pipelines; the quantile-normalization trick that made RNA transfer is not obviously valid for discrete-ish copy-number calls. The cross-cohort validation is therefore **weaker** than dmoi's RNA result — reported, not hidden.
- **No clean matched external at identical processing** for CNV — a recorded limitation.
- **Modality addition ≠ accuracy win.** The honest hypothesis is that CNV helps amplicon-driven axes and adds little for others. The ablation is the deliverable; a modest or null delta is acceptable (the dmoi v0.13 posture).
- The unit tests run on **small synthetic fixtures**; real cohorts are downloaded by the user (`scripts/download_*.sh`). This repo demonstrates the method and the engineering, not a benchmark claim about real cohorts.

## Substrate

Each run appends to a hash-chained NDJSON audit trail; MLflow logging is a no-op unless a server is configured; a deterministic canary smoke runs in under a second. ruff + pytest + English-only checks gate every change.

See [`docs/architecture.md`](docs/architecture.md) for the design and [`docs/what-is-out-of-scope.md`](docs/what-is-out-of-scope.md) for the boundaries. The build sequence is in [`ROADMAP.md`](ROADMAP.md).
