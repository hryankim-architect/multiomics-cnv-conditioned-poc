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

- **Cross-cohort CNV is a platform jump.** TCGA GISTIC2 and METABRIC SNP6 CNA are different copy-number platforms; gene-level z-scoring is the bridge, not a validated normalization. Measured in v0.3 (below): standalone, CNV transfers *comparably to or better than* RNA across platforms, but it adds **nothing incremental** over RNA+meth — the two modalities are redundant on the HER2 axis. Reported, not hidden.
- **No clean matched external at identical processing** for CNV — a recorded limitation.
- **Modality addition ≠ accuracy win.** The honest hypothesis is that CNV helps amplicon-driven axes and adds little for others. The ablation is the deliverable; a modest or null delta is acceptable (the dmoi v0.13 posture).
- The unit tests run on **small synthetic fixtures**; real cohorts are downloaded by the user (`scripts/download_*.sh`). This repo demonstrates the method and the engineering, not a benchmark claim about real cohorts.

## Demo results — TCGA modality ablation (v0.2)

Same `MultiOmicsModel`, 5-fold StratifiedKFold (seed 42), v0.6 sizing. CNV is gene-level GISTIC2 restricted to the amplicon pole masks, z-scored. Reproduce: `python scripts/eval_ablation_real.py` (needs `dmoi-brca-poc` alongside for RNA/meth).

| Axis | n | rna+meth AUROC | rna+meth+cnv AUROC | CNV Δ | CNV IG top-5 |
|---|---|---|---|---|---|
| HER2-vs-Luminal | 421 | 0.809 ± 0.118 | **0.934 ± 0.066** | **+0.125** | PGAP3, ERBB2, GRB7, MIEN1, STARD3 |
| LumA-vs-LumB | 403 | 0.855 ± 0.030 | 0.848 ± 0.054 | −0.007 | PVT1, CASC8, MYC, POU5F1B, FGF19 |

**Reading.** CNV adds real signal on the amplicon-driven **HER2** axis (+0.125 AUROC), and the CNV attribution lands squarely on the **17q12 ERBB2 amplicon** — every one of the top-5 attributed genes is a co-amplified 17q12 gene, exactly the biology HER2 is defined by. On **LumA-vs-LumB**, which is not amplicon-defined, CNV adds nothing (−0.007, within noise) — the honest null. The third modality helps where the biology says it should and not elsewhere; that per-axis honesty is the deliverable, not a uniform accuracy bump. (Synthetic-fixture demo: `scripts/run_ablation_synth.py`.)

## Cross-cohort transfer — TCGA → METABRIC (v0.3)

Does the CNV branch survive a **platform jump**? Train on TCGA (RNA-seq + HM450 + GISTIC2 CNV, HER2 class-weighted), score METABRIC (RNA microarray quantile-normalized to TCGA, methylation silenced, **SNP6** CNA — a *different* copy-number platform — amplicon-masked + z-scored), HER2-vs-Luminal. Reproduce: `python scripts/eval_cross_cohort.py`.

| Setting (cross-cohort, METABRIC n=1399) | AUROC |
|---|---|
| RNA only | 0.684 |
| **CNV only** (SNP6, amplicon) | **0.762** |
| RNA + meth(silenced) — baseline | 0.770 |
| + CNV (full) | 0.752 |
| **CNV delta (full − base)** | **−0.018** |

**Reading.** The per-modality table separates two honest findings a single delta would hide. **(1)** Standalone, the SNP6 amplicon **CNV transfers across platforms at least as well as RNA** (0.762 vs 0.684) — copy-number amplification is a discrete, platform-robust event (GISTIC2 or SNP6, an amplification is an amplification), whereas microarray→RNA-seq expression transfer stays noisier even after quantile normalization. **(2)** In the full model the **CNV delta is null** (−0.018): CNV and RNA are *redundant* on the HER2 axis (ERBB2 amplification drives ERBB2 over-expression), so once RNA+meth reaches 0.770 the amplicon CNV adds nothing incremental. The within-cohort **+0.125** (v0.2, same-platform) does not reappear cross-cohort — **not because CNV fails to transfer, but because the baseline already carries the shared signal.** The CNV attribution still keys on the 17q12 ERBB2 amplicon (STARD3, PGAP3, MIEN1, ERBB2) cross-platform. Full audit, including the sub-chance-baseline trap we caught and fixed: [`audit/cross_cohort_v0.3.md`](audit/cross_cohort_v0.3.md).

## Second amplicon axis — amplicon-general or HER2-specific? (v0.4)

Does the v0.3 result generalize to another amplicon? v0.4 reruns the per-modality cross-cohort transfer on a second axis — **LumA-vs-LumB** (proliferation; the MYC 8q24 / CCND1 11q13 amplicon pole), built from METABRIC's PAM50 split. It also **drops the silenced-meth branch** (METABRIC has no methylation), so the baseline is RNA-only and the delta is uncontaminated. Reproduce: `python scripts/build_metabric_cohort_v2.py && python scripts/eval_cross_cohort_v0.4.py`.

| Axis (cross-cohort, METABRIC) | RNA-only | CNV-only | RNA+CNV | CNV Δ (vs RNA) |
|---|---|---|---|---|
| HER2-vs-Luminal | 0.684 | **0.762** | 0.786 | **+0.101** |
| LumA-vs-LumB | **0.922** | 0.686 | 0.723 | **−0.199** |

**Reading.** Both amplicons transfer across platforms standalone (CNV-only 0.762 and 0.686, both > chance) — copy-number **amplicon transfer is general**, not an ERBB2 quirk. But **CNV's *value* is axis-specific**: it helps where the amplicon *defines* the axis (HER2, +0.101, CNV-only ≥ RNA) and *hurts* where RNA defines it (LumB, −0.199 — adding the weaker CNV dilutes a strong 0.922 RNA signal in a plain concat fusion). So v0.3's "CNV ≥ RNA" is HER2/definitional-specific; the underlying transfer is general.

**Attribution diagnostic** ([`audit/ig_within_vs_cross_v0.4.md`](audit/ig_within_vs_cross_v0.4.md)). CNV attribution is **platform-stable** (within-cohort TCGA ≈ cross-cohort METABRIC IG, both axes — we tested for a platform effect and there is none). A CNV-only HER2 classifier keys on the co-amplified **proliferation** loci (8q24/11q13), not ERBB2 — on both platforms — because, copy-number-only, proliferation amplification separates HER2 from a LumA-dominated Luminal group better than ERBB2 does. The v0.2 *full* model's CNV branch keyed on ERBB2 because RNA carried the rest, so which locus CNV attributes to is **model-composition-dependent**. The HER2 copy-number signature is broader than ERBB2 alone.

**Method note — per-modality beats a single delta, twice.** A base-vs-full delta misled in *both* cross-cohort versions: v0.3's RNA+meth baseline went sub-chance from class imbalance (fixed with `pos_weight`), and the meth-silenced baseline went sub-chance again on LumB (fixed by dropping meth). Each time the per-modality columns (RNA-only, CNV-only) stayed clean and carried the real finding. Report the transfer profile, not one delta.

## v0.5 — three follow-ups (strength law · fusion null · calibration)

Three probes of the v0.4 result (CNV's value is axis-specific), in order.

**(1) Amplicon strength predicts transfer — a continuous law.** Per amplicon gene: within-cohort discriminative strength (single-gene AUROC on TCGA) vs cross-cohort single-gene transfer (TCGA→METABRIC). **Pooled Spearman ρ = +0.836** over 40 gene-points (HER2 +0.88, LumB +0.39). The ERBB2/17q12 block is individually strongest (within 0.95 → transfer 0.74); diffuse proliferation genes are weaker. This also resolves the v0.4 IG puzzle: ERBB2/17q12 is the strongest block but is near-perfectly co-amplified (collinear), so multi-gene attribution dilutes each member — the univariate view is collinearity-free and puts ERBB2 on top. ([`audit/amplicon_strength_v0.5.md`](audit/amplicon_strength_v0.5.md); sklearn, no torch.)

**(2) Gated fusion does *not* rescue the dilution — a null.** An input-conditioned softmax gate over modalities collapsed to CNV on both axes (gate ≈ 0.99), making the gated model ≈ CNV-only and *worsening* LumB (−0.249 vs concat −0.199). The gate learns a TCGA modality preference (trust the clean 20-gene CNV over 20530-dim RNA) that does not transfer — so the dilution is not a simple concat artifact, and *which modality to trust* does not transfer either. Plain concat beats the gate on both axes. ([`audit/gated_fusion_v0.5.md`](audit/gated_fusion_v0.5.md).)

**(3) Calibration tracks the axis-specific value.** Cross-cohort Brier/ECE: where CNV helps (HER2) it also improves calibration (RNA-only ECE 0.327 → RNA+CNV 0.152); where it hurts (LumB) it worsens it (0.093 → 0.166). A modality that genuinely carries the axis is both better-ranking and better-calibrated; the wrong one degrades both. Notably RNA-only on HER2 ranks adequately (0.684) but is badly miscalibrated (ECE 0.327) — CNV fixes the probabilities, not just the ranking. ([`audit/calibration_v0.5.md`](audit/calibration_v0.5.md).)

## Substrate

Each run appends to a hash-chained NDJSON audit trail; MLflow logging is a no-op unless a server is configured; a deterministic canary smoke runs in under a second. ruff + pytest + English-only checks gate every change.

See [`docs/architecture.md`](docs/architecture.md) for the design and [`docs/what-is-out-of-scope.md`](docs/what-is-out-of-scope.md) for the boundaries. The build sequence is in [`ROADMAP.md`](ROADMAP.md).
