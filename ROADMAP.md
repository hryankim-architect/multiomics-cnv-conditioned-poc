# Roadmap — `multiomics-cnv-conditioned-poc`

Public build plan: a clean-room three-modality (RNA + meth + **CNV**) extension of
the DMOI architecture, delivering an honest per-axis **modality ablation**. Short
and dated so anyone clicking the repo sees exactly where the work stands.

Design rationale: `docs/architecture.md`. Spin-off decision + honest limits: the
private plan `~/Downloads/AI/cnv-multimodal-spinoff-plan.md`.

---

## v0.0 — scaffold + data layer (2026-06-04)

**Goal**: substrate-clean repo with the CNV datasets downloading, so the modeling
work can start from a known-good base.

- [x] Inherit the shared scaffold (audit / canary / tracking, ruff + pytest + English-only)
- [x] Package `mocnv`, `pyproject.toml`, canary smoke + audit-chain tests green
- [x] `README.md` (honest-scope + capability claim), `docs/architecture.md`, `docs/what-is-out-of-scope.md`, this ROADMAP
- [x] `data/manifest.yaml` + `data/README.md` (provenance + clean-room note)
- [x] Data-download scripts: `scripts/download_tcga_cnv.sh` (GISTIC2 / UCSC Xena), `scripts/download_metabric_cna.sh` (SNP6 CNA / cBioPortal)
- [ ] Push to GitHub; scaffold CI green
- [ ] v0.0 tag

---

## v0.1 — CNV data layer + pole masks (2026-06-05)

**Goal**: load and harmonize the CNV modality; define the amplicon-centered pole masks.

- [x] `cohort.py` — cohort-table builder with the 3-modality intersection (`has_rna`/`has_meth`/`has_cnv`); real TCGA/METABRIC clinical→group assignment reuses dmoi at v0.2
- [x] `cnv.py` — GISTIC2 gene-level ingest (real file: 1080 × 24,776), per-gene harmonization (the cross-platform decision, with its recorded limit), gene alignment, TCGA barcode normalization
- [x] `priors.py` — CNV pole masks: HER2 → ERBB2 amplicon (17q12, 8 genes); proliferation → MYC (8q24) + CCND1 (11q13, 12 genes). All present in the real gene universe (8/8, 12/12)
- [x] `synth.py` — synthetic 3-modality fixtures with planted amplicon structure (HER2 axis CNV-informative; Luminal axis CNV-flat) for tests
- [x] Tests: GISTIC2 load, alignment/harmonization, mask coverage, synthetic ablation sanity (20 tests green)
- [ ] v0.1 tag

---

## v0.2 — three-modality model + ablation (the headline)

**Goal**: add the CNV branch, run the modality ablation, attribute it.

- [x] `encoder.py` — `ModalityEncoder` / `make_encoder` (RNA/meth/CNV); `model.py` — `MultiOmicsModel`, per-modality encoders + concat-fuse + head, configurable modality subset for the ablation, v0.6 sizing
- [x] `eval_ablation.py` — RNA+meth vs RNA+meth+CNV on identical splits (`run_ablation`/`fit_model`); `attribution.py` — lightweight Integrated Gradients (no Captum dep) + gene ranking
- [x] Tests (torch-guarded): 3-modality forward shape, modality-subset, CNV-helps-HER2-more-than-Luminal, IG-keys-on-ERBB2-amplicon (synthetic)
- [x] `scripts/run_ablation_synth.py` — synthetic ablation demo + `audit/ablation_synth_v0.2.md`
- [x] **Real-cohort ablation** — `scripts/eval_ablation_real.py`: TCGA GISTIC2 + dmoi RNA/meth, 5-fold. **HER2-vs-Luminal CNV Δ +0.125 (0.809→0.934; ERBB2 amplicon leads the CNV IG); LumA-vs-LumB Δ −0.007 (honest null).** See `audit/ablation_v0.2.md`
- [x] v0.2 tag

---

## v0.3 — cross-cohort honest limit

**Goal**: the METABRIC external check — does CNV survive the GISTIC2→SNP6 platform jump? Reported honestly, including the sub-chance-baseline trap and the RNA/CNV redundancy finding.

- [x] `cnv.load_cbioportal_cna` (SNP6 `data_CNA.txt` loader) + tests; `scripts/eval_cross_cohort.py` (TCGA-train → METABRIC-score: RNA QN→TCGA + meth silenced + SNP6 CNV amplicon-masked, HER2 class-weighted via opt-in `pos_weight`); download source fixed (git-LFS media + fallbacks)
- [x] **Run** (Mac): `audit/cross_cohort_v0.3.md`. Per-modality cross-cohort AUROC — **RNA-only 0.684 | CNV-only 0.762 | RNA+meth 0.770 | +CNV 0.752 (delta −0.018)**. Honest finding: CNV transfers *better* than RNA standalone, but is **redundant** with RNA on HER2 (null full-model delta); ERBB2 amplicon leads CNV IG cross-platform. (An unweighted baseline scored sub-chance and inflated the delta to +0.4 — caught + corrected with `pos_weight`.)
- [x] README climax: cross-cohort per-modality table
- [x] v0.3 tag + release notes

---

## v0.4 — second amplicon axis: amplicon-general or HER2-specific?

**Goal**: test whether the v0.3 result generalizes to a second amplicon (MYC 8q24 / CCND1 11q13) on the LumA-vs-LumB axis, with a meth-free (uncontaminated) cross-cohort baseline.

- [x] `cohort.build_metabric_cohort_v2` + `scripts/build_metabric_cohort_v2.py` — split METABRIC's Luminal lump into LumA/LumB via PAM50 `CLAUDIN_SUBTYPE` (LumA 700 / LumB 475); tests
- [x] `scripts/eval_cross_cohort_v0.4.py` — two-axis per-modality cross-cohort, **meth-free** (RNA+CNV only; drops the silenced-meth branch that broke the LumB baseline). pos_weight + RNA-only/CNV-only/RNA+CNV
- [x] **Run** (Mac): `audit/cross_cohort_v0.4.md`. **HER2 RNA 0.684 | CNV 0.762 | RNA+CNV 0.786 (Δ +0.101)**; **LumB RNA 0.922 | CNV 0.686 | RNA+CNV 0.723 (Δ −0.199)**. Finding: amplicon transfer is **general** (both CNV-only > chance), but CNV's *value* is **axis-specific** — helps where the amplicon defines the axis (HER2), dilutes where RNA does (LumB)
- [x] **IG diagnostic** — `scripts/ig_within_vs_cross_v0.4.py` / `audit/ig_within_vs_cross_v0.4.md`: attribution is **platform-stable** (within ≈ cross); the CNV-only HER2 model keys on co-amplified proliferation loci (not ERBB2) — **model-composition-dependent**, not a platform effect (refuted that hypothesis)
- [x] README v0.4 section (two-axis table + diagnostic + "per-modality beats a single delta, twice")
- [x] v0.4 tag + release notes

---

## v0.5 — three follow-ups: strength law, fusion null, calibration

**Goal**: probe the v0.4 axis-specific result three ways — is amplicon strength a continuous predictor, is the LumB dilution a fixable fusion artifact, and does CNV's value extend to calibration?

- [x] **(1) Amplicon strength law** — `mocnv/strength.py` + `scripts/amplicon_strength_v0.5.py` (sklearn, no torch) + tests. Per-gene within-cohort strength vs cross-cohort transfer: **pooled Spearman ρ = +0.836** (HER2 +0.88, LumB +0.39). Resolves the v0.4 IG puzzle (ERBB2/17q12 strongest but collinear -> diluted marginal attribution). `audit/amplicon_strength_v0.5.md`
- [x] **(2) Gated fusion (null)** — `MultiOmicsModel(gated=True)` input-conditioned softmax gate + `scripts/eval_gated_fusion_v0.5.py` + test. Gate **collapsed to CNV** (~0.99) on both axes -> worsened LumB (−0.249); the TCGA modality preference does not transfer, plain concat wins. Dilution is not a naive-fusion artifact. `audit/gated_fusion_v0.5.md`
- [x] **(3) Cross-cohort calibration** — `mocnv/calibration.py` (Brier + ECE) + `scripts/eval_calibration_v0.5.py` + tests. Calibration tracks value: CNV improves HER2 (ECE 0.327 -> 0.152), worsens LumB (0.093 -> 0.166); RNA-only HER2 ranks ok but is badly miscalibrated. `audit/calibration_v0.5.md`
- [x] README v0.5 section (three follow-ups)
- [x] v0.5 tag + release notes

---

## v0.6 — robustness + external validation (and a second fusion null)

**Goal**: stress-test the v0.5 results — can the gate collapse be regularized away, is the strength law metric-robust, and does it hold on an independent third cohort?

- [x] **(1) Modality-dropout gate (null)** — `fit_model(modality_dropout=)` + `scripts/eval_gated_fusion_v0.6.py`. Dropout (p=0.5) left the CNV gate at ~0.99 on both axes; gated+dropout still below concat (LumB −0.236, HER2 +0.075). Collapse is not a regularization problem — plain concat stays the default. `audit/gated_fusion_v0.6.md`
- [x] **(2) Strength-law robustness** — `mocnv/strength.py::per_gene_amplitude` + `scripts/amplicon_strength_v0.6.py` + test. Raw GISTIC2 amplitude (model-free) gives pooled **ρ = +0.790** vs AUROC-based +0.836 — the law is metric-robust. `audit/amplicon_strength_v0.6.md`
- [x] **(3) Third-cohort external validation** — `scripts/download_thirdcohort_probe.sh` (probes 6 cBioPortal studies) found **MBC Project** (`brca_mbcproject_wagle_2017`, genome-wide CNA + reported HER2); `scripts/external_validation_v0.6.py`. Strength→transfer **ρ = +0.889** (n=218, HER2+=76), matching METABRIC's +0.877 — the law generalizes to an independent metastatic cohort. `audit/external_validation_v0.6.md`
- [x] README v0.6 section
- [ ] v0.6 tag + release notes

---

## Why this sequence

CNV data harmonization (v0.1) lands before the model (v0.2) because the
cross-platform decision (GISTIC2 vs SNP6) is the real risk, and the masks must be
locked before the branch is wired. Cross-cohort (v0.3) is last because it depends
on both the model and the harmonization, and is where the honest limit lives.
