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

**Goal**: the METABRIC external check, with its weaker-than-RNA caveat reported, not hidden.

- [x] `cnv.load_cbioportal_cna` (SNP6 `data_CNA.txt` loader) + tests; `scripts/eval_cross_cohort.py` (TCGA-train → METABRIC-score: RNA QN→TCGA + meth silenced + SNP6 CNV amplicon-masked); download source fixed (git-LFS media + fallbacks)
- [ ] **Run** (Mac, after `scripts/download_metabric_cna.sh`): produces `audit/cross_cohort_v0.3.md` + cross-cohort CNV delta + METABRIC CNV IG (expected weaker than the TCGA +0.125 — honest cross-platform limit)
- [ ] README climax: cross-cohort row
- [ ] v0.3 tag + release notes

---

## Why this sequence

CNV data harmonization (v0.1) lands before the model (v0.2) because the
cross-platform decision (GISTIC2 vs SNP6) is the real risk, and the masks must be
locked before the branch is wired. Cross-cohort (v0.3) is last because it depends
on both the model and the harmonization, and is where the honest limit lives.
