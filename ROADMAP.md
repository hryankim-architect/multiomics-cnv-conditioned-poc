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

## v0.1 — CNV data layer + pole masks (next)

**Goal**: load and harmonize the CNV modality; define the amplicon-centered pole masks.

- [ ] `cohort.py` — TCGA-BRCA + METABRIC cohort loaders (reuse dmoi cohort schema), aligned to the existing RNA/meth sample IDs
- [ ] `cnv.py` — GISTIC2 gene-level ingest + cross-platform harmonization decision (and its recorded limits), gene→locus alignment
- [ ] `priors.py` — CNV pole masks: HER2 → ERBB2 amplicon (17q12); proliferation → MYC (8q24) / CCND1 (11q13)
- [ ] `synth.py` — small synthetic 3-modality fixtures with planted amplicon structure (HER2 axis CNV-informative; Luminal axis CNV-flat) for tests
- [ ] Tests: harmonization round-trip, mask coverage, synthetic ablation sanity
- [ ] v0.1 tag

---

## v0.2 — three-modality model + ablation (the headline)

**Goal**: add the CNV branch, run the modality ablation, attribute it.

- [ ] `encoder.py` — `CNVEncoder` (mirrors RNA/Meth encoders); `model.py` — 3-perspective pole fusion (extend dmoi `PoleFuser`), v0.6 hyperparameters unchanged
- [ ] `eval_ablation.py` — RNA+meth vs RNA+meth+CNV, 5-fold CV AUROC/bacc per axis; per-pole IG (Captum) showing CNV keys on the expected amplicons
- [ ] `audit/ablation_v0.2.md` — the honest per-axis verdict (including any null)
- [ ] Tests: 3-modality forward shape, ablation determinism, IG-on-amplicon assertion (synthetic)
- [ ] v0.2 tag

---

## v0.3 — cross-cohort honest limit

**Goal**: the METABRIC external check, with its weaker-than-RNA caveat reported, not hidden.

- [ ] METABRIC SNP6 CNA external scoring (RNA+meth+CNV), cross-platform caveat documented
- [ ] `audit/cross_cohort_v0.3.md` — the honest cross-cohort CNV result (expected weaker than dmoi's RNA transfer)
- [ ] README climax: the modality-ablation table + cross-cohort row
- [ ] v0.3 tag + release notes summarizing the three-modality claim and its limits

---

## Why this sequence

CNV data harmonization (v0.1) lands before the model (v0.2) because the
cross-platform decision (GISTIC2 vs SNP6) is the real risk, and the masks must be
locked before the branch is wired. Cross-cohort (v0.3) is last because it depends
on both the model and the harmonization, and is where the honest limit lives.
