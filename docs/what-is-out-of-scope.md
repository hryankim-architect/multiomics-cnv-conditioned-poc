# What is out of scope

This file is **required** in every repo from the scaffold. The CI lint verifies it
exists; the PR template references it. This repo's value comes from being *small and
complete* — a single honest claim (does a third modality help?) with its limits on
record. If a PR proposes something on this list, it must answer one question:

> Why is this still out of scope?

## Default out-of-scope items

- **Statistical-power claims.** The demo uses small public subsets; effect sizes are
  illustrative, not conclusive.
- **Full-cohort reproduction.** Beyond the manifest cap requires editing both
  `data/manifest.yaml` and the README's subset claim.
- **Production hardening** (HA, RBAC, multi-tenant). The substrate provides the
  foundation; this repo does not re-implement it.
- **Cloud-cost optimization.** The demo runs on a single workstation.

## Per-project out-of-scope items

The capability is a *three-modality, CNV-conditioned* ablation. Anything that does
not serve that single demonstration is out of scope.

- **Beating dmoi's accuracy.** This repo does not chase a higher headline AUROC than
  `dmoi-brca-poc`. *The claim is "the architecture extends to a third modality, and
  here is honestly where CNV does and does not add signal" — not a new accuracy
  record.*
- **Strong cross-cohort CNV transfer.** *GISTIC2 (TCGA) and SNP6 (METABRIC) are
  different platforms; cross-platform copy-number harmonization is hard and the v0.3
  external check is expected to be weaker than dmoi's RNA transfer. That weakness is
  a reported result, not a bug to engineer away here.*
- **Segment-level / allele-specific CNV, focal-vs-arm calling, purity/ploidy
  correction.** *Gene-level copy number is the intentional representation; finer CNV
  modeling is a different capability.*
- **Additional modalities** (proteomics/RPPA, miRNA, ATAC, imaging). *Three
  modalities is the claim; more is a different repo.*
- **New task axes beyond those inherited from dmoi.** *The axes are reused so the
  ablation is comparable; inventing axes here would confound the modality effect.*
- **Wet-lab or causal validation of attributed amplicons.** *Integrated-Gradients
  attribution is interpretive evidence the model uses sensible biology, not a
  mechanistic claim.*
- **Real patient data beyond the public TCGA/METABRIC subsets.** *Governance; the
  unit tests use synthetic fixtures only.*

## How to add an item

Open a PR that places the item in the right section, states the reason in one
italicized sentence, and points to the originating issue. A slow list stays accurate.
