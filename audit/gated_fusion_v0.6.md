# Modality-dropout vs the gate collapse — v0.6 (1)

v0.5 (2) showed the softmax gate collapsed to CNV (~0.99) and worsened LumB. v0.6 (1) adds **modality-dropout** training (p=0.5: per epoch, with this probability one modality's input is zeroed) to keep each modality usable and the gate adaptive. Cross-cohort, deltas vs RNA-only; values in parentheses are the mean CNV gate.

| Axis | RNA-only | concat delta | gated delta (CNV gate) | gated+dropout delta (CNV gate) |
|---|---|---|---|---|
| HER2 | 0.684 | +0.101 | +0.068 (0.99) | +0.075 (0.99) |
| LumB | 0.922 | -0.199 | -0.249 (0.99) | -0.236 (0.99) |

Honest reading: if **gated+dropout** moves the CNV gate off ~0.99 and lifts LumB's delta toward >= 0 (the gate stops discarding the cross-platform-robust RNA) while HER2 keeps its gain, modality-dropout fixes the v0.5 collapse — the dilution was a regularization problem. If the gate still collapses and LumB stays negative, the failure is deeper than regularization (the gate has no held-out signal to learn cross-cohort modality trust), and plain concat remains the honest default. Either outcome is reported.
