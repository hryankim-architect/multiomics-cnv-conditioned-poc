"""CNV data layer: GISTIC2 gene-level ingest, harmonization, gene alignment.

GISTIC2 ``all_data_by_genes`` is a genes x samples table of continuous copy-number
scores (log2-ratio-like). We load it as samples x genes, collapse duplicate gene
symbols, optionally harmonize per gene (the cross-platform decision; see the honest
limit in ``docs/what-is-out-of-scope.md``), and align to a shared gene universe.

Cross-platform note: TCGA GISTIC2 and METABRIC SNP6 CNA are different platforms.
Per-gene z-scoring (``harmonize_gene_level``) puts each cohort on a comparable
internal scale, but it is **not** a validated cross-platform normalization — the
v0.3 cross-cohort check is explicitly weaker for CNV than for RNA.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class CNVMatrix:
    values: np.ndarray        # samples x genes, float32
    sample_ids: list[str]
    gene_names: list[str]

    @property
    def n_samples(self) -> int:
        return self.values.shape[0]

    @property
    def n_genes(self) -> int:
        return self.values.shape[1]


def load_gistic2(path: str | Path) -> CNVMatrix:
    """Load a GISTIC2 gene-level table (genes x samples) as samples x genes.

    Duplicate gene symbols are collapsed by mean (a harmonization choice recorded
    here, not hidden). Reads .gz transparently.
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if df.index.has_duplicates:
        df = df.groupby(level=0).mean()
    gene_names = [str(g) for g in df.index]
    sample_ids = [str(s) for s in df.columns]
    values = np.nan_to_num(df.to_numpy(dtype=np.float32).T, nan=0.0)  # missing -> neutral 0
    return CNVMatrix(values=values, sample_ids=sample_ids, gene_names=gene_names)


def load_cbioportal_cna(path: str | Path, *, sample_ids: set[str] | None = None) -> CNVMatrix:
    """Load a cBioPortal `data_CNA.txt` (METABRIC SNP6) as samples x genes.

    Format: `Hugo_Symbol <tab> Entrez_Gene_Id <tab> <sample columns…>`, identical
    to the METABRIC mRNA matrix. The Entrez column is dropped; duplicate Hugo
    symbols are collapsed by mean. Optionally restrict to ``sample_ids``.

    Note: METABRIC CNA is a different platform (SNP6) from TCGA GISTIC2 — aligning
    both at gene level is the cross-platform decision whose limits v0.3 reports.
    """
    df = pd.read_csv(path, sep="\t", low_memory=False)
    gene_col = df.columns[0]  # Hugo_Symbol
    sample_cols = [c for c in df.columns[2:] if sample_ids is None or c in sample_ids]
    sub = df[[gene_col, *sample_cols]].rename(columns={gene_col: "gene"}).set_index("gene")
    sub.index = sub.index.astype(str)
    if sub.index.has_duplicates:
        sub = sub.groupby(level=0).mean()
    values = np.nan_to_num(sub.to_numpy(dtype=np.float32).T, nan=0.0)  # missing CNA -> neutral 0
    return CNVMatrix(
        values=values,
        sample_ids=[str(s) for s in sub.columns],
        gene_names=[str(g) for g in sub.index],
    )


def barcode_to_sample(sample_ids: list[str]) -> list[str]:
    """Normalize TCGA barcodes to the sample level (TCGA-XX-XXXX-NN).

    The 4th field carries the 2-digit sample-type code plus an optional vial
    letter (e.g. ``01A``); the sample-level barcode keeps only the 2-digit code,
    so aliquot-/vial-level IDs collapse onto the same sample as the RNA/meth tables.
    """
    out = []
    for sid in sample_ids:
        parts = sid.split("-")
        if len(parts) >= 4:
            out.append("-".join(parts[:3]) + "-" + parts[3][:2])
        else:
            out.append(sid)
    return out


def align_to_genes(
    values: np.ndarray,
    genes: list[str],
    target_genes: list[str],
    *,
    fill_value: float = 0.0,
) -> np.ndarray:
    """Reorder/intersect CNV columns onto ``target_genes``; absent genes are filled.

    Returns a (n_samples x len(target_genes)) array. The fill_value (0.0 = neutral
    copy number on the GISTIC2 log-ratio scale) is used for genes not measured.
    """
    idx = {g: i for i, g in enumerate(genes)}
    out = np.full((values.shape[0], len(target_genes)), fill_value, dtype=np.float32)
    for j, tg in enumerate(target_genes):
        i = idx.get(tg)
        if i is not None:
            out[:, j] = values[:, i]
    return out


def harmonize_gene_level(
    values: np.ndarray, *, method: str = "zscore", eps: float = 1e-8
) -> np.ndarray:
    """Put CNV on a comparable internal scale per gene.

    ``zscore`` (default): per-gene standardization across samples. ``none``:
    pass-through (keep raw GISTIC2 scores). Cross-platform validity is a recorded
    limit, not a claim.
    """
    # Missing copy-number calls (NaN) -> neutral 0 before standardizing, so one
    # absent value can't poison a whole gene column (SNP6 CNA carries some NaN).
    values = np.nan_to_num(np.asarray(values, dtype=np.float32), nan=0.0)
    if method == "none":
        return values
    if method == "zscore":
        mu = values.mean(axis=0, keepdims=True)
        sd = values.std(axis=0, keepdims=True)
        return ((values - mu) / (sd + eps)).astype(np.float32)
    raise ValueError(f"unknown harmonize method {method!r}; use 'zscore' or 'none'")


def gene_overlap(genes_a: list[str], genes_b: list[str]) -> dict[str, int]:
    """Shared-gene stats between a CNV table and a target universe (honest coverage)."""
    a, b = set(genes_a), set(genes_b)
    return {"n_a": len(a), "n_b": len(b), "n_shared": len(a & b)}
