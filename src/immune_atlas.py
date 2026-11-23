"""Reusable helpers for the PBMC atlas and the CD20 dropout audit.

The functions here deliberately separate two things that are often mixed up:

1. building an *unsupervised* atlas (no cell-type labels are used), and
2. auditing a marker gene *after* the atlas exists (labels are used only here).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# Marker panels used only for post-hoc naming of unsupervised clusters.
MARKER_PANELS = {
    "T cells": ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "TRAC"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD74", "CD19", "CD22"],
    "Monocytes": ["LYZ", "CST3", "S100A8", "S100A9", "FCN1", "VCAN"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "Dendritic": ["FCER1A", "CST3", "CLEC10A"],
}

B_PANEL = ["CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD74", "HLA-DRA", "PAX5"]
PLASMA_PANEL = ["PRDM1", "XBP1", "JCHAIN", "MZB1", "IGHG1", "IGKC"]


def load_data():
    """Load the 700-cell PBMC dataset shipped with Scanpy."""
    return sc.datasets.pbmc68k_reduced()


def expression_frame(adata) -> pd.DataFrame:
    """Return a dense gene x cell expression table (raw layer when present)."""
    source = adata.raw if adata.raw is not None else adata
    matrix = source.X
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(
        np.asarray(matrix), index=adata.obs_names, columns=source.var_names
    )


def cluster(adata, resolution: float = 1.0, seed: int = 0):
    """Neighbours -> Leiden -> UMAP using the supplied PCA, ignoring all labels."""
    out = adata.copy()
    sc.pp.neighbors(out, use_rep="X_pca", random_state=seed)
    try:  # Scanpy >= 1.10
        sc.tl.leiden(
            out, resolution=resolution, key_added="leiden",
            flavor="igraph", n_iterations=2, random_state=seed,
        )
    except TypeError:  # older Scanpy
        sc.tl.leiden(out, resolution=resolution, key_added="leiden", random_state=seed)
    sc.tl.umap(out, random_state=seed)
    return out


def annotate_by_markers(expr: pd.DataFrame, leiden: pd.Series) -> dict[str, str]:
    """Name each cluster by its highest mean marker-panel score."""
    scores = {}
    for name, genes in MARKER_PANELS.items():
        present = [g for g in genes if g in expr.columns]
        if not present:
            continue
        # mean of per-gene z-scores, so highly expressed genes don't dominate
        block = expr[present]
        z = (block - block.mean()) / block.std().replace(0, 1)
        scores[name] = z.mean(axis=1).groupby(leiden).mean()
    table = pd.DataFrame(scores)
    return table.idxmax(axis=1).to_dict()


def b_cell_mask(adata) -> pd.Series:
    """B cells from the coarse labels, used ONLY for the post-hoc marker audit."""
    labels = adata.obs["bulk_labels"].astype(str)
    # Handles both "B-cells" and "CD19+ B" naming across Scanpy versions.
    mask = labels.str.contains(r"\bB\b|B-cells|CD19", regex=True, case=False)
    mask &= ~labels.str.contains(r"CD4|CD8|CD14|CD34|NK|Dendritic", case=False)
    return mask


def panel_score(expr: pd.DataFrame, genes: list[str], cells=None) -> pd.Series:
    """Mean z-score of a gene panel across the given cells.

    Computed from the raw expression layer, because the reduced dataset's
    `X` matrix is centred/scaled and contains negative values.
    """
    rows = expr if cells is None else expr.loc[cells]
    present = [g for g in genes if g in rows.columns]
    block = rows[present]
    z = (block - block.mean()) / block.std().replace(0, 1)
    return z.mean(axis=1)


def detection_table(expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    present = [g for g in genes if g in expr.columns]
    return pd.DataFrame({
        "gene": present,
        "detected": [int((expr[g] > 0).sum()) for g in present],
        "fraction": [float((expr[g] > 0).mean()) for g in present],
    })


def dropout_audit(adata) -> dict:
    """Quantify how often the canonical CD20 gene MS4A1 reads zero in B cells."""
    expr = expression_frame(adata)
    mask = b_cell_mask(adata)
    b = expr.loc[mask[mask].index]
    targets = [g for g in B_PANEL + PLASMA_PANEL if g in b.columns]

    detection = detection_table(b, targets)
    ms4a1 = b["MS4A1"] if "MS4A1" in b.columns else None
    if ms4a1 is None:
        raise ValueError("MS4A1 not present in this dataset")

    zero = ms4a1 <= 0
    b_score = panel_score(expr, B_PANEL, b.index)
    plasma_score = panel_score(expr, PLASMA_PANEL, b.index)
    depth = (adata.obs.loc[b.index, "n_counts"].astype(float)
             if "n_counts" in adata.obs else b.sum(axis=1))

    # Empirical context: genes with similar mean expression have similar
    # detection rates. This avoids inventing a capture-efficiency parameter.
    means = b.mean(axis=1)
    detected = (b > 0).mean(axis=1)
    expressed = means[means > 0.01].index

    return {
        "expr": expr,
        "b_expr": b,
        "mask": mask,
        "detection": detection,
        "n_b": len(b),
        "n_ms4a1_zero": int(zero.sum()),
        "zero_fraction": float(zero.mean()),
        "zero_cells": b.index[zero],
        "pos_cells": b.index[~zero],
        "b_score": b_score,
        "plasma_score": plasma_score,
        "depth": depth,
        "means": means,
        "detected_genes": detected,
        "expressed_genes": expressed,
        "labels": adata.obs.loc[mask, "bulk_labels"].astype(str),
    }


def binomial_dropout(mean_molecules: float, capture: float) -> float:
    """P(see zero molecules) under simple binomial sampling.

    Illustrative only: 10 molecules at 10% capture gives ~0.35.
    """
    return float((1 - capture) ** mean_molecules)
