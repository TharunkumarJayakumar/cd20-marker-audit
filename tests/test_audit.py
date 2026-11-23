"""Tests that pin the biological expectations of the CD20 audit."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from immune_atlas import (  # noqa: E402
    b_cell_mask, binomial_dropout, dropout_audit, load_data,
)


def test_dataset_shape_and_b_cell_count():
    adata = load_data()
    assert (adata.n_obs, adata.n_vars) == (700, 765)
    assert int(b_cell_mask(adata).sum()) == 95


def test_b_cell_mask_excludes_other_lineages():
    adata = load_data()
    mask = b_cell_mask(adata)
    labels = adata.obs.loc[mask, "bulk_labels"].astype(str)
    assert not labels.str.contains(r"CD4|CD8|CD14|NK|Dendritic", case=False).any()


def test_expected_ms4a1_zero_count():
    audit = dropout_audit(load_data())
    assert audit["n_b"] == 95
    assert audit["n_ms4a1_zero"] == 34
    assert round(audit["zero_fraction"], 3) == 0.358


def test_zero_cells_keep_b_lineage_and_show_plasma_marker():
    audit = dropout_audit(load_data())
    expr = audit["b_expr"]
    zero, pos = audit["zero_cells"], audit["pos_cells"]
    assert int((expr.loc[zero, "CD79A"] > 0).sum()) == 33
    assert int((expr.loc[zero, "MZB1"] > 0).sum()) == 30
    assert int((expr.loc[pos, "CD79A"] > 0).sum()) == 61
    assert audit["depth"].loc[zero].mean() > audit["depth"].loc[pos].mean()


def test_binomial_illustration():
    assert abs(binomial_dropout(10, 0.10) - 0.9 ** 10) < 1e-12
