"""Project 2: audit the canonical B-cell marker CD20 (MS4A1), cell by cell."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immune_atlas import (
    B_PANEL, PLASMA_PANEL, audit_figure, binomial_dropout, dropout_audit, load_data,
)

FIG = ROOT / "figures"
RES = ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

adata = load_data()
audit = dropout_audit(adata)

print(f"B cells: {audit['n_b']}")
print(f"MS4A1-zero B cells: {audit['n_ms4a1_zero']} ({audit['zero_fraction']:.1%})")
print("\nDetection by gene in B cells:")
print(audit["detection"].to_string(index=False))

zero = audit["zero_cells"]
pos = audit["pos_cells"]
print("\nMean B-cell panel score — MS4A1-zero vs MS4A1-positive:")
print(f"  zero: {audit['b_score'].loc[zero].mean():+.3f}")
print(f"  pos : {audit['b_score'].loc[pos].mean():+.3f}")
print("\nMean UMI counts:")
print(f"  zero: {audit['depth'].loc[zero].mean():.0f}")
print(f"  pos : {audit['depth'].loc[pos].mean():.0f}")
print("\nMZB1 (plasma signature) detection:")
print(f"  zero: {int((audit['b_expr'].loc[zero, 'MZB1'] > 0).sum())}/{len(zero)}")
print(f"  pos : {int((audit['b_expr'].loc[pos, 'MZB1'] > 0).sum())}/{len(pos)}")
print("\nCD79A (B-lineage) detection:")
print(f"  zero: {int((audit['b_expr'].loc[zero, 'CD79A'] > 0).sum())}/{len(zero)}")
print(f"  pos : {int((audit['b_expr'].loc[pos, 'CD79A'] > 0).sum())}/{len(pos)}")
print(f"\nIllustrative binomial dropout: 10 molecules at 10% capture -> "
      f"{binomial_dropout(10, 0.10):.1%}")
print("\nInterpretation: the zero cells retain CD79A, have HIGHER UMI counts and")
print("are enriched for MZB1 -> mostly genuine CD20-low plasmablasts, not dropout.")

# --- tables ------------------------------------------------------------------
audit["detection"].to_csv(RES / "b_cell_gene_detection.tsv", sep="\t", index=False)
summary = audit["b_score"].groupby(audit["b_expr"]["MS4A1"] > 0).agg(["count", "mean", "std"])
summary.to_csv(RES / "ms4a1_zero_vs_positive.tsv", sep="\t")

# --- figure ------------------------------------------------------------------
fig = audit_figure(audit)
fig.savefig(FIG / "03_cd20_dropout.png", dpi=180, bbox_inches="tight")
print(f"\nSaved figure to {FIG / '03_cd20_dropout.png'}")
