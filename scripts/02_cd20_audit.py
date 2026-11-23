"""Project 2: audit the canonical B-cell marker CD20 (MS4A1), cell by cell."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immune_atlas import (
    B_PANEL, PLASMA_PANEL, binomial_dropout, dropout_audit, load_data,
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
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("Why 34 of 95 B cells read zero for CD20 (MS4A1)",
             fontsize=15, fontweight="bold")

# A: detection rates across marker panels
ax = axes[0, 0]
det = audit["detection"].set_index("gene")["fraction"]
order = [g for g in B_PANEL + PLASMA_PANEL if g in det.index]
colours = ["#E45756" if g == "MS4A1" else "#22A7B8" for g in order]
ax.barh(order[::-1], det[order][::-1], color=colours[::-1])
ax.set_xlim(0, 1.02)
ax.set_xlabel("Fraction of B cells with a non-zero count")
ax.set_title("A. One gene is not enough")

# B: detection rate vs mean expression (empirical dropout curve)
ax = axes[0, 1]
ax.scatter(audit["means"][audit["expressed_genes"]],
           audit["detected_genes"][audit["expressed_genes"]],
           s=9, alpha=0.35, color="#4C78A8")
if "MS4A1" in audit["means"].index:
    mx, my = audit["means"]["MS4A1"], audit["detected_genes"]["MS4A1"]
    ax.scatter([mx], [my], s=130, color="#E45756", edgecolor="black", zorder=5)
    ax.annotate("MS4A1", (mx, my), xytext=(8, -14), textcoords="offset points",
                fontsize=10, fontweight="bold", color="#B91C1C")
ax.set_xlabel("Mean expression in B cells")
ax.set_ylabel("Fraction detected")
ax.set_title("B. Sampling predicts this zero — but is it the cause?")

# C: what the zero cells DO express
ax = axes[1, 0]
genes_c = [g for g in ["CD79A", "MS4A1", "MZB1"] if g in audit["b_expr"].columns]
zero_rate = [(audit["b_expr"].loc[zero, g] > 0).mean() for g in genes_c]
pos_rate = [(audit["b_expr"].loc[pos, g] > 0).mean() for g in genes_c]
xpos = np.arange(len(genes_c))
ax.bar(xpos - 0.19, zero_rate, 0.38, label="MS4A1 = 0", color="#E45756")
ax.bar(xpos + 0.19, pos_rate, 0.38, label="MS4A1 > 0", color="#22A7B8")
ax.set_xticks(xpos, genes_c)
ax.set_ylim(0, 1.05)
ax.set_ylabel("Fraction of cells with non-zero count")
ax.set_title("C. CD79A stays on; MZB1 is enriched")
ax.legend(fontsize=9)

# D: sequencing depth by MS4A1 status
ax = axes[1, 1]
groups = [audit["depth"].loc[zero], audit["depth"].loc[pos]]
ax.boxplot(groups, tick_labels=["MS4A1 = 0", "MS4A1 > 0"], patch_artist=True,
           boxprops=dict(facecolor="#FDE68A"), medianprops=dict(color="#B45309"))
ax.set_ylabel("UMI counts per cell")
ax.set_title("D. More RNA, not less — argues against dropout")

fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(FIG / "03_cd20_dropout.png", dpi=180, bbox_inches="tight")
print(f"\nSaved figure to {FIG / '03_cd20_dropout.png'}")
