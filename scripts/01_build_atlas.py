"""Project 1: unsupervised atlas of 700 PBMCs. No labels are used here."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immune_atlas import annotate_by_markers, cluster, expression_frame, load_data, MARKER_PANELS

FIG = ROOT / "figures"
RES = ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

adata = load_data()
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
print("Label columns present but WITHHELD from clustering:",
      [c for c in adata.obs.columns if "label" in c.lower() or c == "louvain"])

atlas = cluster(adata, resolution=1.0, seed=0)
expr = expression_frame(atlas)
names = annotate_by_markers(expr, atlas.obs["leiden"])
atlas.obs["cell_type"] = atlas.obs["leiden"].map(names).astype(str)

# UMAP coloured by unsupervised clusters, then by marker-based annotation.
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, colour, title in zip(
    axes,
    ["leiden", "cell_type"],
    ["Unsupervised clusters (no labels)", "Clusters named by marker panels"],
):
    sc.pl.umap(atlas, color=colour, ax=ax, show=False, title=title,
               legend_loc="on data" if colour == "leiden" else "right margin")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
fig.tight_layout()
fig.savefig(FIG / "01_atlas_umap.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# Marker dot plot: evidence used to name clusters.
markers = list(dict.fromkeys(g for genes in MARKER_PANELS.values() for g in genes))
sc.pl.dotplot(atlas, var_names=[g for g in markers if g in atlas.var_names],
              groupby="cell_type", standard_scale="var", show=False)
plt.gcf().savefig(FIG / "02_marker_dotplot.png", dpi=180, bbox_inches="tight")
plt.close("all")

# Cross-tab is validation only, computed AFTER annotation.
agreement = pd.crosstab(atlas.obs["cell_type"], atlas.obs["bulk_labels"].astype(str))
agreement.to_csv(RES / "withheld_label_agreement.tsv", sep="\t")
pd.Series(names).rename("cell_type").to_csv(RES / "cluster_annotation.tsv", sep="\t")

print("\nCluster -> marker-based name:")
for cluster_id, name in sorted(names.items(), key=lambda kv: int(kv[0])):
    print(f"  {cluster_id}: {name}")
print("\nAgreement with withheld coarse labels (validation only):")
print(agreement.to_string())
print(f"\nSaved figures to {FIG}")
