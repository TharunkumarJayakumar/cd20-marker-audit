"""Record exactly which data and software produced the results.

The fingerprint is computed from the matrix contents (cell names, gene names
and raw counts), so anyone with the same public dataset can confirm a match.
"""
from __future__ import annotations

import hashlib
import importlib.metadata as metadata
import json
import platform
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from immune_atlas import expression_frame, load_data  # noqa: E402

RES = ROOT / "results"
RES.mkdir(exist_ok=True)

adata = load_data()
expr = expression_frame(adata)

digest = hashlib.sha256()
digest.update("\x00".join(adata.obs_names).encode())
digest.update(b"\x00")
digest.update("\x00".join(adata.var_names).encode())
digest.update(b"\x00")
digest.update(np.ascontiguousarray(expr.to_numpy(dtype=np.float32)).tobytes())

packages = ["scanpy", "anndata", "numpy", "pandas", "scipy",
            "matplotlib", "scikit-learn", "leidenalg", "umap-learn"]

record = {
    "dataset": "Scanpy pbmc68k_reduced — 700-cell reduced PBMC subset",
    "access": "scanpy.datasets.pbmc68k_reduced()",
    "origin": "Zheng et al. 2017, Nature Communications 8:14049, doi:10.1038/ncomms14049",
    "shape": list(adata.shape),
    "n_b_cells": int((adata.obs["bulk_labels"] == "CD19+ B").sum()),
    "content_sha256": digest.hexdigest(),
    "seeds": {"cluster": 0, "umap": 0},
    "runtime": {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    },
    "versions": {p: metadata.version(p) for p in packages},
}

out = RES / "dataset_provenance.json"
out.write_text(json.dumps(record, indent=2) + "\n")
print(f"Wrote {out}")
print(f"Dataset SHA-256: {record['content_sha256']}")
