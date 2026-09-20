"""
GSE303448 - Step 03: filter, normalise, cluster.

What this script does, in plain terms:

  1. Throws away cells that are too empty to be informative, and cells
     that look like two cells stuck together.
  2. Puts every cell on a comparable scale, so a cell that happened to
     be measured more deeply does not look more active.
  3. Groups cells that behave similarly into clusters.
  4. Draws a map (UMAP) where similar cells sit near each other.

Filter choices come from figures/violin_01_qc.png:
  - min 500 genes: below this there is too little material to trust.
  - max 30000 counts: the spikes far above the typical ~7000 are almost
    certainly two cells captured as one (doublets).
  - NO mitochondrial filter: this experiment used a fixed probe panel
    that contains only 13 mitochondrial genes, so the usual dead-cell
    check reads ~0.1% for every cell regardless of health. Filtering on
    a metric that is structurally zero would do nothing.

REMEMBER: cells were sorted into 4 populations and pooled at a fixed
1:3:3:3 ratio. Cluster SIZES are set by that ratio, not by biology.
Never say "IFNa increased cell type X" from this data.
"""

import scanpy as sc
import numpy as np
from pathlib import Path

sc.settings.verbosity = 2
sc.settings.figdir = Path("figures")
sc.settings.set_figure_params(dpi=100, figsize=(6, 5))

MIN_GENES = 500
MAX_COUNTS = 30_000

lung = sc.read_h5ad("results/01_lung_merged_unfiltered.h5ad")
print(f"Loaded: {lung.n_obs} cells")

# --- 1. filter --------------------------------------------------------------
before = lung.n_obs

sc.pp.filter_cells(lung, min_genes=MIN_GENES)
lung = lung[lung.obs["total_counts"] < MAX_COUNTS].copy()

# a gene seen in fewer than 3 cells across 50k cells is noise
sc.pp.filter_genes(lung, min_cells=3)

print(f"\nFiltered: {before} -> {lung.n_obs} cells "
      f"({100 * lung.n_obs / before:.1f}% kept)")
print(lung.obs["sample"].value_counts().to_string())

# keep the raw counts before we transform anything - differential
# expression later needs the originals, not the scaled version
lung.layers["counts"] = lung.X.copy()

# --- 2. normalise -----------------------------------------------------------
# Put every cell on the same total (10,000 molecules), so cells measured
# more deeply do not look more active than cells measured less deeply.
sc.pp.normalize_total(lung, target_sum=1e4)

# Compress the range. Counts span 0 to thousands; without this, a handful
# of very highly expressed genes would dominate everything downstream.
sc.pp.log1p(lung)

lung.raw = lung  # keep the full normalised gene set for plotting later

# --- 3. pick the informative genes -----------------------------------------
# Most genes are expressed at a similar level in every cell and tell you
# nothing about which cell is which. Keep the ones that vary.
sc.pp.highly_variable_genes(lung, n_top_genes=2000, batch_key="sample")
print(f"\nHighly variable genes: {int(lung.var['highly_variable'].sum())}")

lung = lung[:, lung.var["highly_variable"]].copy()

sc.pp.scale(lung, max_value=10)

# --- 4. reduce and cluster --------------------------------------------------
# PCA: compress 2000 gene measurements per cell down to 50 numbers that
# capture most of the variation. Everything below works on those 50.
sc.tl.pca(lung, n_comps=50, svd_solver="arpack")

sc.pl.pca_variance_ratio(lung, n_pcs=50, show=False, save="_03_pca_variance.png")

# Build a network where each cell is linked to the cells most similar
# to it, then find densely connected communities in that network.
sc.pp.neighbors(lung, n_neighbors=15, n_pcs=30)
sc.tl.leiden(lung, resolution=0.5, key_added="leiden", flavor="igraph", n_iterations=2)

print(f"\nClusters found: {lung.obs['leiden'].nunique()}")
print(lung.obs["leiden"].value_counts().sort_index().to_string())

# UMAP: a 2D map for looking at. Distances on it are rough - trust the
# clusters, not the exact gaps between them.
sc.tl.umap(lung)

# --- 5. plots ---------------------------------------------------------------
sc.pl.umap(lung, color="leiden", legend_loc="on data",
           title="Clusters", show=False, save="_03_clusters.png")

sc.pl.umap(lung, color="condition", title="IFNa vs PBS",
           show=False, save="_03_condition.png")

sc.pl.umap(lung, color="sample", title="By mouse",
           show=False, save="_03_sample.png")

# --- 6. check whether the mice mixed ---------------------------------------
# If each mouse forms its own separate island, the differences between
# clusters are technical, not biological, and we need batch correction.
# If mice are spread across every cluster, we are fine.
print("\n=== how the 4 mice distribute across clusters (%) ===")
xtab = (
    lung.obs.groupby(["leiden", "sample"], observed=True).size().unstack(fill_value=0)
)
print((100 * xtab.T / xtab.sum(axis=1)).T.round(1).to_string())

lung.write("results/03_lung_clustered.h5ad")
print("\nSaved results/03_lung_clustered.h5ad")
print("Look at figures/umap_03_clusters.png and umap_03_sample.png next.")