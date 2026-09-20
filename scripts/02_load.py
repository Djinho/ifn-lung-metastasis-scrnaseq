"""
GSE303448 - Step 02: load the four samples into one AnnData object.


IFNa vs PBS treated mouse lung, 24h after IV injection of MMTV-PyMT cells.

DESIGN WARNING - read before interpreting anything:
  Cells were FACS-sorted into four populations and then POOLED AT A FIXED
  1:3:3:3 RATIO before sequencing:
      cancer (GFP+) : immune (CD45+) : epithelial (CD326+) : mesenchymal
  Cell type proportions are therefore an artefact of the sorting gate,
  not lung biology. Never interpret cluster sizes as composition.
  Only compare expression WITHIN a cell type, across conditions.

Design is 2 conditions x 2 biological replicates (one mouse each).
Keep 'replicate' attached - it is the unit for differential expression,
not the individual cell.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io
import scipy.sparse as sp
from pathlib import Path

sc.settings.verbosity = 2
sc.settings.figdir = Path("figures")

RAW = Path("data/raw")
SAMPLES = ["IFNa1", "IFNa2", "PBS1", "PBS2"]

Path("figures").mkdir(exist_ok=True)
Path("results").mkdir(exist_ok=True)


def load_sample(sample: str) -> ad.AnnData:
    """Read one 10x mtx triplet into AnnData (cells x genes)."""
    # 10x stores genes x cells, AnnData wants cells x genes -> transpose
    X = scipy.io.mmread(RAW / f"GSM9126968_{sample}_matrix.mtx.gz")
    X = sp.csr_matrix(X.T)

    barcodes = pd.read_csv(
        RAW / f"GSM9126968_{sample}_barcodes.tsv.gz", sep="\t", header=None
    )
    features = pd.read_csv(
        RAW / f"GSM9126968_{sample}_features.tsv.gz", sep="\t", header=None
    )

    obs = pd.DataFrame(index=barcodes[0].astype(str).values)
    var = pd.DataFrame(index=features[1].astype(str).values)
    var["gene_id"] = features[0].astype(str).values
    if features.shape[1] >= 3:
        var["feature_type"] = features[2].astype(str).values

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()

    adata.obs["sample"] = sample
    adata.obs["condition"] = "IFNa" if sample.startswith("IFNa") else "PBS"
    adata.obs["replicate"] = sample[-1]

    print(f"{sample}: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


print("\n=== loading ===")
adatas = {s: load_sample(s) for s in SAMPLES}

# --- merge on shared genes only --------------------------------------------
shared = set(adatas[SAMPLES[0]].var_names)
for s in SAMPLES[1:]:
    shared &= set(adatas[s].var_names)
shared = sorted(shared)
print(f"\nGenes shared across all 4 samples: {len(shared)}")

adatas = {s: a[:, shared].copy() for s, a in adatas.items()}

lung = ad.concat(
    list(adatas.values()),
    label="batch",
    keys=SAMPLES,
    index_unique="-",
)
lung.obs_names_make_unique()

print(f"\nMerged: {lung.n_obs} cells, {lung.n_vars} genes")
print(lung.obs["sample"].value_counts().to_string())
print(lung.obs["condition"].value_counts().to_string())

# --- QC metrics -------------------------------------------------------------
# mouse mitochondrial genes are lowercase: mt-Nd1, mt-Co1 ...
# (the human uppercase "MT-" would match nothing and silently give 0%)
lung.var["mt"] = lung.var_names.str.startswith("mt-")
print(f"\nMitochondrial genes found: {int(lung.var['mt'].sum())}")

lung.var["ribo"] = lung.var_names.str.match(r"^Rp[sl]")

sc.pp.calculate_qc_metrics(
    lung, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
)

print("\n=== QC summary by sample ===")
summary = (
    lung.obs.groupby("sample", observed=True)
    .agg(
        cells=("n_genes_by_counts", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    .round(1)
)
print(summary.to_string())

# --- QC plots: look before filtering ----------------------------------------
sc.pl.violin(
    lung,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="sample",
    jitter=False,
    multi_panel=True,
    show=False,
    save="_01_qc.png",
)

sc.pl.scatter(
    lung, x="total_counts", y="pct_counts_mt", color="sample",
    show=False, save="_01_counts_vs_mt.png",
)

lung.write("results/01_lung_merged_unfiltered.h5ad")

print("\nSaved results/01_lung_merged_unfiltered.h5ad")
print("Look at figures/violin_01_qc.png and set filter cutoffs from what")
print("you see - do not copy thresholds from a tutorial.")