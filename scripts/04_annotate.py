"""
GSE303448 - Step 04: work out what each cluster is, and characterise the
IFNa-specific clusters (3 and 19).

Two questions:
  1. Which cell type is each of the 30 clusters?
  2. What defines clusters 3 and 19, which are ~97% IFNa?

Method:
  - rank_genes_groups: for each cluster, find genes expressed much more
    there than elsewhere. Those are the cluster's identity.
  - score_genes: for each known cell type, take its textbook marker set
    and give every cell a score for how strongly it expresses them.
    High score = that cell is probably that type.

REMINDER: sorting was 1:3:3:3 between the four gated populations, so
sizes BETWEEN broad types are artificial. Comparisons WITHIN a type
across conditions are fine.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

sc.settings.verbosity = 2
sc.settings.figdir = Path("figures")

lung = sc.read_h5ad("results/03_lung_clustered.h5ad")
print(f"Loaded: {lung.n_obs} cells, {lung.obs['leiden'].nunique()} clusters")

# The clustered object was cut down to 2000 variable genes. Marker genes
# may not be among them, so work from .raw, which kept every gene.
full = lung.raw.to_adata()
full.obs = lung.obs.copy()
full.obsm = lung.obsm.copy()

# --- 1. top genes per cluster ----------------------------------------------
print("\n=== finding marker genes per cluster ===")
sc.tl.rank_genes_groups(full, "leiden", method="wilcoxon", n_genes=25)

top = pd.DataFrame(full.uns["rank_genes_groups"]["names"]).head(12)
print("\n=== top 12 genes per cluster ===")
print(top.to_string())

top.to_csv("results/04_cluster_top_genes.csv")

# --- 2. score cells against known cell types -------------------------------
# Mouse lung marker sets. Each is a small set of genes that, together,
# identify one cell type.
MARKERS = {
    "Epithelial_general":  ["Epcam", "Krt8", "Krt18", "Cdh1"],
    "AT2_alveolar":        ["Sftpc", "Sftpb", "Lamp3", "Napsa"],
    "AT1_alveolar":        ["Ager", "Pdpn", "Hopx", "Akap5"],
    "Club_airway":         ["Scgb1a1", "Scgb3a2", "Cyp2f2"],
    "Ciliated_airway":     ["Foxj1", "Tppp3", "Dnah5"],
    "Endothelial":         ["Pecam1", "Cdh5", "Cldn5", "Tek"],
    "Fibroblast_mesench":  ["Col1a1", "Col1a2", "Pdgfra", "Dcn"],
    "SmoothMuscle":        ["Acta2", "Myh11", "Tagln"],
    "Immune_general":      ["Ptprc"],
    "AlveolarMacrophage":  ["Marco", "Siglecf", "Ear2", "Itgax"],
    "Monocyte_Macrophage": ["Cd68", "Csf1r", "Lyz2", "Fcgr1"],
    "Neutrophil":          ["S100a8", "S100a9", "Retnlg", "Mmp9"],
    "T_cell":              ["Cd3e", "Cd3d", "Trbc2", "Lck"],
    "NK_cell":             ["Nkg7", "Gzma", "Klrb1c", "Ncr1"],
    "B_cell":              ["Cd79a", "Ms4a1", "Ighm"],
    "DendriticCell":       ["Cd209a", "Flt3", "Batf3"],
    "Proliferating":       ["Mki67", "Top2a", "Ccnb1", "Birc5"],
    # MMTV-PyMT breast cancer cells: mammary epithelium, NOT lung.
    # Expect Epcam/Krt8/Krt18 positive but Sftpc/Scgb1a1 negative.
    "Mammary_cancer":      ["Krt14", "Krt5", "Wnt6", "Cldn3", "Cldn4"],
    # Interferon-stimulated genes - the response the paper is about
    "ISG_interferon":      ["Isg15", "Ifit1", "Ifit3", "Irf7", "Oasl2",
                            "Rsad2", "Stat1", "Bst2", "Cxcl10", "Usp18",
                            "Ifi27l2a", "Oas3", "Mx1", "Ifi44"],
}

print("\n=== scoring cells against marker sets ===")
for name, genes in MARKERS.items():
    present = [g for g in genes if g in full.var_names]
    missing = [g for g in genes if g not in full.var_names]
    if not present:
        print(f"  {name}: SKIPPED - no marker genes in this probe panel")
        continue
    if missing:
        print(f"  {name}: using {len(present)}/{len(genes)} "
              f"(not in panel: {', '.join(missing)})")
    sc.tl.score_genes(full, present, score_name=f"score_{name}")

score_cols = [c for c in full.obs.columns if c.startswith("score_")]

# Mean score per cluster - the highest score names the cluster
print("\n=== mean marker score per cluster ===")
by_cluster = full.obs.groupby("leiden", observed=True)[score_cols].mean().round(2)
print(by_cluster.to_string())
by_cluster.to_csv("results/04_cluster_scores.csv")

print("\n=== best-matching cell type per cluster ===")
cell_type_cols = [c for c in score_cols if c != "score_ISG_interferon"]
best = by_cluster[cell_type_cols].idxmax(axis=1).str.replace("score_", "")
sizes = full.obs["leiden"].value_counts().sort_index()
summary = pd.DataFrame({
    "n_cells": sizes,
    "best_match": best,
    "ISG_score": by_cluster["score_ISG_interferon"],
})
print(summary.to_string())
summary.to_csv("results/04_cluster_labels_draft.csv")

# --- 3. the IFNa-specific clusters -----------------------------------------
print("\n" + "=" * 70)
print("CLUSTERS 3 AND 19 - present almost only in IFNa mice")
print("=" * 70)

for cl in ["3", "19"]:
    if cl not in full.obs["leiden"].values:
        continue
    print(f"\n--- cluster {cl} ---")
    names = full.uns["rank_genes_groups"]["names"][cl][:20]
    lfc = full.uns["rank_genes_groups"]["logfoldchanges"][cl][:20]
    padj = full.uns["rank_genes_groups"]["pvals_adj"][cl][:20]
    print(pd.DataFrame({
        "gene": names,
        "log2FC": np.round(lfc, 2),
        "padj": ["%.1e" % p for p in padj],
    }).to_string(index=False))
    print(f"\nmarker scores for cluster {cl}:")
    print(by_cluster.loc[cl].sort_values(ascending=False).head(6).to_string())

# --- 4. ISG score across all clusters, split by condition ------------------
# The real question: does the interferon response appear in every cell
# type, or only some? This is the paper's claim about epithelium.
print("\n=== ISG score by cluster and condition ===")
isg = (
    full.obs.groupby(["leiden", "condition"], observed=True)["score_ISG_interferon"]
    .mean().unstack().round(3)
)
isg["difference"] = (isg["IFNa"] - isg["PBS"]).round(3)
print(isg.sort_values("difference", ascending=False).to_string())
isg.to_csv("results/04_isg_by_cluster_condition.csv")

# --- 5. plots ---------------------------------------------------------------
sc.pl.umap(full, color="score_ISG_interferon", cmap="viridis",
           title="Interferon-stimulated gene score",
           show=False, save="_04_isg_score.png")

sc.pl.umap(full, color=["score_Epithelial_general", "score_Endothelial",
                        "score_Immune_general", "score_Fibroblast_mesench"],
           ncols=2, show=False, save="_04_broad_types.png")

sc.pl.rank_genes_groups_dotplot(
    full, n_genes=3, groupby="leiden",
    show=False, save="_04_cluster_markers.png",
)

full.write("results/04_lung_annotated.h5ad")

print("\nSaved results/04_lung_annotated.h5ad")
print("Key outputs:")
print("  results/04_cluster_labels_draft.csv   - what each cluster probably is")
print("  results/04_isg_by_cluster_condition.csv - where the IFN response lands")
print("  figures/umap_04_isg_score.png         - IFN response on the map")