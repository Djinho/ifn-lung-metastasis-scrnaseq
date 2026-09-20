# Interferon-α and the alveolar environment — a single-cell reproduction

A reproduction of the published single-cell RNA-seq analysis from
**GSE303448**, re-run from the deposited count matrices.

No new biological claims are made here. The aim was to work through a
published dataset end to end and see whether the reported cell-state
changes fall out of an independent analysis.

**Source study:** Farias A, Bridgeman V, Rodrigues F, Puttur F, Owen A,
Ruhland S, Ferreira R, Mack M, Malanchi I, Johansson C. *Type I
interferons induced upon respiratory viral infection impair lung
metastatic initiation.* PNAS, 21 April 2026.
DOI: [10.1073/pnas.2412919123](https://doi.org/10.1073/pnas.2412919123)

**Data:** [GSE303448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303448)


## The question

Breast cancer commonly spreads to the lung. The source study asks whether
a lung already responding to type I interferon — the signal released
during a respiratory viral infection — is a less hospitable place for
arriving cancer cells to seed.

The experiment: mice received intranasal IFN-α or PBS, then 600,000
primary MMTV-PyMT breast cancer cells intravenously 18 h later. Lungs were
harvested 24 h after that, dissociated, sorted, and sequenced.


## Data and design

| | |
|---|---|
| Organism | *Mus musculus* |
| Samples | 4 — IFNa1, IFNa2, PBS1, PBS2 (2 conditions × 2 mice) |
| Chemistry | 10x Fixed RNA Profiling (Flex), Chromium Mouse Transcriptome Probe Set v1.0 |
| Reference | refdata-gex-mm10-2020-A |
| Features | 19,059 (probe panel, not whole transcriptome) |
| Cells loaded | 50,313 |
| Cells after QC | 46,533 (92.5%) |

The four samples were multiplexed in a single library (probe barcodes
BC001–BC004) and are deposited pre-demultiplexed.


## Three things about this dataset that shape the analysis

**1. Cell-type proportions are set by the sorting, not by biology.**
Cells were FACS-sorted into four populations — cancer (GFP⁺), immune
(CD45⁺), epithelial (CD326⁺) and mesenchymal (CD45⁻CD31⁻) — then pooled at
a fixed **1:3:3:3 ratio** before sequencing. Cluster sizes therefore say
nothing about lung composition. All comparisons here are of cell *states*
within a cell type across conditions, never of abundance between types.

**2. The usual dead-cell filter does not work here.** The probe panel
contains only 13 mitochondrial genes, so mitochondrial fraction reads
~0.1% for every cell regardless of health. Filtering on it would remove
nothing. Cells were filtered on gene count (≥500) and total counts
(<30,000, to drop likely doublets) instead.

**3. The cancer cells carry no marker in the reference.** Although cancer
cells were sorted on GFP, the reference is stock mm10 with no transgene,
so neither GFP nor PyMT was quantified. The injected cells were identified
by expression profile instead — Epcam⁺ Krt8⁺ Krt14⁺ but Sftpc⁻ Scgb1a1⁻,
i.e. mammary rather than lung epithelium.


## Pipeline

```
scripts/
  01_inspect.py         inspect raw files; confirm shared feature space
  02_load.py            load 4 matrices -> AnnData, QC metrics, QC plots
  03_filter_cluster.py  filter, normalise, HVG, PCA, neighbours, Leiden, UMAP
  04_annotate.py        marker ranking, cell-type scoring, ISG scoring
```

Run from the project root:

```bash
python scripts/01_inspect.py
python scripts/02_load.py
python scripts/03_filter_cluster.py
python scripts/04_annotate.py
```

Parameters: 2,000 highly variable genes (batch-aware by sample), 50 PCs
computed / 30 used for the neighbour graph, k = 15, Leiden resolution 0.5
-> 30 clusters.


## What the reproduction shows

**Cell types recovered.** Clusters annotated by scoring against canonical
mouse-lung marker sets: alveolar type 2 and type 1 epithelium, club and
ciliated airway cells, endothelium, fibroblasts, smooth muscle, alveolar
and monocyte-derived macrophages, neutrophils, T, B and NK cells, and a
mammary-epithelial population corresponding to the injected MMTV-PyMT
cells.

**The alveolar epithelium separates by condition.** Two clusters score as
AT2 cells. Cluster 2 (n = 6,072) carries a low interferon-stimulated gene
score (−0.15) and is drawn mostly from control mice. Cluster 3
(n = 2,882) carries a high ISG score (0.89) and is ~97% IFN-α, with both
treated mice contributing equally (49.2% / 48.3%). This is consistent with
the source study's account of an ISG-driven alveolar response.

**The response is lung-wide but uneven.** Every cluster shows a positive
IFN-α − PBS shift in ISG score. By that measure the largest shifts in this
reanalysis fall in the monocyte–macrophage clusters (1.13, 1.01), followed
by alveolar epithelium (0.93, 0.92) and fibroblasts (0.88); endothelial
clusters sit lower (0.38–0.74).

**A proliferating NK population appears only under IFN-α.** Cluster 19
(n = 175) is Nkg7⁺ Gzma⁺ Gzmb⁺ Klrk1⁺ with strong cell-cycle expression
(Mki67, Top2a, Ccna2), and is ~95% IFN-α.

Per-cluster tables are in `results/`; figures in `figures/`.


## Limitations

- Single-cell data from a probe panel, not whole-transcriptome capture:
  absence of a gene may reflect absence of a probe.
- Two biological replicates per condition. Cluster-level differences are
  descriptive; formal differential expression should treat the mouse, not
  the cell, as the unit of replication.
- No doublet-detection tool was run; a total-count ceiling was used as a
  proxy.
- Cell-type labels come from marker-set scoring and have not been
  cross-checked against a reference atlas.


## Environment

Python 3.11, scanpy, anndata, scipy, pandas, numpy.

```bash
pip install scanpy anndata pandas numpy scipy python-igraph leidenalg
```

Raw data is not tracked in this repository. Download `GSE303448_RAW.tar`
from GEO and extract into `data/raw/`.


## Author

Djinho Itshary — [github.com/Djinho](https://github.com/Djinho)
