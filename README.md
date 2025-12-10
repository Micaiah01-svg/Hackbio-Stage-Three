**SARS-CoV-2 Infection Dynamics - Single-Cell RNA-seq Analysis**

**Reference:** Ravindra _et al._, _PLOS Biology_ (2021) - _Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium_  
<https://doi.org/10.1371/journal.pbio.3001143>

**Data:** GEO GSE166766 (10x mtx/tsv)  
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766>

**🎯 Project Objective**

Reproduce **neighborhood clustering**, **cell-type annotation**, and **trajectory/pseudotime analysis** of SARS-CoV-2 infection in HBECs (human bronchial epithelial cells). Focus on **mock, 1 dpi, 2 dpi, and 3 dpi**, targeting reproduction of figures similar to **1G(i-iii), 3A/B, 4A/B**.

Special emphasis:

- Comprehensive **pseudotime analysis** to capture cell state transitions over infection.
- Clear annotation of **ionocytes** as a rare but biologically relevant population.

**🔬 Biological Overview**

- Eight epithelial cell types identified: **ciliated, basal, club, BC/club (basal→club), goblet, neuroendocrine, ionocytes, tuft**.
- **Ciliated cells** are the primary targets of early infection. Infection spreads to **basal, club, BC/club**, and **ionocytes** at later stages.
- **ACE2 expression** indicates cell-type susceptibility but is a poor per-cell predictor of infection.
- **ENO2** marks neuroendocrine/metabolic state changes and cellular stress, offering complementary insights into infected vs. bystander cells.

**✅ Analysis Workflow (Modular Approach Recommended)**

**1\. Data Loading**

def load_adata(path, var_names='gene_symbols'):

"""Load 10x matrix into AnnData object with error handling."""

import scanpy as sc

import os

if not os.path.exists(path):

raise FileNotFoundError(f"Data path {path} not found.")

adata = sc.read_10x_mtx(path, var_names=var_names, cache=True)

return adata

- Load each condition (mock, 1 dpi, 2 dpi, 3 dpi) separately.
- Avoid system-specific download commands; use Python requests or GEOquery to fetch data.

**2\. Quality Control**

- Filter cells with **<200 genes** and genes present in **<3 cells**.
- Remove cells with **\>15% mitochondrial reads**.
- Use modular QC functions with inline comments for clarity.

**3\. Normalization & Feature Selection**

- Normalize total counts → CPM → log1p.
- Identify **highly variable genes** (HVGs, n=2000).
- Encapsulate preprocessing in reusable functions:

def preprocess_adata(adata):

import scanpy as sc

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

return adata\[:, adata.var\['highly_variable'\]\]

**4\. Dimensionality Reduction & Batch Correction**

- **PCA** → **BB-kNN** (preserves neighborhood across batches).
- **UMAP & PHATE** embeddings for visualization and pseudotime.

**5\. Clustering & Annotation**

- Louvain clustering → manual annotation using canonical markers:
  - FOXJ1 → ciliated
  - KRT5 → basal
  - SCGB1A1 → club
  - MUC5AC → goblet
  - PAX6 / ENO2 → neuroendocrine
  - POU2F3 → tuft
  - FOXI1 → ionocytes
- Modular functions for annotation allow future reuse.

**6\. Infection Flagging**

- Cells with **≥10 viral transcripts** flagged as infected.
- ACE2 vs ENO2 expression analyzed per cluster for mechanistic insights.

**7\. Pseudotime / Trajectory Analysis**

- Compute diffusion pseudotime (DPT) using PHATE embeddings.
- Overlay **pseudotime on UMAP/PHATE** to interpret infection progression.
- Optional: RNA velocity with **scVelo** if spliced/unspliced counts are available.

def compute_pseudotime(adata, n_dcs=10):

import scanpy as sc

sc.tl.diffmap(adata)

sc.tl.dpt(adata, n_dcs=n_dcs)

return adata

**🔍 Key Observations**

- **Cell types at each infection stage:**
  - **Mock:** All eight types, including **ionocytes**, present; no viral signal.
  - **1 dpi:** Infection primarily in **ciliated cells**.
  - **2-3 dpi:** Spread to **basal, club, BC/club, and ionocytes**.
- **Correlation with infection:**
  - Tropism is **cell type-driven** via ACE2 and TMPRSS2.
  - Pseudotime analysis reveals transitions in cell states and temporal infection patterns.
- **ACE2 vs ENO2:**
  - **ACE2:** Entry receptor; indicates cell susceptibility.
  - **ENO2:** Metabolic/neuroendocrine stress marker; highlights state changes in infected/bystander cells.
- **ACE2 enrichment (3 dpi):**
  - Highest in **ciliated cells**, overlaps with infected populations.
  - Suggests **primary entry sites** and persistent susceptibility.

**📦 Deliverables**

- README.md - this improved overview
- notebook_clustering.ipynb - modular, annotated preprocessing, clustering, and annotation
- notebook_pseudotime.ipynb - detailed trajectory/pseudotime analysis
- figures/ - UMAP/PHATE, heatmaps, violin plots (reproducing 1G, 3A/B, 4A/B)
- results/ - per-cluster ACE2 expression & infected cell counts

**📚 References**

- Ravindra NG, Alfajaro MM, Gasque V, et al. _Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium._ PLOS Biology (2021). [DOI](https://doi.org/10.1371/journal.pbio.3001143)
- GEO: GSE166766 - 10x single-cell RNA-seq [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766)
- ENO2 / NSE context: Cione _et al._, 2021; Moreno Jr _et al._, 2022
