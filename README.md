**SARS-CoV-2 Infection Dynamics Analysis**

**Project Overview**

This project reproduces the trajectory analysis of SARS-CoV-2-infected human bronchial epithelial cells using single-cell RNA-seq (scRNA-seq) data. The analysis follows the workflow in the study: [PLOS Biology, 2020](https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3001143&utm_source=chatgpt.com).

We focus on:

- Neighborhood clustering
- Cell type identification
- Pseudotime ordering of differentiation
- Marker gene analysis (ACE2 and ENO2)

The dataset used was downloaded from GEO (GSE166766), consisting of four time points: **mock, 1 dpi, 2 dpi, 3 dpi**.

**Code Workflow Explained**

**1\. Setup and Dependencies**

The first step installs and imports the required Python packages:

!pip install scanpy anndata igraph celltypist decoupler fa2-modified louvain scvelo

import scanpy as sc

import anndata as ad

import numpy as np

import decoupler as dc

These packages handle scRNA-seq preprocessing, clustering, visualization, and cell type annotation.

**2\. Data Download**

Data for each time point is downloaded and stored in separate folders:

!mkdir -p GSM5082289_Mock GSM5082290_1dpi GSM5082291_2dpi GSM5082292_3dpi

\# Example for mock

!wget &lt;url_to_mock_barcodes&gt; -O GSM5082289_Mock/barcodes.tsv.gz

!wget &lt;url_to_mock_features&gt; -O GSM5082289_Mock/features.tsv.gz

!wget &lt;url_to_mock_matrix&gt; -O GSM5082289_Mock/matrix.mtx.gz

**Key idea:** load each dataset **separately** to reduce RAM usage.

**3\. Loading Data into AnnData**

Each time point is loaded individually:

mock_adata = sc.read_10x_mtx('/content/GSM5082289_Mock/')

dpi1_adata = sc.read_10x_mtx('/content/GSM5082290_1dpi/')

dpi2_adata = sc.read_10x_mtx('/content/GSM5082291_2dpi/')

dpi3_adata = sc.read_10x_mtx('/content/GSM5082292_3dpi/')

Each AnnData object stores gene expression counts and metadata.

**4\. Quality Control (QC)**

We calculate metrics like the proportion of mitochondrial genes, ribosomal genes, and hemoglobin genes:

adata.var\['MT'\] = adata.var_names.str.startswith("MT-")

adata.var\['RIBO'\] = adata.var_names.str.startswith(("RPS", "RPL"))

adata.var\['HB'\] = adata.var_names.str.contains(r"^HB\[^P\]")

sc.pp.calculate_qc_metrics(adata, qc_vars=\["MT", "RIBO", "HB"\], inplace=True, log1p=True)

Cells with low gene counts or high mitochondrial percentages are filtered out.

**5\. Normalization and Feature Selection**

We normalize counts, log-transform the data, and select highly variable genes:

adata.layers\["counts"\] = adata.X.copy()

sc.pp.normalize_total(adata)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=1000)

This prepares the data for dimensionality reduction and clustering.

**6\. Dimensionality Reduction and Clustering**

We perform PCA and UMAP, then identify clusters using Leiden algorithm:

sc.tl.pca(adata)

sc.pp.neighbors(adata)

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5)

Clusters are visualized using UMAP:

sc.pl.umap(adata, color=\["leiden_res0_5"\])

**7\. Marker Gene Analysis**

We visualize key genes ACE2 and ENO2 to understand infection dynamics:

sc.pl.umap(adata, color=\["ACE2", "ENO2"\])

- **ACE2**: SARS-CoV-2 entry receptor
- **ENO2**: marker for neuronal or stress-related response

**8\. Cell Type Annotation**

We annotate clusters using **Decoupler** and **PanglaoDB markers**:

markers = dc.op.resource(name="PanglaoDB", organism="human")

markers = markers\[markers\["organ"\] == 'Lungs'\]

markers = markers\[~markers.duplicated(\["cell_type", "genesymbol"\])\]

markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})\[\["source","target"\]\]

dc.mt.ulm(data=adata, net=markers, tmin=3)

score = dc.pp.get_obsm(adata, key="score_ulm")

gene_rank = dc.tl.rankby_group(score, groupby="leiden_res0_5", reference="rest", method="t-test_overestim_var")

top_cell_type_per_group = gene_rank.groupby('group')\['name'\].apply(lambda x: x.head(1))

adata.obs\["cell_type"\] = adata.obs\["leiden_res0_5"\].map(top_cell_type_per_group.to_dict())

UMAP plots with annotated cell types are generated for interpretation.

**Questions and Answers**

**1\. What cell types did you identify at the different stages of infection?**

- **Mock**: ciliated epithelial cells, basal cells, club cells
- **1 dpi**: increase in infected ciliated cells, some secretory cells
- **2 dpi**: more immune-related cells start appearing (early infection response)
- **3 dpi**: ciliated epithelial cells dominate ACE2-high clusters; immune response more pronounced

**2\. Why do these cell types correlate with COVID-19 infection?**

- SARS-CoV-2 infects cells expressing **ACE2**, mainly ciliated and secretory epithelial cells in the airway.
- Immune-related clusters appear later as the host responds to infection.
- The shift in cell type proportions mirrors viral infection dynamics in bronchial tissue.

**3\. Is ACE2 a good marker for tracking COVID-19 infection rate (based on this dataset)?**

- Yes, ACE2 expression correlates with infected cells (especially ciliated epithelial cells).
- UMAP visualizations show ACE2-high clusters expanding over dpi 1-3.

**4\. What is the difference between ENO2 and ACE2 as biomarkers?**

- **ACE2**: indicates susceptibility to viral entry; receptor for SARS-CoV-2
- **ENO2**: marks stress response or metabolic changes; not directly related to viral entry

**5\. Which cell cluster has the highest ACE2 at 3 dpi and biological interpretation?**

- **Ciliated epithelial cells** (identified via Leiden clustering + PanglaoDB annotation)
- **Biological meaning**: These cells are primary viral targets; high ACE2 explains why SARS-CoV-2 replicates efficiently in airway epithelium.
- Visual interpretation: UMAP shows ACE2-high spots localized to the ciliated cell cluster.

FIG 1: [VISUAL REPRESENTATION OF CELL CLUSTER WITH THE HIGHEST ACE2 AT 3 DPI.](https://drive.google.com/file/d/1GAhyQyKkaiGt34v6aYTJY1S1c0IS9w9a/view?usp=sharing)

**Notes on Code Execution**

- Each dataset is loaded **separately** to reduce RAM usage.
- QC metrics are calculated per dataset to ensure clean data.
- Layers in AnnData store original counts before normalization.
- Leiden clustering resolution and UMAP visualization can be adjusted to refine cluster separation.

**References**

- PLOS Biology. SARS-CoV-2 infection dynamics in human airway epithelial cells. 2020. [Link](https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3001143&utm_source=chatgpt.com)
- GEO Accession: [GSE166766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766&utm_source=chatgpt.com)
- Scanpy Documentation: <https://scanpy.readthedocs.io>
- Decoupler Documentation: <https://saezlab.github.io/decoupler>
