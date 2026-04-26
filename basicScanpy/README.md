# Spatial Analysis: 10x Visium + Scanpy

Analyze 10x Genomics Visium Human Lymph Node. Map transcripts to H&E tissue. 

## Dataset
Human Lymph Node (V1_Human_Lymph_Node). 4,035 spots. 36,601 genes. H&E image included.

## Setup
```bash
pip install scanpy seaborn igraph leidenalg
```

## Workflow

### 1. Initialize
Set plot style and logging level.
```python
import scanpy as sc
import matplotlib.pyplot as plt

# Setup visuals and logging
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3 # High detail logging
```

### 2. Load & QC Metrics
Load Visium data. Flag mitochondrial genes. Calculate per-spot metrics.
```python
# Load dataset from 10x
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()

# Identify MT genes and calculate QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
```

### 3. Filtering
Remove low-quality spots and noise genes.
```python
# Filter spots by transcript count
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)

# Remove high MT spots (cell death indicator)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

# Remove genes seen in < 10 spots
sc.pp.filter_genes(adata, min_cells=10)
```

### 4. Normalization
Scale data for comparison. Select High Variable Genes (HVG).
```python
# Normalize depth and log transform
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# Find top 2000 genes for PCA
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```

### 5. Dimensionality & Clustering
Compute PCA, neighbors, and UMAP. Group spots via Leiden.
```python
# Build graph and clusters
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", flavor="igraph")
```

### 6. Spatial Mapping
Overlay clusters on H&E image.
```python
# Plot Leiden clusters on tissue
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

# Zoom into specific tissue regions
sc.pl.spatial(adata, img_key="hires", color="clusters",
              groups=["5", "9"], crop_coord=[7000, 10000, 0, 6000])
```

### 7. Differential Expression
Identify markers for spatial niches.
```python
# Statistical test for cluster markers
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")

# Visualize marker spatial spread (CR2 follicle marker)
sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])
```
