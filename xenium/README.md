# Squidpy Xenium Data Analysis

## Overview
This is a 10x Genomics Xenium spatial transcriptomics dataset of Human Lung Cancer providing true single-cell spatial resolution. [Source dataset](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html).

## Setup
Install required packages for spatial data and download the dataset.

```python
# Install required spatial libraries
!pip install -q scanpy squidpy igraph leidenalg scikit-image dask spatialdata spatialdata-io spatialdata-plot seaborn

# Download the official Xenium Human Lung test dataset from 10x Genomics
!curl -O https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Lung_2fov/Xenium_V1_human_Lung_2fov_outs.zip

# Extract it into a folder named Xenium
!unzip -q -o Xenium_V1_human_Lung_2fov_outs.zip -d Xenium
```

## Steps

Load Xenium data, convert to Zarr format, and extract the count matrix.
```python
# Import core packages
import spatialdata as sd
import spatialdata_plot
from spatialdata_io import xenium
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns

# Load raw Xenium data, save to Zarr, and read back
xenium_path = "./Xenium"
zarr_path = "./Xenium.zarr"
sdata_raw = xenium(xenium_path, cells_as_circles=True)
sdata_raw.write(zarr_path)

# Extract the AnnData object containing expression and spatial coordinates
sdata = sd.read_zarr(zarr_path)
adata = sdata.tables["table"]
```

Calculate quality control metrics and plot basic distributions.
```python
# Calculate basic QC metrics in place
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

# Calculate percentage of total counts from negative controls
cprobes = adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
cwords  = adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100

print(f"Negative DNA probe %: {cprobes:.2f}")
print(f"Negative decoding %:  {cwords:.2f}")

# Plot transcripts, genes, cell area, and nucleus ratio
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0]).set_title("Total transcripts")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1]).set_title("Unique transcripts")
sns.histplot(adata.obs["cell_area"], kde=False, ax=axs[2]).set_title("Cell area")
sns.histplot(adata.obs["nucleus_area"] / adata.obs["cell_area"], kde=False, ax=axs[3]).set_title("Nucleus ratio")
plt.savefig("xenium_qc_histograms.png", bbox_inches="tight", dpi=300)
plt.show()
```

Filter low quality cells and genes based on metrics.
```python
# Remove cells with too few counts and genes found in too few cells
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)
```

Preprocess counts, build neighbor graph, and cluster cells.
```python
# Keep raw counts, normalize, log transform, and run PCA
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)

# Build neighbor graph in PCA space and generate Leiden clusters
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor="igraph", directed=False, n_iterations=2)
```

Visualize transcriptional clusters on a UMAP.
```python
# Plot UMAP colored by QC metrics and Leiden clusters
sc.pl.umap(
    adata,
    color=["total_counts", "n_genes_by_counts", "leiden"],
    wspace=0.4
)
plt.savefig("xenium_umap_clusters.png", bbox_inches="tight", dpi=300)
plt.show()
```

Build a spatial neighborhood graph using true cellular geometry.
```python
# Connect each cell to nearest spatial neighbors using Delaunay triangulation
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```

Determine the structural importance of each cell cluster.
```python
# Compute centrality scores for Leiden clusters
sq.gr.centrality_scores(adata, cluster_key="leiden")

# Plot degree, closeness, and average centrality
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))
plt.savefig("xenium_centrality_scores.png", bbox_inches="tight", dpi=300)
plt.show()
```

Find out which clusters physically colocalize at different distances.
```python
# Subsample data to 50% to speed up computation
sdata.tables["subsample"] = sc.pp.subsample(adata, fraction=0.5, copy=True)
adata_subsample = sdata.tables["subsample"]

# Calculate co-occurrence across distance bins
sq.gr.co_occurrence(adata_subsample, cluster_key="leiden")

# Plot spatial co-occurrence for cluster 1 against all others
sq.pl.co_occurrence(adata_subsample, cluster_key="leiden", clusters="1", figsize=(10, 10))
plt.savefig("xenium_co_occurrence.png", bbox_inches="tight", dpi=300)
plt.show()
```

Calculate and visualize neighborhood enrichment between clusters.
```python
# Compute neighborhood enrichment scores
sq.gr.nhood_enrichment(adata, cluster_key="leiden")

# Plot enrichment heatmap and spatial scatter side-by-side
fig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(adata, cluster_key="leiden", figsize=(8, 8), title="Neighborhood enrichment", ax=ax[0])
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])
plt.savefig("xenium_nhood_enrichment.png", bbox_inches="tight", dpi=300)
plt.show()
```

Identify genes that are spatially clustered across the tissue.
```python
# Re-compute spatial graph on subsample and calculate Moran's I
sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(adata_subsample, mode="moran", n_perms=100, n_jobs=1)

# Display top spatially variable genes
display(adata_subsample.uns["moranI"].head())

# Plot selected spatially autocorrelated genes
sq.pl.spatial_scatter(
    adata_subsample,
    library_id="spatial",
    color=["AREG", "MET"],
    shape=None,
    size=2,
    img=False
)
plt.savefig("xenium_spatial_scatter_genes.png", bbox_inches="tight", dpi=300)
plt.show()
```

Render highly autocorrelated genes directly over the tissue morphology image.
```python
# Render gene expression using static spatialdata plotting
gene_names = ["AREG", "MET"]
for name in gene_names:
    sdata.pl.render_images("morphology_focus").pl.render_shapes(
        "cell_circles", color=name, table_name="table"
    ).pl.show(
        title=f"{name} expression over Morphology image", 
        coordinate_systems="global", 
        figsize=(10, 5)
    )
```

Launch an interactive viewer for deep data exploration.
```python
# Open napari interactive viewer (works best locally)
from napari_spatialdata import Interactive
Interactive(sdata)
```

## Results

| Output Plot | Description |
| :--- | :--- |
| `xenium_qc_histograms.png` | Four histograms showing QC metric distributions. |
| `xenium_umap_clusters.png` | UMAP scatter plots colored by counts, genes, and Leiden clusters. |
| `xenium_centrality_scores.png` | Scatter plots showing degree, closeness, and average centrality by cell cluster. |
| `xenium_co_occurrence.png` | Line plot showing probability of cluster 1 co-occurring with other clusters across distances. |
| `xenium_nhood_enrichment.png` | Heatmap of neighborhood enrichment next to a spatial scatter plot. |
| `xenium_spatial_scatter_genes.png` | 2D scatter plots mapping top spatially variable genes (AREG, MET). |
| Inline Rendering | Morphology image overlay with AREG and MET expression using spatialdata-plot. |

## Reference
Tutorial source: [Analyze Xenium data — Squidpy](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html)
