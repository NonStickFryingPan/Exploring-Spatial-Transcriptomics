# Squidpy MERFISH Data Analysis

## Overview
This is a pre-processed 3D MERFISH dataset of the mouse hypothalamic preoptic region from [Moffitt et al., 2018](https://science.sciencemag.org/content/362/6416/eaau5324). 

## Setup
Install required packages for Colab environment.

```python
# Install spatial analysis tools
!pip install -q numpy pandas anndata scanpy squidpy matplotlib
```

## Steps

Import packages and load the dataset.
```python
# Load required libraries
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

# Fetch the pre-processed MERFISH data
adata = sq.datasets.merfish()
```

Visualize the cell classes across all 3D slices.
```python
# Plot all cell classes in 3D space
sc.pl.embedding(
    adata, 
    basis="spatial3d", 
    projection="3d", 
    color="Cell_class", 
    show=False
)
plt.savefig("merfish_3d_cell_class.png", bbox_inches="tight", dpi=300)
plt.show()
```

Extract and visualize a single 2D slice.
```python
# Filter data for Bregma -9 and plot 2D scatter
sq.pl.spatial_scatter(
    adata[adata.obs.Bregma == -9], 
    shape=None, 
    color="Cell_class", 
    size=1
)
plt.savefig("merfish_2d_slice_bregma_minus9.png", bbox_inches="tight", dpi=300)
plt.show()
```

Compute and plot neighborhood enrichment in 3D.
```python
# Build spatial graph using 3D coordinates
sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial3d")

# Calculate neighborhood enrichment score
sq.gr.nhood_enrichment(adata, cluster_key="Cell_class")

# Plot heatmap of enrichment scores
sq.pl.nhood_enrichment(
    adata, 
    cluster_key="Cell_class", 
    method="single", 
    cmap="inferno", 
    vmin=-50, 
    vmax=100,
    show=False
)
plt.savefig("merfish_3d_nhood_enrichment.png", bbox_inches="tight", dpi=300)
plt.show()
```

Highlight specific co-enriched cell clusters in 3D.
```python
# Plot target cell classes and make others transparent
sc.pl.embedding(
    adata,
    basis="spatial3d",
    groups=["OD Mature 1", "OD Mature 2", "OD Mature 4"],
    na_color=(1, 1, 1, 0),
    projection="3d",
    color="Cell_class",
    show=False
)
plt.savefig("merfish_3d_coenriched_clusters.png", bbox_inches="tight", dpi=300)
plt.show()
```

Find marker genes and show expression in 3D.
```python
# Rank genes by cell class
sc.tl.rank_genes_groups(adata, groupby="Cell_class", method="t-test")

# Plot Gad1 and Mlc1 expression on 3D scatter
sc.pl.embedding(
    adata, 
    basis="spatial3d", 
    projection="3d", 
    color=["Gad1", "Mlc1"], 
    show=False
)
plt.savefig("merfish_3d_gene_expression.png", bbox_inches="tight", dpi=300)
plt.show()
```

Find spatially variable genes on a 2D slice using Moran's I.
```python
# Copy single slice data
adata_slice = adata[adata.obs.Bregma == -9].copy()

# Build 2D spatial graph
sq.gr.spatial_neighbors(adata_slice, coord_type="generic")

# Calculate spatial autocorrelation
sq.gr.spatial_autocorr(adata_slice, mode="moran")

# Plot top spatially variable genes
sq.pl.spatial_scatter(
    adata_slice, 
    shape=None, 
    color=["Cd24a", "Necab1", "Mlc1"], 
    size=3
)
plt.savefig("merfish_2d_slice_morans_I_genes.png", bbox_inches="tight", dpi=300)
plt.show()
```

## Results

| Output Plot | Description |
| :--- | :--- |
| `merfish_3d_cell_class.png` | 3D scatter plot of all cells colored by class. |
| `merfish_2d_slice_bregma_minus9.png` | 2D scatter plot of single slice (Bregma -9). |
| `merfish_3d_nhood_enrichment.png` | Heatmap showing which cell classes neighbor each other. |
| `merfish_3d_coenriched_clusters.png` | 3D scatter highlighting only OD Mature 1, 2, and 4. |
| `merfish_3d_gene_expression.png` | 3D scatter plots colored by Gad1 and Mlc1 expression. |
| `merfish_2d_slice_morans_I_genes.png` | 2D maps of top spatially autocorrelated genes. |

## Reference
Tutorial source: [Analyze Merfish data — Squidpy](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_merfish.html)
