# Visium H&E Spatial Graph Analysis with Squidpy

Statistical modeling of spatial organization and cell-cell co-occurrence in mouse brain H&E sections. [Source Data](https://support.10xgenomics.com/spatial-gene-expression/datasets).

## Setup
```python
import scanpy as sc
import squidpy as sq
import numpy as np
```

## Steps

Build a spatial graph using coordinates to define neighbors for every spot.
```python
# Compute adjacency matrix based on spatial proximity
sq.gr.spatial_neighbors(adata)
```

Run a permutation test to see which clusters are statistically likely to be adjacent.
```python
# Calculate neighborhood enrichment scores
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
```

Analyze the conditional probability of finding specific clusters at increasing distances.
```python
# Calculate co-occurrence across distance bins
sq.gr.co_occurrence(adata, cluster_key="cluster")
```

Identify spatially variable genes using the Moran's I global autocorrelation statistic.
```python
# Find genes with non-random spatial patterns
sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100)
```

## Results

| Image | Description |
| :--- | :--- |
| **1.spatial_clusters.jpeg** | H&E section with annotated anatomical clusters. |
| **2.neighborhood_enrichment.jpeg** | Heatmap proving the Hippocampus and Pyramidal layers are significant neighbors. |
| **3.co_occurrence_hippocampus.jpeg** | Plot showing how cluster proximity decays over distance. |
