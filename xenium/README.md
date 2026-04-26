# Xenium Single-Cell Analysis with SpatialData

Subcellular transcript mapping and single-cell graph centrality analysis of human lung cancer. [Data Info](https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-replicate-1-standard).

## Setup
```python
import spatialdata as sd
import squidpy as sq
import scanpy as sc
```

## Steps

Read the Xenium output and convert to Zarr for high-performance handling of millions of transcripts.
```python
# Load Xenium data into SpatialData object
sdata = sd.read_zarr("Xenium.zarr")
adata = sdata.tables["table"]
```

Filter cells by transcript counts and normalize for single-cell clustering.
```python
# Standard single-cell QC and clustering
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.tl.leiden(adata)
```

Build a Delaunay triangulation graph to map the irregular physical neighbors of individual cells.
```python
# Connect dots between cell centers
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```

Calculate centrality metrics to identify "hub" cells in the tumor microenvironment.
```python
# Compute closeness and degree centrality
sq.gr.centrality_scores(adata, cluster_key="leiden")
```

## Results

| Image | Description |
| :--- | :--- |
| **1.spatial_clusters.jpeg** | Thousands of individual cells colored by Leiden cluster at single-cell resolution. |
| **2.centrality_scores.jpeg** | Plot identifying which clusters act as structural hubs in the tumor. |
| **3.co_occurrence.jpeg** | Spatial niche analysis showing tumor cell aggregation patterns. |
