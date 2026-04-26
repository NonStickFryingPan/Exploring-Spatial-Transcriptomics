# Visium Fluorescence Analysis with Squidpy

Analysis of mouse brain Visium data using fluorescence channels (DAPI/NEUN/GFAP) for image-based feature extraction. [Source Data](https://support.10xgenomics.com/spatial-gene-expression/datasets).

## Setup
```python
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
```

## Steps

Load the fluorescence image and the pre-processed AnnData object.
```python
# Load cropped fluorescence dataset
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()
```

Smooth the image and run watershed segmentation on the DAPI channel to identify nuclei.
```python
# Pre-process and segment
sq.im.process(img=img, layer="image", method="smooth")
sq.im.segment(img=img, layer="image_smooth", method="watershed", channel=0)
```

Calculate image-derived features like cell counts and channel intensity per spot.
```python
# Extract features based on segmentation
sq.im.calculate_image_features(
    adata, img, features="segmentation", 
    layer="image", key_added="features_segmentation"
)
```

Cluster the spots based on image features alone to compare with gene expression clusters.
```python
# Scale, PCA, and cluster image features
sc.pp.scale(features_adata)
sc.pp.pca(features_adata)
sc.tl.leiden(features_adata)
```

## Results

| Image | Description |
| :--- | :--- |
| **1.spatial_clusters.jpeg** | Initial gene-expression clusters overlaid on the tissue. |
| **2.segmentation_comparison.jpeg** | Raw DAPI signal vs. computer-generated nuclei masks. |
| **3.image_features_vs_clusters.jpeg** | Comparison showing image features provide finer structural detail in the Hippocampus. |
