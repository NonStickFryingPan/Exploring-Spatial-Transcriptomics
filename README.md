<div align="center">

# 🔭 Spatial Transcriptomics
**An exploration of 'Spatial Biology' using 10x Genomics Visium & Xenium data.**

[![Scanpy](https://img.shields.io/badge/Analysis-Scanpy-1f6feb?style=flat-square)](https://scanpy.readthedocs.io/en/latest/)
[![Squidpy](https://img.shields.io/badge/Spatial-Squidpy-16a34a?style=flat-square)](https://squidpy.readthedocs.io/)
[![SpatialData](https://img.shields.io/badge/Format-SpatialData%2FZarr-7c3aed?style=flat-square)](https://spatialdata.scverse.org/en/latest/)
[![10x Visium](https://img.shields.io/badge/Platform-10x%20Visium-f59e0b?style=flat-square)](https://www.10xgenomics.com/spatial-gene-expression)
[![10x Xenium](https://img.shields.io/badge/Platform-10x%20Xenium-e11d48?style=flat-square)](https://www.10xgenomics.com/platforms/xenium)

[Overview](#overview) · [Sections](#sections) · [Results](#results) · [Definitions](#definitions)

</div>

---

## Overview

This repo aims to explore spatial transcriptomics across two 10x Genomics scales: Visium (55 µm spots) and Xenium (subcellular single-cell resolution) using [Scanpy](https://scanpy.readthedocs.io/en/latest/) for clustering, [Squidpy](https://squidpy.readthedocs.io/) for image features and spatial statistics (Moran's I, ligand–receptor, neighborhoods), and the [SpatialData/Zarr](https://spatialdata.scverse.org/en/latest/) ecosystem for cell-level analysis.

---

## What is Spatial Transcriptomics?

For decades we've ground up tissues, sequenced them in bulk, or dissociated them into single cells, throwing away the one piece of information that often defines function: *location*. But cells are not isolated beings; they work together in tandem. Spatial biology is a rejection of the idea that a cell's identity can be separated from its location. *It matters where a cell is located.*

> This repo works with 10x Genomics Visium & Xenium data.

---

## Sections
Check each folder link below for its own respective guide on how to replicate it.

| # | Section | Dataset | Focus |
|---|----------|---------|-------|
| 1 | [Basic Scanpy Spatial](basicScanpy/) | Human Lymph Node (Visium) | QC, clustering, spatial visualization, marker genes |
| 2 | [Visium Fluorescence](visFluorescence/) | Mouse Brain (Visium + fluorescence) | Image segmentation, segmentation features |
| 3 | [Visium H&E](visHnE) | Mouse Brain (Visium + H&E) | Neighborhood enrichment, co-occurrence, ligand-receptor |
| 4 | [Xenium](xenium/) | Human Lung Cancer (Xenium) | Single-cell spatial, centrality, co-occurrence, Moran's I |

## Tutorials Followed

**1: [Basic Scanpy Spatial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)**: Standard scRNA-seq pipeline applied to spatial spots. QC → normalize → HVGs → PCA → UMAP → Leiden clusters. The spatial twist: `pl.spatial()` projects clusters back onto H&E tissue. First time seeing *where* clusters live.

**2: [Visium Fluorescence](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html)**: Moved from transcriptomics into image space. Segmented nuclei from DAPI via Watershed, extracted per-spot features (cell count, marker intensity, texture). Ran PCA + Leiden on image features — no sequencing data used. Key finding: image features sub-divided hippocampus into structural layers that gene expression couldn't distinguish.

**3: [Visium H&E](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)**: Shifted from *description* to *proof*. Multi-scale morphology → neighborhood enrichment (permutation test) → co-occurrence curves → spatially-filtered ligand-receptor (ligrec/OmniPath) → Moran's I for spatially variable genes. Linked morphology + proximity + function — the "holy trinity."

**4: [Xenium](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html)**: Single-molecule FISH platform, true single-cell resolution (not spots), transcripts localized within cell boundaries.

---

## Results

### 1: Basic Scanpy

~4,000 spots processed. The key takeaway was the high correlation between transcriptional clusters and the physical H&E image. Cluster 9 was identified as a distinct follicular region, confirmed by high expression of the marker gene *CR2*. These results prove that tissue structures can be digitally identified without a pathologist manually labelling every region.

### 2: Visium Fluorescence

Used DAPI and antibody staining (NEUN/GFAP). Watershed segmentation estimated actual cell counts per spot, which varied significantly across tissue layers. Image-based texture and intensity features provided *finer* detail than genes — gene clusters grouped the hippocampus as one large block, while image features correctly subdivided it into known anatomical sub-layers.

### 3: Visium H&E

**Neighborhood Enrichment** scores showed that the Pyramidal layers and the Hippocampus are statistically significant neighbors, not just transcriptionally similar. **Moran's I** identified genes like *Olfm1* and *Plp1* as highly spatially structured (non-random). Specific ligand-receptor pairs were also found, identifying potential protein-level interactions at cluster borders.

### 4: Xenium

Over 11,000 individual cells analyzed in a lung cancer sample. **Centrality scores** identified which cell types act as "hubs" in the tumor microenvironment. **Delaunay Triangulation** proved that tumor cells were not randomly distributed but formed tight, high-density networks distinct from the surrounding healthy stroma.

---

## Definitions

* **Spots**: Physical circles on the glass slide, 55 µm diameter. Each spot catches RNA from a small group (1–10) of cells. Think of them as "pixels" of gene data.

* **Clusters**: Groups of spots that speak the same "molecular language." If two spots have nearly identical gene expression, the algorithm (Leiden) puts them in the same cluster.

* **AnnData**: The data structure that keeps everything together: the gene counts, the cell names, and the spatial coordinates.

* **H&E Staining**: The "pink and purple" picture. A standard chemical dye: purple (Hematoxylin) shows cell nuclei (DNA), pink (Eosin) shows proteins and structures.

* **PCA & UMAP**: Humans can't visualize 20,000 genes at once. These tools compress thousands of dimensions into a 2D map so we can see which cells are similar.

* **Moran's I**: A score from −1 to 1. A high score means a gene isn't just floating everywhere — it's forming a specific, meaningful pattern or clump in the tissue.**Ligand-Receptor** — A ligand is a signal sent by one cell; a receptor is the ear that hears it. Finding these pairs tells us if two cell types are actually communicating.

* **Watershed**: A segmentation method. Imagine the image as a landscape — the computer "fills the valleys" with water until the "lakes" (cells) touch. Where they touch is the boundary.

* **Delaunay Triangulation**: Creates a network of triangles between cells to determine exactly who is a neighbor to whom.
