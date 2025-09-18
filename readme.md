# Single-Cell Analysis Pipeline (Streamlit + Scanpy)

This project provides a simple interactive pipeline for preprocessing and analyzing single-cell RNA-seq data (`.h5ad` format) using [Scanpy](https://scanpy.readthedocs.io/) and a [Streamlit](https://streamlit.io/) web interface.

---
## Overview

This project is a Streamlit-based interactive pipeline for single-cell RNA-seq analysis. It allows users to upload .h5ad files and perform common preprocessing steps such as filtering, normalization, variable gene selection, PCA, and UMAP visualization. The app also supports Harmony batch integration for correcting batch effects and includes a differential expression analysis module to identify marker genes between groups of cells. By combining the power of Scanpy with a simple web interface, the tool makes single-cell data analysis more accessible and user-friendly.
## Features

- **Upload `.h5ad` files** directly in the browser.  
- **Data summary**: shows dataset shape and metadata.  
- **Gene filtering**: automatically removes ERCC and mitochondrial (MT) genes.  
- **Preprocessing steps**:
  - Cell and gene filtering
  - Normalization & log-transformation
  - Highly variable gene selection
  - Scaling & PCA
  - Optional Harmony batch integration (if a `batch` column exists in `.obs`)
  - Neighbor graph construction & UMAP embedding
- **Differential expression analysis**:
  - Choose grouping column (`.obs` metadata)
  - Run DE analysis with the Wilcoxon test
  - Results stored inside the AnnData object

---
