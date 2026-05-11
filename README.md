# Harmony Batch Correction for Single-cell and Spatial Transcriptomics

## Project Overview

This project performs batch correction and benchmarking on single-cell RNA-seq and spatial transcriptomics melanoma datasets using Harmony. The main aim is to check whether Harmony can reduce batch effects while preserving biological clustering patterns.

The analysis was performed on both:

1. Single-cell RNA-seq data
2. Spatial transcriptomics Visium data

The workflow includes data loading, quality control, preprocessing, normalization, PCA, UMAP visualization before and after Harmony, clustering, cell type annotation using SingleR, and evaluation using multiple metrics.

---

## Project Title

**Benchmarking Harmony Batch Correction on Single-cell and Spatial Transcriptomics Data**

---
## Dataset

The dataset used in this project is from:

**GSE207592**

### Single-cell RNA-seq Samples

The scRNA-seq analysis was performed using multiple melanoma single-cell samples.

Example samples used:

- GSM6300689
- GSM6300693
- GSM6300697
- GSM6300698
- GSM6300699

### Spatial Transcriptomics Samples

The spatial transcriptomics analysis was performed using Visium melanoma samples:

- GSM6319696 - Visium melanoma_2
- GSM6319697 - Visium melanoma_3
- GSM6319698 - Visium melanoma_4

---

## Important Note About Dataset

The raw dataset is **not uploaded** in this GitHub repository because the files are large.

To run the analysis, download the raw data from GEO using accession ID:

```text
GSE207592
