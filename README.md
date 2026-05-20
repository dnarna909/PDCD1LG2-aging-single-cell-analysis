# PDCD1LG2 Aging Single-Cell Analysis

This repository contains R scripts used to analyze age-associated expression of PDCD1LG2 and related genes in human single-cell or single-nucleus RNA-seq datasets.

## Project goals

1. Evaluate whether PDCD1LG2 expression changes with age across all cells.
2. Evaluate whether PDCD1LG2 expression changes with age within individual cell types.
3. Compare expression of selected genes, including MMP2, MMP9, CDKN1A, CDKN2A, IFNG, and TNF, between PDCD1LG2-positive and PDCD1LG2-negative cells.

## Input data

The analysis uses `.h5ad` files downloaded from public CELLxGENE datasets. (https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289?explainNewTabFat)

Large data files are not included in this repository. Please download the datasets separately and place them in a local data directory.

## Main script

The main analysis script is located at:

```text
scripts/PDCD1LG2_age_analysis.R
```text

## How to use this repository

1. Clone or download this repository:

```bash
git clone https://github.com/dnarna909/PDCD1LG2-aging-single-cell-analysis.git
cd PDCD1LG2-aging-single-cell-analysis
```bash

2. Download the required .h5ad file from CELLxGENE.

3. Place the .h5ad file in the data/ folder:
```text
data/Healthy pediatric and adult human liver tissue.h5ad
```text

4. Install the required R packages:
```text
install.packages(c(
  "Seurat",
  "SeuratDisk",
  "dplyr",
  "ggplot2",
  "stringr",
  "fs",
  "tidyr",
  "broom",
  "scales",
  "ggpubr",
  "remotes",
  "BiocManager"
))

remotes::install_github("mojaveazure/seurat-disk")
BiocManager::install(c("rhdf5", "anndataR"))
```text

5. Run the analysis:
```text
source("scripts/PDCD1LG2_age_analysis.R")
```text

Or from Terminal:
```bash
Rscript scripts/PDCD1LG2_age_analysis.R
```bash

6. Output files will be saved in the results/ folder.
```text
One small correction: you wrote **H5D**, but for this script the expected file type is **`.h5ad`**.
::contentReference[oaicite:1]{index=1}
```text
