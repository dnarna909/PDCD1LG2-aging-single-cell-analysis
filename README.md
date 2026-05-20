
# PDCD1LG2 Aging Single-Cell Analysis

This repository contains R scripts used to analyze age-associated expression of PDCD1LG2 and related genes in human single-cell or single-nucleus RNA-seq datasets.

## Project goals

1. Evaluate whether PDCD1LG2 expression changes with age across all cells.
2. Evaluate whether PDCD1LG2 expression changes with age within individual cell types.
3. Compare expression of selected genes, including MMP2, MMP9, CDKN1A, CDKN2A, IFNG, and TNF, between PDCD1LG2-positive and PDCD1LG2-negative cells.

## Input data

The analysis uses `.h5ad` files downloaded from public CELLxGENE datasets. 

Dataset collection: [Healthy pediatric and adult human liver tissue - CELLxGENE] (https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289?explainNewTab)

Large data files are not included in this repository. Please download the datasets separately and place them in a local data directory.

## Main script

The main analysis script is located at:

```text
scripts/PDCD1LG2_age_analysis.R
```

## How to use this repository

### 1. Clone or download this repository
Run in Terminal:
```bash
git clone https://github.com/dnarna909/PDCD1LG2-aging-single-cell-analysis.git
cd PDCD1LG2-aging-single-cell-analysis
```

Alternatively, download the repository as a ZIP file from GitHub.

### 2. Download the required `.h5ad` file

Download the required `.h5ad` file from CELLxGENE.

### 3. Place the `.h5ad` file in the `data/` folder

The expected file path is:

```text
data/Healthy pediatric and adult human liver tissue.h5ad
```

If your downloaded file has a different name, either rename the file or update the `file_list` section in the R script.

### 4. Install the required R packages

Run the following commands in R or RStudio:

```R
install.packages(c(
  "Seurat",
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

BiocManager::install(c(
  "rhdf5",
  "anndataR"
))
```

### 5. Run the analysis

From R or RStudio:

```r
source("scripts/PDCD1LG2_age_analysis.R")
```

Or from Terminal:

```bash
Rscript scripts/PDCD1LG2_age_analysis.R
```

### 6. View output files

Output files will be saved in the `results/` folder.

The script generates:

Donor-level age association plots
Cell type-specific age association plots
Target gene-positive versus target gene-negative comparison plots
CSV files containing statistical summaries
A `sessionInfo.txt` file for reproducibility
Output structure

Example output structure:

```text
results/
└── Healthy pediatric and adult human liver tissue/
    ├── PDCD1LG2/
    ├── MMP2/
    ├── MMP9/
    ├── CDKN1A/
    ├── CDKN2A/
    ├── IFNG/
    └── TNF/
```

## Notes

The `.h5ad` data files are intentionally excluded from this repository because they are large. Users should download the data directly from CELLxGENE and place the file in the `data/` folder before running the analysis.

## Code availability

The analysis code is available at:

https://github.com/dnarna909/PDCD1LG2-aging-single-cell-analysis
