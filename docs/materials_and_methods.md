Materials and Methods


Data acquisition and preprocessing
Single-cell or single-nucleus RNA-sequencing datasets in .h5ad format were obtained from publicly available CELLxGENE resources (https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289?explainNewTabFat). The analysis workflow was designed to evaluate age-associated expression patterns of PDCD1LG2 and selected inflammation-, senescence-, and extracellular matrix-associated genes across human tissue datasets. The primary dataset analyzed in the current workflow was the Healthy pediatric and adult human liver tissue dataset. 
Data were imported into R as Seurat objects using the anndataR and Seurat packages. The RNA assay was set as the default assay. When multiple expression layers were present, the X layer was used preferentially; otherwise, the first available RNA layer was selected for downstream analysis. Cell metadata were extracted from the Seurat object and merged with gene expression values using cell barcodes.
Cell type annotations were obtained from the cell_type metadata field. To reduce potential confounding from low-relevance or contaminating cell populations, erythrocytes, red blood cells, platelets, thrombocytes, dendritic cells, conventional dendritic cells, and plasmacytoid dendritic cells were excluded prior to analysis. Dimensionality reduction plots were generated using the available UMAP embedding to inspect cell type distributions after filtering.


Gene selection and annotation
The primary gene of interest was PDCD1LG2. Additional genes included in the analysis were MMP2, MMP9, CDKN1A, CDKN2A, IFNG, and TNF. These genes were selected to evaluate whether PDCD1LG2-positive cells showed increased expression of genes associated with tissue remodeling, cellular senescence, and inflammatory signaling.
Because feature annotation differed across datasets, gene identifiers were evaluated before expression extraction. If the majority of feature names corresponded to Ensembl gene identifiers, Ensembl IDs were used for gene matching. Otherwise, HGNC gene symbols were used. The following gene symbol-to-Ensembl mappings were applied: PDCD1LG2: ENSG00000197646; MMP2: ENSG00000087245; MMP9: ENSG00000100985; CDKN1A: ENSG00000124762; CDKN2A: ENSG00000147889; IFNG: ENSG00000111537; and TNF: ENSG00000232810. Genes not detected in the expression matrix were excluded from the corresponding analysis.


Metadata harmonization and age extraction
Cell-level expression data were merged with donor and sample metadata. Required metadata fields included donor identifier, developmental stage or age annotation, sex, tissue annotation, and cell type. When a tissue column was not available, the tissue_type field was used as the tissue annotation.
Numeric age was extracted from the development_stage metadata field using regular expression matching. Cells without valid numeric age information were excluded from age-associated analyses. All downstream age analyses were performed using donor-level summaries to reduce the effect of unequal cell numbers across donors.


Donor-level analysis of gene expression and age
For each target gene, expression values were summarized at the donor level across all retained cells. Two donor-level metrics were calculated: mean target gene expression and the percentage of target gene-positive cells. Cells were classified as target gene-positive when the normalized expression value of the target gene was greater than zero.
Linear regression models were used to evaluate the association between donor age and target gene expression. Separate models were fitted for mean expression and percentage of positive cells. When available and sufficiently variable across donors, sex and tissue were included as covariates. The general model structures were:
Mean target gene expression ~ age + sex + tissue
Percentage of target gene-positive cells ~ age + sex + tissue
For each model, the coefficient of determination and age-associated P value were extracted. Scatterplots were generated with donor-level values, linear regression fits, confidence intervals, and statistical annotations. Figures were exported as publication-quality PDF files and high-resolution PNG files.


Cell type-specific age association analysis
To determine whether age-associated expression patterns differed by cell type, donor-level summaries were calculated separately within each annotated cell type. For each donor-cell type pair, mean target gene expression, percentage of target gene-positive cells, and the number of contributing cells were calculated. Donor-cell type groups with fewer than 20 cells were excluded.
For each cell type represented by at least four donors, linear regression models were fitted to test the association between age and either mean target gene expression or the percentage of target gene-positive cells:
Mean target gene expression within cell type ~ age
Percentage of target gene-positive cells within cell type ~ age
P values were adjusted for multiple testing using the Benjamini-Hochberg procedure. Cell type-specific trends were classified as increased with age, decreased with age, or not significantly changed based on the adjusted P value and direction of the age coefficient. Cell type-specific summary tables and regression plots were generated for downstream interpretation.


Identification and characterization of target gene-positive cells
To evaluate whether target gene-positive cells expressed higher levels of selected comparison genes, cells were stratified into target gene-positive and target gene-negative populations. A cell was defined as target gene-positive if the expression value of the target gene was greater than zero and target gene-negative if expression was equal to zero.
For each target gene, the remaining genes in the selected gene panel were analyzed as comparison genes. For example, when PDCD1LG2 was treated as the target gene, expression of MMP2, MMP9, CDKN1A, CDKN2A, IFNG, and TNF was compared between PDCD1LG2-positive and PDCD1LG2-negative cells.
Cell-level boxplots were generated to visualize expression differences between target gene-positive and target gene-negative cells. These plots were used for visualization only, because cell-level statistical tests can be biased by large and unequal cell numbers across donors.


Donor-level comparison of target gene-positive and target gene-negative populations
Statistical testing of comparison gene expression between target gene-positive and target gene-negative populations was performed at the donor level. For each donor and target gene status group, mean expression of each comparison gene was calculated. Donor-status groups represented by fewer than 20 cells were excluded.
For each comparison gene, donor-level mean expression was compared between target gene-positive and target gene-negative groups using Wilcoxon rank-sum tests. P values were adjusted using the Benjamini-Hochberg method. For each gene, the analysis reported mean expression, median expression, log2 fold-change between positive and negative populations, raw P value, adjusted P value, and direction of change.
Donor-level boxplots were generated for each comparison gene, with individual points representing donor-level means. Adjusted P values were displayed on the plots.


Cell type-specific comparison of target gene-positive and target gene-negative populations
The positive versus negative comparison was repeated within each cell type. For each donor, cell type, and target gene status group, mean expression of each comparison gene was calculated. Groups with fewer than 10 cells were excluded. Cell type-gene combinations were retained only when both target gene-positive and target gene-negative populations were present and at least four donors were available.
Within each cell type, Wilcoxon rank-sum tests were used to compare donor-level mean expression between target gene-positive and target gene-negative populations. P values were adjusted using the Benjamini-Hochberg method within each cell type. Results were summarized in cell type-specific tables and visualized as publication-ready boxplots, with each point representing a donor-level mean.


Statistical analysis and visualization
All analyses were performed in R using the Seurat, anndataR, dplyr, tidyr, stringr, broom, ggplot2, ggpubr, scales, and fs packages. Linear regression was used to assess associations between age and gene expression metrics. Wilcoxon rank-sum tests were used for donor-level comparisons between target gene-positive and target gene-negative groups. Multiple testing correction was performed using the Benjamini-Hochberg method. Adjusted P values less than 0.05 were considered statistically significant.
Figures were generated using ggplot2 with publication-oriented formatting. Outputs included regression plots, donor-level boxplots, cell type-specific plots, and CSV summary tables containing model statistics and comparison results.

