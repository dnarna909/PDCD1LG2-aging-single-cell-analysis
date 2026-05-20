# ============================================================
# Goal:data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
#   1. Correlation "PDCD1LG2" with age to see if it is positive or negative or no change in all Cells
# 2.Correlation "PDCD1LG2" with age to see if it is positive or negative or no change in all cell types
# 3. The "PDCD1LG2" positive population are also express high of others genes (MMP2, MMP9, CDKN1A, CDKN2A, IFN, TNFA)
# ============================================================

# install if needed -------------------------------------------------
# install.packages("Seurat")
# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")
# BiocManager::install("rhdf5")
# BiocManager::install("anndataR")

# library ---------------------------------------------
library(anndataR)
library(Seurat)
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(fs)
library(tidyr)

# mapping table -----------------------------------------------------------------------
# library(biomaRt)
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# genes_symbol <- c(target_gene_symbol, other_genes_symbol)
# 
# mapping <- getBM(
#   attributes = c("hgnc_symbol", "ensembl_gene_id"),
#   filters = "hgnc_symbol",
#   values = genes_symbol,
#   mart = mart
# )
# mapping

gene_map <- data.frame(
  hgnc_symbol = c("PDCD1LG2", "MMP2", "MMP9", "CDKN1A", "CDKN2A", "IFNG", "TNF"),
  ensembl_gene_id = c(
    "ENSG00000197646", # PDCD1LG2
    "ENSG00000087245", # MMP2
    "ENSG00000100985", # MMP9
    "ENSG00000124762", # CDKN1A
    "ENSG00000147889", # CDKN2A
    "ENSG00000111537", # IFNG
    "ENSG00000232810"  # TNF
  )
)
gene_map

# Define values --------------------------------------------------
# Fat: https://cellxgene.cziscience.com/e/5e5e7a2f-8f1c-42ac-90dc-b4f80f38e84c.cxg/ # not expression
# Tabula Sapiens - Fat.h5ad
# Liver: https://cellxgene.cziscience.com/e/6d41668c-168c-4500-b06a-4674ccf3e19d.cxg/ # not expression
# Tabula Sapiens - Liver.h5ad
# Liver: https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289?explainNewTabFat: # not expression
# Healthy pediatric and adult human liver tissue.h5ad
# Kidney: https://cellxgene.cziscience.com/collections/9c9d04c4-8899-417f-bb6f-6107dcadf14f?explainNewTab
# sn RNA-seq of the Adult Human Kidney.h5ad  

root.dir <- "/media/jianie/ExtremeSSD/2026-04-26 PDCD1LG2 gene for Selim Project/"

# Define genes 2 -----------------------------------------------------------------
# target_gene_symbol <- "PDCD1LG2" #
# other_gene_symbols <- c(
#   "MMP2",
#   "MMP9",
#   "CDKN1A",
#   "CDKN2A",
#   "IFNG",
#   "TNF"
# )  

analysis_gene_symbols <- c(
  "PDCD1LG2",
  "MMP2",
  "MMP9",
  "CDKN1A",
  "CDKN2A",
  "IFNG",
  "TNF"
)

file_list <- list(
  list(names = "Healthy pediatric and adult human liver tissue.h5ad",
       vars = c( "development_stage", "donor_age", "AgeGroup" , "age_condition", "cell_type"  ,"sex" ,"tissue_type"   ))
  # list(names = "sn RNA-seq of the Adult Human Kidney.h5ad" ,
  #      vars = c(""))
)

# file loop -----------------------------------------------------------------------------------------------------
for (ff in seq_along(file_list)) {
  # ff = 3
  # file_name = file_list[[3]][["names"]]
  
  # 1. Load data: Load as a Seurat object directly --------------------------------------------------
  file_name = file_list[[ff]][["names"]]
  out_dir = paste0(root.dir, tools::file_path_sans_ext(file_name))
  dir_create(out_dir)
  
  seurat_obj <- tryCatch(
    {
      read_h5ad(file.path(root.dir, file_name), as = "Seurat")
    },
    error = function(e) {
      warning(paste0("Skipping ", file_name, ": ", e$message))
      return(NULL)
    }
  )
  
  if (is.null(seurat_obj)) {
    next
  }
  DefaultAssay(seurat_obj) <- "RNA"
  Layers(seurat_obj[["RNA"]])
  DimPlot(seurat_obj, reduction = "X_umap", label = TRUE)
  colnames(seurat_obj@meta.data)
  
  
  # Choose the cell type annotation column
  # Change this if your dataset uses another name, e.g. "cell_type", "celltype",
  # "annotation", "predicted_cell_type", "major_cell_type"
  celltype_col <- "cell_type"
  
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0(
      "Cell type column '", celltype_col, "' not found. Available columns are: ",
      paste(colnames(seurat_obj@meta.data), collapse = ", ")
    ))
  }
  
  # Exclude unwanted cell types
  exclude_celltypes <- c(
    "erythrocyte",
    "erythrocytes",
    "red blood cell",
    "red blood cells",
    "RBC",
    "platelet",
    "platelets",
    "thrombocyte",
    "thrombocytes",
    "dendritic cell",
    "dendritic cells",
    "DC",
    "cDC",
    "pDC"
  )
  
  # Case-insensitive filtering
  celltype_values <- seurat_obj@meta.data[[celltype_col]]
  
  keep_cells <- !tolower(celltype_values) %in% tolower(exclude_celltypes)
  
  seurat_obj <- subset(
    seurat_obj,
    cells = colnames(seurat_obj)[keep_cells]
  )
  
  message("Remaining cells after filtering: ", ncol(seurat_obj))
  
  Layers(seurat_obj[["RNA"]])
  
  DimPlot(
    seurat_obj,
    reduction = "X_umap",
    group.by = celltype_col,
    label = TRUE
  )
  
  layer_use <- if ("X" %in% Layers(seurat_obj[["RNA"]])) {
    "X"
  } else {
    Layers(seurat_obj[["RNA"]])[1]
  }
  
  expr_mat <- LayerData(seurat_obj, assay = "RNA", layer = layer_use)
  
  # 2. Check which genes exist ----------------------------------------------------
  feature_names <- rownames(seurat_obj)
  if (mean(grepl("^ENSG", feature_names)) > 0.5) {
    message("This Seurat object likely uses Ensembl gene IDs.")
    target_gene = target_gene_Ensembl
    other_genes = other_gene_Ensembls
    all_genes = all_gene_Ensembls
  } else {
    message("This Seurat object likely uses gene symbols.")
    target_gene = target_gene_symbol
    other_genes = other_gene_symbols
    all_genes = all_gene_symbols
  }
  target_gene
  other_genes
  all_genes
  
  gene_check <- data.frame(
    gene = all_genes,
    present = all_genes %in% rownames(expr_mat)
  )
  print(gene_check)
  
  genes_found <- all_genes[all_genes %in% rownames(expr_mat)]
  if (!target_gene %in% genes_found) {
    stop(paste0(target_gene_symbol," not found in expression matrix. Check gene names or Ensembl IDs."))
  }
  
  
  # 3. Extract expression + metadata ------------------------------------------------
  expr_df <- as.data.frame(
    t(as.matrix(expr_mat[genes_found, colnames(seurat_obj), drop = FALSE]))
  )
  colnames(expr_df)[colnames(expr_df) %in% target_gene] <- "target_gene"
  expr_df$cell_barcode <- rownames(expr_df)
  
  meta_df <- seurat_obj@meta.data
  meta_df$cell_barcode <- rownames(meta_df)
  
  df <- meta_df %>%
    left_join(expr_df, by = "cell_barcode")
  
  if (!"tissue" %in% colnames(df) && "tissue_type" %in% colnames(df)) {
    df$tissue <- df$tissue_type
  }
  
  required_cols <- c("donor_id", "development_stage", "sex", "tissue", "cell_type")
  
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "Skipping ", file_name,
        " because missing metadata columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
    next
  }
  
  
  # 4. Extract numeric age -------------------------------------------
  unique(df$development_stage)
  
  df <- df %>%
    mutate(age = as.numeric(str_extract(development_stage, "\\d+"))) %>%
    filter(!is.na(age))
  summary(df$age)
  
  # loop over each gene as target gene -------------------------------
  for (target_gene_symbol in analysis_gene_symbols) {
    
    message("===================================================")
    message("Running analysis for target gene: ", target_gene_symbol)
    message("===================================================")
    
    other_gene_symbols <- setdiff(analysis_gene_symbols, target_gene_symbol)
    
    target_gene_Ensembl <- gene_map$ensembl_gene_id[
      match(target_gene_symbol, gene_map$hgnc_symbol)
    ]
    
    other_gene_Ensembls <- gene_map$ensembl_gene_id[
      match(other_gene_symbols, gene_map$hgnc_symbol)
    ]
    
    if (mean(grepl("^ENSG", feature_names)) > 0.5) {
      message("This Seurat object uses Ensembl gene IDs.")
      
      target_gene <- target_gene_Ensembl
      other_genes <- other_gene_Ensembls
      all_genes <- c(target_gene, other_genes)
      
      gene_label_map <- setNames(
        other_gene_symbols,
        other_genes
      )
      
    } else {
      message("This Seurat object uses gene symbols.")
      
      target_gene <- target_gene_symbol
      other_genes <- other_gene_symbols
      all_genes <- c(target_gene, other_genes)
      
      gene_label_map <- setNames(
        other_gene_symbols,
        other_genes
      )
    }
    
    gene_check <- data.frame(
      target_gene_symbol = target_gene_symbol,
      gene_symbol = c(target_gene_symbol, other_gene_symbols),
      gene_id_used = all_genes,
      present = all_genes %in% rownames(expr_mat)
    )
    
    print(gene_check)
    
    genes_found <- all_genes[all_genes %in% rownames(expr_mat)]
    
    if (!target_gene %in% genes_found) {
      warning(
        paste0(
          "Skipping target gene ",
          target_gene_symbol,
          " because it was not found in expression matrix."
        )
      )
      next
    }
    
    other_genes_found <- other_genes[other_genes %in% rownames(expr_mat)]
    
    if (length(other_genes_found) == 0) {
      warning(
        paste0(
          "Skipping target gene ",
          target_gene_symbol,
          " because no other genes were found."
        )
      )
      next
    }
    
    # Gene-specific output directory
    gene_out_dir <- paste0(out_dir, "/", target_gene_symbol)
    dir_create(gene_out_dir)
    
    # Extract expression for this target + its other genes
    genes_for_this_analysis <- c(target_gene, other_genes_found)
    
    expr_df <- as.data.frame(
      t(as.matrix(expr_mat[genes_for_this_analysis, colnames(seurat_obj), drop = FALSE]))
    )
    
    # Rename target gene column to generic target_gene
    colnames(expr_df)[colnames(expr_df) == target_gene] <- "target_gene"
    
    expr_df$cell_barcode <- rownames(expr_df)
    
    meta_df <- seurat_obj@meta.data
    meta_df$cell_barcode <- rownames(meta_df)
    
    df <- meta_df %>%
      left_join(expr_df, by = "cell_barcode")
    
    if (!"tissue" %in% colnames(df) && "tissue_type" %in% colnames(df)) {
      df$tissue <- df$tissue_type
    }
    
    df <- df %>%
      mutate(age = as.numeric(str_extract(development_stage, "\\d+"))) %>%
      filter(!is.na(age))
    
    genes_for_positive_analysis <- other_genes_found[other_genes_found %in% colnames(df)]
    
    if (length(genes_for_positive_analysis) == 0) {
      warning(
        paste0(
          "Skipping positive-vs-negative analysis for ",
          target_gene_symbol,
          " because none of the other genes were found in df."
        )
      )
    }
    
    # 5. Goal 1: target_gene vs age in all cells, donor-level ---------------------------------------
    library(dplyr)
    library(ggplot2)
    library(broom)
    library(scales)
    
    donor_expr <- df %>%
      group_by(donor_id, age, sex, tissue) %>%
      summarise(
        mean_target_gene = mean(.data[["target_gene"]], na.rm = TRUE),
        pct_target_gene_pos = mean(.data[["target_gene"]] > 0, na.rm = TRUE) * 100,
        n_cells = n(),
        .groups = "drop"
      )
    
    covars <- c("sex", "tissue")
    
    valid_covars <- covars[
      sapply(covars, function(x) {
        x %in% colnames(donor_expr) &&
          n_distinct(na.omit(donor_expr[[x]])) >= 2
      })
    ]
    
    formula_mean <- as.formula(
      paste(
        "mean_target_gene ~ age",
        paste(valid_covars, collapse = " + "),
        sep = ifelse(length(valid_covars) > 0, " + ", "")
      )
    )
    
    formula_pct <- as.formula(
      paste(
        "pct_target_gene_pos ~ age",
        paste(valid_covars, collapse = " + "),
        sep = ifelse(length(valid_covars) > 0, " + ", "")
      )
    )
    
    model_all_mean <- lm(formula_mean, data = donor_expr)
    model_all_pct  <- lm(formula_pct, data = donor_expr)
    
    # Extract R2 and P value for age
    get_model_label <- function(model) {
      model_summary <- summary(model)
      
      r2 <- model_summary$r.squared
      
      age_p <- broom::tidy(model) %>%
        filter(term == "age") %>%
        pull(p.value)
      
      p_label <- ifelse(
        age_p < 0.001,
        "P < 0.001",
        paste0("P = ", signif(age_p, 3))
      )
      
      paste0(
        "R² = ", round(r2, 3),
        "\n", p_label
      )
    }
    
    label_mean <- get_model_label(model_all_mean)
    label_pct  <- get_model_label(model_all_pct)
    
    # Publication-ready theme
    # Better annotation function
    add_stats_label <- function(label) {
      annotate(
        "label",
        x = -Inf,
        y = Inf,
        label = label,
        hjust = -0.05,
        vjust = 3.5,
        size = 4,
        fontface = "bold",
        label.size = NA,
        fill = "white",
        alpha = 0.8
      )
    }
    
    theme_pub <- theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        axis.line = element_line(linewidth = 0.6, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        plot.margin = margin(12, 20, 12, 20)
      )
    
    # Mean expression plot
    p_all_mean <- ggplot(donor_expr, aes(x = age, y = mean_target_gene)) +
      geom_point(
        aes(size = n_cells),
        alpha = 0.75,
        color = "black"
      ) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#2C7BB6",
        fill = "#ABD9E9",
        linewidth = 1
      ) +
      add_stats_label(label_mean) +
      scale_size_continuous(
        name = "Cells per donor",
        range = c(2.5, 7)
      ) +
      labs(
        title = paste0(target_gene_symbol, " mean expression vs age"),
        x = "Age",
        y = paste0("Mean ", target_gene_symbol, " expression per donor")
      ) +
      theme_pub
    
    print(p_all_mean)
    
    # Percent-positive plot
    p_all_pct <- ggplot(donor_expr, aes(x = age, y = pct_target_gene_pos)) +
      geom_point(
        aes(size = n_cells),
        alpha = 0.75,
        color = "black"
      ) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#D7191C",
        fill = "#FDAE61",
        linewidth = 1
      ) +
      add_stats_label(label_pct) +
      scale_size_continuous(
        name = "Cells per donor",
        range = c(2.5, 7)
      ) +
      labs(
        title = paste0(target_gene_symbol, "-positive cell percentage vs age"),
        x = "Age",
        y = paste0("% ", target_gene_symbol, "+ cells per donor")
      ) +
      theme_pub
    
    print(p_all_pct)
    
    # Optional: inspect model summaries
    summary(model_all_mean)
    summary(model_all_pct)
    
    # Save plots: vector PDF + high-resolution PNG
    ggsave(
      filename = paste0(gene_out_dir, "/", target_gene_symbol, "_vs_age_all_donors_mean.pdf"),
      plot = p_all_mean,
      width = 6,
      height = 5,
      device = cairo_pdf
    )
    
    ggsave(
      filename = paste0(gene_out_dir, "/", target_gene_symbol, "_vs_age_all_donors_percent_positive.pdf"),
      plot = p_all_pct,
      width = 6,
      height = 5,
      device = cairo_pdf
    )
    
    ggsave(
      filename = paste0(gene_out_dir, "/", target_gene_symbol, "_vs_age_all_donors_mean.png"),
      plot = p_all_mean,
      width = 6,
      height = 5,
      dpi = 600
    )
    
    ggsave(
      filename = paste0(gene_out_dir, "/", target_gene_symbol, "_vs_age_all_donors_percent_positive.png"),
      plot = p_all_pct,
      width = 6,
      height = 5,
      dpi = 600
    )
    
    # 6. Goal 2: target_gene vs age by cell type -------------------------------------------------------
    library(dplyr)
    library(ggplot2)
    library(broom)
    library(stringr)
    library(fs)
    
    celltype_expr <- df %>%
      group_by(donor_id, age, sex, tissue, cell_type) %>%
      summarise(
        mean_target_gene = mean(.data[["target_gene"]], na.rm = TRUE),
        pct_target_gene_pos = mean(.data[["target_gene"]] > 0, na.rm = TRUE) * 100,
        n_cells = n(),
        .groups = "drop"
      ) %>%
      filter(n_cells >= 20)
    
    # Create output folders
    celltype_plot_dir <- paste0(gene_out_dir, "/", target_gene_symbol, "_celltype_plots")
    dir_create(celltype_plot_dir)
    
    mean_plot_dir <- paste0(celltype_plot_dir, "/mean_expression")
    pct_plot_dir  <- paste0(celltype_plot_dir, "/percent_positive")
    
    dir_create(mean_plot_dir)
    dir_create(pct_plot_dir)
    
    
    # Publication-ready theme
    theme_pub <- theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold"),
        axis.line = element_line(linewidth = 0.6, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        plot.margin = margin(12, 16, 12, 16)
      )
    
    # Helper: safe file names
    safe_filename <- function(x) {
      x %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("_+", "_") %>%
        str_replace_all("^_|_$", "")
    }
    
    # Helper: format p values
    format_p <- function(p) {
      case_when(
        is.na(p) ~ "P = NA",
        p < 0.001 ~ "P < 0.001",
        TRUE ~ paste0("P = ", signif(p, 3))
      )
    }
    
    # Helper: get R2 and P for age
    get_lm_label <- function(data, y_var) {
      model <- lm(as.formula(paste0(y_var, " ~ age")), data = data)
      
      r2 <- summary(model)$r.squared
      
      p_age <- broom::tidy(model) %>%
        filter(term == "age") %>%
        pull(p.value)
      
      paste0(
        "R² = ", round(r2, 3),
        "\n", format_p(p_age)
      )
    }
    
    # Helper: add label to upper-left corner
    add_stats_label <- function(label) {
      annotate(
        "label",
        x = -Inf,
        y = Inf,
        label = label,
        hjust = -0.05,
        vjust = 1.1,
        size = 4,
        fontface = "bold",
        label.size = NA,
        fill = "white",
        alpha = 0.85
      )
    }
    
    # Helper: make one publication-ready plot
    make_celltype_plot <- function(data, celltype_name, y_var, y_label, title_text, line_color, fill_color) {
      label_text <- get_lm_label(data, y_var)
      
      ggplot(data, aes(x = age, y = .data[[y_var]])) +
        geom_point(
          aes(size = n_cells),
          alpha = 0.75,
          color = "black"
        ) +
        geom_smooth(
          method = "lm",
          se = TRUE,
          linewidth = 1,
          color = line_color,
          fill = fill_color
        ) +
        add_stats_label(label_text) +
        scale_size_continuous(
          name = "Cells per donor",
          range = c(2.5, 7)
        ) +
        labs(
          title = title_text,
          x = "Age",
          y = y_label
        ) +
        theme_pub
    }
    
    ## Fit statistics by cell type -------------------------------------------------
    
    celltype_stats <- celltype_expr %>%
      group_by(cell_type) %>%
      filter(n_distinct(donor_id) >= 4) %>%
      summarise(
        n_donors = n_distinct(donor_id),
        
        beta_age_mean = coef(lm(mean_target_gene ~ age))[["age"]],
        r2_mean = summary(lm(mean_target_gene ~ age))$r.squared,
        p_age_mean = summary(lm(mean_target_gene ~ age))$coefficients["age", "Pr(>|t|)"],
        
        beta_age_pct = coef(lm(pct_target_gene_pos ~ age))[["age"]],
        r2_pct = summary(lm(pct_target_gene_pos ~ age))$r.squared,
        p_age_pct = summary(lm(pct_target_gene_pos ~ age))$coefficients["age", "Pr(>|t|)"],
        
        .groups = "drop"
      ) %>%
      mutate(
        p_adj_mean = p.adjust(p_age_mean, method = "BH"),
        p_adj_pct = p.adjust(p_age_pct, method = "BH"),
        
        direction_mean = case_when(
          p_adj_mean < 0.05 & beta_age_mean > 0 ~ "Increase with age",
          p_adj_mean < 0.05 & beta_age_mean < 0 ~ "Decrease with age",
          TRUE ~ "No significant change"
        ),
        
        direction_pct = case_when(
          p_adj_pct < 0.05 & beta_age_pct > 0 ~ "Increase with age",
          p_adj_pct < 0.05 & beta_age_pct < 0 ~ "Decrease with age",
          TRUE ~ "No significant change"
        ),
        
        file_name = file_name
      ) %>%
      arrange(p_adj_mean)
    
    print(celltype_stats)
    
    write.csv(
      celltype_stats,
      paste0(gene_out_dir, "/", target_gene_symbol, "_age_correlation_by_celltype.csv"),
      row.names = FALSE
    )
    
    ## Save one figure per cell type ----------------------------------------------
    
    celltypes_to_plot <- celltype_stats %>%
      filter(n_donors >= 4) %>%
      pull(cell_type)
    
    for (ct in celltypes_to_plot) {
      
      plot_data <- celltype_expr %>%
        filter(cell_type == ct)
      
      ct_safe <- safe_filename(ct)
      
      p_mean <- make_celltype_plot(
        data = plot_data,
        celltype_name = ct,
        y_var = "mean_target_gene",
        y_label = paste0("Mean ", target_gene_symbol, " expression per donor"),
        title_text = paste0(target_gene_symbol, " mean expression vs age\n", ct),
        line_color = "#2C7BB6",
        fill_color = "#ABD9E9"
      )
      
      p_pct <- make_celltype_plot(
        data = plot_data,
        celltype_name = ct,
        y_var = "pct_target_gene_pos",
        y_label = paste0("% ", target_gene_symbol, "+ cells per donor"),
        title_text = paste0(target_gene_symbol, "-positive cell percentage vs age\n", ct),
        line_color = "#D7191C",
        fill_color = "#FDAE61"
      )
      
      print(p_mean)
      print(p_pct)
      
      ggsave(
        filename = paste0(mean_plot_dir, "/", target_gene_symbol, "_", ct_safe, "_mean_expression_vs_age.pdf"),
        plot = p_mean,
        width = 6.5,
        height = 5,
        device = cairo_pdf
      )
      
      ggsave(
        filename = paste0(mean_plot_dir, "/", target_gene_symbol, "_", ct_safe, "_mean_expression_vs_age.png"),
        plot = p_mean,
        width = 6.5,
        height = 5,
        dpi = 600
      )
      
      ggsave(
        filename = paste0(pct_plot_dir, "/", target_gene_symbol, "_", ct_safe, "_percent_positive_vs_age.pdf"),
        plot = p_pct,
        width = 6.5,
        height = 5,
        device = cairo_pdf
      )
      
      ggsave(
        filename = paste0(pct_plot_dir, "/", target_gene_symbol, "_", ct_safe, "_percent_positive_vs_age.png"),
        plot = p_pct,
        width = 6.5,
        height = 5,
        dpi = 600
      )
    }
    
    # 7. Goal 3: Are target_gene+ cells high for other genes? --------------------------------------------
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
    library(fs)
    library(ggpubr)
    
    genes_for_positive_analysis <- other_genes_found[other_genes_found %in% colnames(df)]
    
    gene_label_map <- setNames(
      other_gene_symbols[match(other_genes_found, other_genes)],
      other_genes_found
    )
    
    positive_label <- paste0(target_gene_symbol, "+")
    negative_label <- paste0(target_gene_symbol, "-")
    
    df_pos <- df %>%
      mutate(
        target_gene_status = ifelse(
          .data[["target_gene"]] > 0,
          positive_label,
          negative_label
        ),
        target_gene_status = factor(
          target_gene_status,
          levels = c(negative_label, positive_label)
        )
      )
    
    # Output folder
    positive_plot_dir <- paste0(gene_out_dir, "/", target_gene_symbol, "_positive_vs_negative_plots")
    dir_create(positive_plot_dir)
    
    ## Publication-ready theme ----------------------------------------------------------------------
    theme_pub <- theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.6, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        plot.margin = margin(12, 16, 12, 16)
      )
    
    format_p <- function(p) {
      case_when(
        is.na(p) ~ "adj. P = NA",
        p < 0.001 ~ "adj. P < 0.001",
        TRUE ~ paste0("adj. P = ", signif(p, 3))
      )
    }
    ## Cell-level visualization ----------------------------------------------------------
    plot_df <- df_pos %>%
      dplyr::select(
        cell_barcode,
        donor_id,
        age,
        cell_type,
        target_gene_status,
        all_of(genes_for_positive_analysis)
      ) %>%
      pivot_longer(
        cols = all_of(genes_for_positive_analysis),
        names_to = "gene",
        values_to = "expression"
      ) %>%
      mutate(
        gene_symbol = gene_label_map[gene]
      )
    
    p_pos_box <- ggplot(
      plot_df,
      aes(x = target_gene_status, y = expression, fill = target_gene_status)
    ) +
      geom_boxplot(
        width = 0.6,
        outlier.size = 0.15,
        linewidth = 0.5
      ) +
      facet_wrap(~ gene_symbol, scales = "free_y") +
      scale_fill_manual(
        values = setNames(
          c("grey80", "#D7191C"),
          c(negative_label, positive_label)
        )
      ) +
      labs(
        title = paste0("Selected gene expression in ", target_gene_symbol, "+ vs ", target_gene_symbol, "- cells"),
        subtitle = "Cell-level visualization only",
        x = "",
        y = "Expression"
      ) +
      theme_pub
    
    print(p_pos_box)
    
    ggsave(
      filename = paste0(positive_plot_dir, "/", target_gene_symbol, "_positive_vs_negative_other_genes_cell_level_boxplot.pdf"),
      plot = p_pos_box,
      width = 10,
      height = 6,
      device = cairo_pdf
    )
    
    ggsave(
      filename = paste0(positive_plot_dir, "/", target_gene_symbol, "_positive_vs_negative_other_genes_cell_level_boxplot.png"),
      plot = p_pos_box,
      width = 10,
      height = 6,
      dpi = 600
    )
    
    ## Donor-level test: safer than raw cell-level test -----------------------------
    donor_pos_expr <- df_pos %>%
      group_by(donor_id, age, sex, tissue, target_gene_status) %>%
      summarise(
        across(
          all_of(genes_for_positive_analysis),
          ~ mean(.x, na.rm = TRUE),
          .names = "mean_{.col}"
        ),
        n_cells = n(),
        .groups = "drop"
      ) %>%
      filter(n_cells >= 20)
    
    donor_pos_long <- donor_pos_expr %>%
      pivot_longer(
        cols = starts_with("mean_"),
        names_to = "gene",
        values_to = "mean_expression"
      ) %>%
      mutate(
        gene = str_remove(gene, "^mean_"),
        gene_symbol = gene_label_map[gene],
        target_gene_status = factor(
          target_gene_status,
          levels = c(negative_label, positive_label)
        )
      )
    
    positive_stats <- donor_pos_long %>%
      group_by(gene, gene_symbol) %>%
      summarise(
        n_negative = sum(target_gene_status == negative_label),
        n_positive = sum(target_gene_status == positive_label),
        
        mean_positive = mean(
          mean_expression[target_gene_status == positive_label],
          na.rm = TRUE
        ),
        mean_negative = mean(
          mean_expression[target_gene_status == negative_label],
          na.rm = TRUE
        ),
        
        median_positive = median(
          mean_expression[target_gene_status == positive_label],
          na.rm = TRUE
        ),
        median_negative = median(
          mean_expression[target_gene_status == negative_label],
          na.rm = TRUE
        ),
        
        log2FC_positive_vs_negative = log2((mean_positive + 0.01) / (mean_negative + 0.01)),
        
        p_value = ifelse(
          n_distinct(target_gene_status) == 2 &&
            n_positive >= 2 &&
            n_negative >= 2,
          wilcox.test(mean_expression ~ target_gene_status, data = cur_data())$p.value,
          NA_real_
        ),
        
        .groups = "drop"
      ) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        p_label = format_p(p_adj),
        direction = case_when(
          p_adj < 0.05 & log2FC_positive_vs_negative > 0 ~ paste0("Higher in ", target_gene_symbol, "+ cells"),
          p_adj < 0.05 & log2FC_positive_vs_negative < 0 ~ paste0("Lower in ", target_gene_symbol, "+ cells"),
          TRUE ~ "No significant difference"
        )
      ) %>%
      arrange(p_adj)
    
    print(positive_stats)
    
    write.csv(
      positive_stats,
      paste0(gene_out_dir, "/", target_gene_symbol, "_positive_population_other_genes_stats.csv"),
      row.names = FALSE
    )
    
    ## Donor-level publication-ready figure with adjusted P values ---------------------------
    label_df <- donor_pos_long %>%
      group_by(gene, gene_symbol) %>%
      summarise(
        y_max = max(mean_expression, na.rm = TRUE),
        y_min = min(mean_expression, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        y_range = ifelse(y_max - y_min == 0, y_max, y_max - y_min),
        y.position = y_max + 0.02 * y_range,
        y.expand = y_max + 0.08 * y_range
      ) %>%
      left_join(
        positive_stats %>% dplyr::select(gene, p_label),
        by = "gene"
      ) %>%
      mutate(
        group1 = negative_label,
        group2 = positive_label
      )
    
    p_pos_donor <- ggplot(
      donor_pos_long,
      aes(x = target_gene_status, y = mean_expression, fill = target_gene_status)
    ) +
      # Forces enough vertical room in every free-y facet
      geom_blank(
        data = label_df,
        aes(x = group1, y = y.expand),
        inherit.aes = FALSE
      ) +
      geom_boxplot(
        width = 0.55,
        outlier.shape = NA,
        linewidth = 0.6,
        alpha = 0.85
      ) +
      geom_jitter(
        aes(size = n_cells),
        width = 0.12,
        alpha = 0.7,
        color = "black"
      ) +
      ggpubr::stat_pvalue_manual(
        label_df,
        label = "p_label",
        y.position = "y.position",
        tip.length = 0.01,
        size = 3.5,
        fontface = "bold",
        bracket.size = 0.4
      ) +
      facet_wrap(~ gene_symbol, scales = "free_y") +
      scale_fill_manual(
        values = setNames(
          c("grey80", "#D7191C"),
          c(negative_label, positive_label)
        )
      ) +
      scale_size_continuous(
        name = "Cells per donor",
        range = c(1.8, 4.8)
      ) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0("Donor-level mean expression in ", target_gene_symbol, "+ vs ", target_gene_symbol, "- cells"),
        subtitle = "Each point represents one donor-level mean; P values are BH-adjusted Wilcoxon tests",
        x = "",
        y = "Mean expression per donor"
      ) +
      theme_pub +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(16, 24, 16, 24),
        panel.spacing = unit(1.2, "lines")
      )
    
    print(p_pos_donor)
    
    ggsave(
      filename = paste0(
        positive_plot_dir,
        "/",
        target_gene_symbol,
        "_positive_vs_negative_other_genes_donor_level.pdf"
      ),
      plot = p_pos_donor,
      width = 10.5,
      height = 7.0,
      device = cairo_pdf
    )
    
    ggsave(
      filename = paste0(
        positive_plot_dir,
        "/",
        target_gene_symbol,
        "_positive_vs_negative_other_genes_donor_level.png"
      ),
      plot = p_pos_donor,
      width = 10.5,
      height = 7.0,
      dpi = 600
    )
    
    # 8. Repeat target_gene+ analysis by cell type ------------------------------------------------------
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(stringr)
    library(fs)
    library(ggpubr)
    
    
    positive_label <- paste0(target_gene_symbol, "+")
    negative_label <- paste0(target_gene_symbol, "-")
    
    # Make sure status labels are publication-friendly and consistent
    df_pos <- df_pos %>%
      mutate(
        target_gene_status = ifelse(
          .data[["target_gene"]] > 0,
          positive_label,
          negative_label
        ),
        target_gene_status = factor(
          target_gene_status,
          levels = c(negative_label, positive_label)
        )
      )
    
    # Output folder
    celltype_pos_plot_dir <- paste0(
      gene_out_dir,
      "/",
      target_gene_symbol,
      "_positive_vs_negative_by_celltype_plots"
    )
    
    dir_create(celltype_pos_plot_dir)
    
    # Publication-ready theme
    theme_pub <- theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        axis.line = element_line(linewidth = 0.6, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        plot.margin = margin(14, 22, 14, 22),
        panel.spacing = unit(1.1, "lines")
      )
    
    ## Helper functions ------------------------------------------------------------
    
    safe_filename <- function(x) {
      x %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("_+", "_") %>%
        str_replace_all("^_|_$", "")
    }
    
    format_p <- function(p) {
      case_when(
        is.na(p) ~ "adj. P = NA",
        p < 0.001 ~ "adj. P < 0.001",
        TRUE ~ paste0("adj. P = ", signif(p, 3))
      )
    }
    
    ## Donor-level summarized data -------------------------------------------------
    
    celltype_pos_long <- df_pos %>%
      group_by(donor_id, age, sex, tissue, cell_type, target_gene_status) %>%
      summarise(
        across(
          all_of(genes_for_positive_analysis),
          ~ mean(.x, na.rm = TRUE),
          .names = "mean_{.col}"
        ),
        n_cells = n(),
        .groups = "drop"
      ) %>%
      filter(n_cells >= 10) %>%
      pivot_longer(
        cols = starts_with("mean_"),
        names_to = "gene",
        values_to = "mean_expression"
      ) %>%
      mutate(
        gene = str_remove(gene, "^mean_"),
        gene_symbol = gene_label_map[gene],
        target_gene_status = factor(
          target_gene_status,
          levels = c(negative_label, positive_label)
        )
      )
    
    ## Statistics by cell type and gene -------------------------------------------
    
    celltype_positive_stats <- celltype_pos_long %>%
      group_by(cell_type, gene, gene_symbol) %>%
      filter(
        n_distinct(target_gene_status) == 2,
        n_distinct(donor_id) >= 4
      ) %>%
      summarise(
        n_donors = n_distinct(donor_id),
        n_negative = sum(target_gene_status == negative_label),
        n_positive = sum(target_gene_status == positive_label),
        
        mean_positive = mean(
          mean_expression[target_gene_status == positive_label],
          na.rm = TRUE
        ),
        mean_negative = mean(
          mean_expression[target_gene_status == negative_label],
          na.rm = TRUE
        ),
        
        median_positive = median(
          mean_expression[target_gene_status == positive_label],
          na.rm = TRUE
        ),
        median_negative = median(
          mean_expression[target_gene_status == negative_label],
          na.rm = TRUE
        ),
        
        log2FC_positive_vs_negative = log2(
          (mean_positive + 0.01) / (mean_negative + 0.01)
        ),
        
        p_value = ifelse(
          n_positive >= 2 && n_negative >= 2,
          wilcox.test(mean_expression ~ target_gene_status, data = cur_data())$p.value,
          NA_real_
        ),
        
        .groups = "drop"
      ) %>%
      group_by(cell_type) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        p_label = format_p(p_adj),
        direction = case_when(
          p_adj < 0.05 & log2FC_positive_vs_negative > 0 ~ paste0("Higher in ", target_gene_symbol, "+ cells"),
          p_adj < 0.05 & log2FC_positive_vs_negative < 0 ~ paste0("Lower in ", target_gene_symbol, "+ cells"),
          TRUE ~ "No significant difference"
        )
      ) %>%
      ungroup() %>%
      arrange(cell_type, p_adj)
    
    print(celltype_positive_stats)
    
    write.csv(
      celltype_positive_stats,
      paste0(gene_out_dir, "/", target_gene_symbol, "_positive_population_other_genes_by_celltype.csv"),
      row.names = FALSE
    )
    
    ## one publication-ready figure per cell type:-------------------------------------
    celltypes_to_plot <- celltype_positive_stats %>%
      distinct(cell_type) %>%
      pull(cell_type)
    
    for (ct in celltypes_to_plot) {
      
      plot_data <- celltype_pos_long %>%
        filter(cell_type == ct) %>%
        semi_join(
          celltype_positive_stats %>% filter(cell_type == ct),
          by = c("cell_type", "gene", "gene_symbol")
        )
      
      stat_data <- celltype_positive_stats %>%
        filter(cell_type == ct)
      
      if (nrow(plot_data) == 0 || nrow(stat_data) == 0) {
        next
      }
      
      label_df <- plot_data %>%
        group_by(gene, gene_symbol) %>%
        summarise(
          y_max = max(mean_expression, na.rm = TRUE),
          y_min = min(mean_expression, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          y_range = ifelse(y_max - y_min == 0, y_max, y_max - y_min),
          y_range = ifelse(y_range == 0, 0.1, y_range),
          y.position = y_max + 0.04 * y_range,
          y.expand = y_max + 0.12 * y_range
        ) %>%
        left_join(
          stat_data %>% dplyr::select(gene, p_label),
          by = "gene"
        ) %>%
        mutate(
          group1 = negative_label,
          group2 = positive_label
        )
      
      p_ct <- ggplot(
        plot_data,
        aes(x = target_gene_status, y = mean_expression, fill = target_gene_status)
      ) +
        geom_blank(
          data = label_df,
          aes(x = group1, y = y.expand),
          inherit.aes = FALSE
        ) +
        geom_boxplot(
          width = 0.55,
          outlier.shape = NA,
          linewidth = 0.6,
          alpha = 0.85
        ) +
        geom_jitter(
          aes(size = n_cells),
          width = 0.12,
          alpha = 0.7,
          color = "black"
        ) +
        ggpubr::stat_pvalue_manual(
          label_df,
          label = "p_label",
          y.position = "y.position",
          tip.length = 0.01,
          size = 3.4,
          fontface = "bold",
          bracket.size = 0.4
        ) +
        facet_wrap(~ gene_symbol, scales = "free_y", ncol = 3) +
        scale_fill_manual(
          values = setNames(
            c("grey80", "#D7191C"),
            c(negative_label, positive_label)
          )
        ) +
        scale_size_continuous(
          name = "Cells per donor",
          range = c(1.8, 4.8)
        ) +
        coord_cartesian(clip = "off") +
        labs(
          title = paste0(
            target_gene_symbol,
            "+ vs ",
            target_gene_symbol,
            "- cells in ",
            ct
          ),
          subtitle = "Each point represents one donor-level mean; P values are BH-adjusted Wilcoxon tests",
          x = "",
          y = "Mean expression per donor",
          fill = "Cell status"
        ) +
        theme_pub +
        theme(
          legend.position = "right"
        )
      
      print(p_ct)
      
      ct_safe <- safe_filename(ct)
      
      ggsave(
        filename = paste0(
          celltype_pos_plot_dir,
          "/",
          target_gene_symbol,
          "_positive_vs_negative_other_genes_",
          ct_safe,
          ".pdf"
        ),
        plot = p_ct,
        width = 10.5,
        height = 6.8,
        device = cairo_pdf
      )
      
      ggsave(
        filename = paste0(
          celltype_pos_plot_dir,
          "/",
          target_gene_symbol,
          "_positive_vs_negative_other_genes_",
          ct_safe,
          ".png"
        ),
        plot = p_ct,
        width = 10.5,
        height = 6.8,
        dpi = 600
      )
    }
    
    gc()
    
    gc()
    
  }
}