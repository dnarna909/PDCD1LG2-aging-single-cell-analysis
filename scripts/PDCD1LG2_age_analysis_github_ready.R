#!/usr/bin/env Rscript

# ============================================================
# PDCD1LG2 aging single-cell analysis
# ============================================================
#
# Goals:
#   1. Test whether target gene expression changes with age across all cells.
#   2. Test whether target gene expression changes with age within cell types.
#   3. Compare selected genes between target gene-positive and target gene-negative cells.
#
# GitHub-ready notes:
#   - Do not hard-code local computer paths.
#   - Place input .h5ad files in the data/ directory, or set DATA_DIR.
#   - Results are written to the results/ directory, or set RESULTS_DIR.
#
# Example:
#   Rscript scripts/PDCD1LG2_age_analysis_github_ready.R
#
# Optional environment variables:
#   DATA_DIR=/path/to/h5ad/files RESULTS_DIR=/path/to/results Rscript scripts/PDCD1LG2_age_analysis_github_ready.R
# ============================================================


# -----------------------------
# 1. Package setup
# -----------------------------

required_packages <- c(
  "anndataR",
  "Seurat",
  "dplyr",
  "ggplot2",
  "stringr",
  "fs",
  "tidyr",
  "broom",
  "scales",
  "ggpubr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(anndataR)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(fs)
  library(tidyr)
  library(broom)
  library(scales)
  library(ggpubr)
})


# -----------------------------
# 2. User configuration
# -----------------------------

# Input and output directories.
# By default, the script expects this repository structure:
#   data/
#   results/
#
# You can override these with environment variables:
#   DATA_DIR=/path/to/data
#   RESULTS_DIR=/path/to/results

data_dir <- Sys.getenv("DATA_DIR", unset = "data")
results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")

fs::dir_create(data_dir)
fs::dir_create(results_dir)

# Input datasets.
# Add more files here if needed.
file_list <- list(
  list(
    name = "Healthy pediatric and adult human liver tissue.h5ad",
    celltype_col = "cell_type",
    age_col = "development_stage",
    required_metadata = c("donor_id", "development_stage", "sex", "cell_type")
  )
  # Example for adding another dataset:
  # list(
  #   name = "sn RNA-seq of the Adult Human Kidney.h5ad",
  #   celltype_col = "cell_type",
  #   age_col = "development_stage",
  #   required_metadata = c("donor_id", "development_stage", "sex", "cell_type")
  # )
)

# Genes to analyze.
analysis_gene_symbols <- c(
  "PDCD1LG2",
  "MMP2",
  "MMP9",
  "CDKN1A",
  "CDKN2A",
  "IFNG",
  "TNF"
)

# Gene symbol to Ensembl ID mapping.
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
  ),
  stringsAsFactors = FALSE
)

# Cell types excluded before analysis.
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

# Minimum cells and donors used in summary analyses.
min_cells_per_donor_celltype <- 20
min_cells_positive_negative_all <- 20
min_cells_positive_negative_celltype <- 10
min_donors_per_celltype <- 4


# -----------------------------
# 3. Helper functions
# -----------------------------

safe_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

format_p <- function(p, prefix = "P") {
  dplyr::case_when(
    is.na(p) ~ paste0(prefix, " = NA"),
    p < 0.001 ~ paste0(prefix, " < 0.001"),
    TRUE ~ paste0(prefix, " = ", signif(p, 3))
  )
}

detect_gene_id_type <- function(feature_names) {
  if (mean(grepl("^ENSG", feature_names)) > 0.5) {
    "ensembl"
  } else {
    "symbol"
  }
}

get_gene_ids <- function(target_gene_symbol, analysis_gene_symbols, gene_map, id_type) {
  other_gene_symbols <- setdiff(analysis_gene_symbols, target_gene_symbol)

  if (id_type == "ensembl") {
    target_gene <- gene_map$ensembl_gene_id[
      match(target_gene_symbol, gene_map$hgnc_symbol)
    ]

    other_genes <- gene_map$ensembl_gene_id[
      match(other_gene_symbols, gene_map$hgnc_symbol)
    ]

    gene_label_map <- setNames(other_gene_symbols, other_genes)
  } else {
    target_gene <- target_gene_symbol
    other_genes <- other_gene_symbols
    gene_label_map <- setNames(other_gene_symbols, other_genes)
  }

  list(
    target_gene = target_gene,
    other_genes = other_genes,
    all_genes = c(target_gene, other_genes),
    other_gene_symbols = other_gene_symbols,
    gene_label_map = gene_label_map
  )
}

get_expression_matrix <- function(seurat_obj) {
  assay_name <- DefaultAssay(seurat_obj)

  available_layers <- Layers(seurat_obj[[assay_name]])

  layer_use <- if ("X" %in% available_layers) {
    "X"
  } else {
    available_layers[1]
  }

  message("Using RNA layer: ", layer_use)

  LayerData(seurat_obj, assay = assay_name, layer = layer_use)
}

get_model_label <- function(model) {
  model_summary <- summary(model)

  r2 <- model_summary$r.squared

  age_p <- broom::tidy(model) %>%
    dplyr::filter(term == "age") %>%
    dplyr::pull(p.value)

  if (length(age_p) == 0) {
    age_p <- NA_real_
  }

  paste0(
    "R\u00b2 = ", round(r2, 3),
    "\n", format_p(age_p, prefix = "P")
  )
}

add_stats_label <- function(label, vjust = 1.1) {
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = label,
    hjust = -0.05,
    vjust = vjust,
    size = 4,
    fontface = "bold",
    label.size = NA,
    fill = "white",
    alpha = 0.85
  )
}

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size - 1, face = "bold"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      plot.margin = margin(12, 18, 12, 18)
    )
}

save_plot_pdf_png <- function(plot, filename_base, width, height, dpi = 600) {
  ggsave(
    filename = paste0(filename_base, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  )

  ggsave(
    filename = paste0(filename_base, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

extract_numeric_age <- function(metadata_df, age_col) {
  metadata_df %>%
    mutate(age = as.numeric(stringr::str_extract(.data[[age_col]], "\\d+"))) %>%
    filter(!is.na(age))
}


# -----------------------------
# 4. Main analysis
# -----------------------------

for (ff in seq_along(file_list)) {

  file_name <- file_list[[ff]]$name
  celltype_col <- file_list[[ff]]$celltype_col
  age_col <- file_list[[ff]]$age_col
  required_metadata <- file_list[[ff]]$required_metadata

  message("===================================================")
  message("Processing dataset: ", file_name)
  message("===================================================")

  input_file <- file.path(data_dir, file_name)

  if (!file.exists(input_file)) {
    warning(
      "Input file not found: ", input_file,
      "\nSkipping this dataset. Place the .h5ad file in data/ or set DATA_DIR."
    )
    next
  }

  dataset_out_dir <- file.path(results_dir, tools::file_path_sans_ext(file_name))
  fs::dir_create(dataset_out_dir)

  # Load .h5ad as Seurat object.
  seurat_obj <- tryCatch(
    {
      anndataR::read_h5ad(input_file, as = "Seurat")
    },
    error = function(e) {
      warning("Skipping ", file_name, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(seurat_obj)) {
    next
  }

  DefaultAssay(seurat_obj) <- "RNA"

  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop(
      "Cell type column '", celltype_col, "' not found. Available metadata columns are: ",
      paste(colnames(seurat_obj@meta.data), collapse = ", ")
    )
  }

  # Exclude selected cell types.
  celltype_values <- seurat_obj@meta.data[[celltype_col]]
  keep_cells <- !tolower(celltype_values) %in% tolower(exclude_celltypes)

  seurat_obj <- subset(
    seurat_obj,
    cells = colnames(seurat_obj)[keep_cells]
  )

  message("Remaining cells after cell type filtering: ", ncol(seurat_obj))

  # Optional UMAP overview, if the expected reduction exists.
  if ("X_umap" %in% Reductions(seurat_obj)) {
    p_umap <- DimPlot(
      seurat_obj,
      reduction = "X_umap",
      group.by = celltype_col,
      label = TRUE
    ) +
      ggtitle(paste0("Cell type overview: ", tools::file_path_sans_ext(file_name)))

    save_plot_pdf_png(
      plot = p_umap,
      filename_base = file.path(dataset_out_dir, "celltype_umap_overview"),
      width = 7,
      height = 6
    )
  }

  expr_mat <- get_expression_matrix(seurat_obj)
  feature_names <- rownames(expr_mat)
  gene_id_type <- detect_gene_id_type(feature_names)

  message("Detected gene identifier type: ", gene_id_type)

  # Harmonize tissue column once in metadata.
  meta_df <- seurat_obj@meta.data
  meta_df$cell_barcode <- rownames(meta_df)

  if (!"tissue" %in% colnames(meta_df) && "tissue_type" %in% colnames(meta_df)) {
    meta_df$tissue <- meta_df$tissue_type
  }

  required_cols <- unique(c(required_metadata, "donor_id", age_col, "sex", "tissue", celltype_col))
  missing_cols <- setdiff(required_cols, colnames(meta_df))

  if (length(missing_cols) > 0) {
    warning(
      "Skipping ", file_name,
      " because the following metadata columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
    next
  }

  # Standardize cell type column name for downstream analysis.
  if (celltype_col != "cell_type") {
    meta_df$cell_type <- meta_df[[celltype_col]]
  }

  if (age_col != "development_stage") {
    meta_df$development_stage <- meta_df[[age_col]]
    age_col <- "development_stage"
  }

  for (target_gene_symbol in analysis_gene_symbols) {

    message("---------------------------------------------------")
    message("Running target gene: ", target_gene_symbol)
    message("---------------------------------------------------")

    gene_info <- get_gene_ids(
      target_gene_symbol = target_gene_symbol,
      analysis_gene_symbols = analysis_gene_symbols,
      gene_map = gene_map,
      id_type = gene_id_type
    )

    target_gene <- gene_info$target_gene
    other_genes <- gene_info$other_genes
    all_genes <- gene_info$all_genes
    other_gene_symbols <- gene_info$other_gene_symbols
    gene_label_map <- gene_info$gene_label_map

    gene_check <- data.frame(
      target_gene_symbol = target_gene_symbol,
      gene_symbol = c(target_gene_symbol, other_gene_symbols),
      gene_id_used = all_genes,
      present = all_genes %in% rownames(expr_mat),
      stringsAsFactors = FALSE
    )

    print(gene_check)

    genes_found <- all_genes[all_genes %in% rownames(expr_mat)]

    if (!target_gene %in% genes_found) {
      warning(
        "Skipping target gene ", target_gene_symbol,
        " because it was not found in the expression matrix."
      )
      next
    }

    other_genes_found <- other_genes[other_genes %in% rownames(expr_mat)]

    if (length(other_genes_found) == 0) {
      warning(
        "No comparison genes were found for target gene ", target_gene_symbol,
        ". Age analysis will proceed, but positive-vs-negative comparison will be skipped."
      )
    }

    gene_out_dir <- file.path(dataset_out_dir, target_gene_symbol)
    fs::dir_create(gene_out_dir)

    genes_for_this_analysis <- c(target_gene, other_genes_found)

    expr_df <- as.data.frame(
      t(as.matrix(expr_mat[genes_for_this_analysis, colnames(seurat_obj), drop = FALSE]))
    )

    colnames(expr_df)[colnames(expr_df) == target_gene] <- "target_gene"
    expr_df$cell_barcode <- rownames(expr_df)

    df <- meta_df %>%
      left_join(expr_df, by = "cell_barcode") %>%
      extract_numeric_age(age_col = age_col)

    # ============================================================
    # Goal 1: Donor-level association with age across all cells
    # ============================================================

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
      vapply(covars, function(x) {
        x %in% colnames(donor_expr) &&
          n_distinct(na.omit(donor_expr[[x]])) >= 2
      }, logical(1))
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
    model_all_pct <- lm(formula_pct, data = donor_expr)

    label_mean <- get_model_label(model_all_mean)
    label_pct <- get_model_label(model_all_pct)

    p_all_mean <- ggplot(donor_expr, aes(x = age, y = mean_target_gene)) +
      geom_point(aes(size = n_cells), alpha = 0.75, color = "black") +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#2C7BB6",
        fill = "#ABD9E9",
        linewidth = 1
      ) +
      add_stats_label(label_mean, vjust = 3.5) +
      scale_size_continuous(name = "Cells per donor", range = c(2.5, 7)) +
      labs(
        title = paste0(target_gene_symbol, " mean expression vs age"),
        x = "Age",
        y = paste0("Mean ", target_gene_symbol, " expression per donor")
      ) +
      theme_pub()

    p_all_pct <- ggplot(donor_expr, aes(x = age, y = pct_target_gene_pos)) +
      geom_point(aes(size = n_cells), alpha = 0.75, color = "black") +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#D7191C",
        fill = "#FDAE61",
        linewidth = 1
      ) +
      add_stats_label(label_pct, vjust = 3.5) +
      scale_size_continuous(name = "Cells per donor", range = c(2.5, 7)) +
      labs(
        title = paste0(target_gene_symbol, "-positive cell percentage vs age"),
        x = "Age",
        y = paste0("% ", target_gene_symbol, "+ cells per donor")
      ) +
      theme_pub()

    save_plot_pdf_png(
      p_all_mean,
      file.path(gene_out_dir, paste0(target_gene_symbol, "_vs_age_all_donors_mean")),
      width = 6,
      height = 5
    )

    save_plot_pdf_png(
      p_all_pct,
      file.path(gene_out_dir, paste0(target_gene_symbol, "_vs_age_all_donors_percent_positive")),
      width = 6,
      height = 5
    )

    write.csv(
      donor_expr,
      file.path(gene_out_dir, paste0(target_gene_symbol, "_donor_level_expression.csv")),
      row.names = FALSE
    )

    # ============================================================
    # Goal 2: Donor-level association with age by cell type
    # ============================================================

    celltype_expr <- df %>%
      group_by(donor_id, age, sex, tissue, cell_type) %>%
      summarise(
        mean_target_gene = mean(.data[["target_gene"]], na.rm = TRUE),
        pct_target_gene_pos = mean(.data[["target_gene"]] > 0, na.rm = TRUE) * 100,
        n_cells = n(),
        .groups = "drop"
      ) %>%
      filter(n_cells >= min_cells_per_donor_celltype)

    celltype_plot_dir <- file.path(gene_out_dir, paste0(target_gene_symbol, "_celltype_plots"))
    mean_plot_dir <- file.path(celltype_plot_dir, "mean_expression")
    pct_plot_dir <- file.path(celltype_plot_dir, "percent_positive")

    fs::dir_create(mean_plot_dir)
    fs::dir_create(pct_plot_dir)

    get_lm_label <- function(data, y_var) {
      model <- lm(as.formula(paste0(y_var, " ~ age")), data = data)
      get_model_label(model)
    }

    make_celltype_plot <- function(data, y_var, y_label, title_text, line_color, fill_color) {
      label_text <- get_lm_label(data, y_var)

      ggplot(data, aes(x = age, y = .data[[y_var]])) +
        geom_point(aes(size = n_cells), alpha = 0.75, color = "black") +
        geom_smooth(
          method = "lm",
          se = TRUE,
          linewidth = 1,
          color = line_color,
          fill = fill_color
        ) +
        add_stats_label(label_text) +
        scale_size_continuous(name = "Cells per donor", range = c(2.5, 7)) +
        labs(title = title_text, x = "Age", y = y_label) +
        theme_pub()
    }

    celltype_stats <- celltype_expr %>%
      group_by(cell_type) %>%
      filter(n_distinct(donor_id) >= min_donors_per_celltype) %>%
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

    write.csv(
      celltype_stats,
      file.path(gene_out_dir, paste0(target_gene_symbol, "_age_correlation_by_celltype.csv")),
      row.names = FALSE
    )

    celltypes_to_plot <- celltype_stats %>%
      filter(n_donors >= min_donors_per_celltype) %>%
      pull(cell_type)

    for (ct in celltypes_to_plot) {
      plot_data <- celltype_expr %>%
        filter(cell_type == ct)

      ct_safe <- safe_filename(ct)

      p_mean <- make_celltype_plot(
        data = plot_data,
        y_var = "mean_target_gene",
        y_label = paste0("Mean ", target_gene_symbol, " expression per donor"),
        title_text = paste0(target_gene_symbol, " mean expression vs age\n", ct),
        line_color = "#2C7BB6",
        fill_color = "#ABD9E9"
      )

      p_pct <- make_celltype_plot(
        data = plot_data,
        y_var = "pct_target_gene_pos",
        y_label = paste0("% ", target_gene_symbol, "+ cells per donor"),
        title_text = paste0(target_gene_symbol, "-positive cell percentage vs age\n", ct),
        line_color = "#D7191C",
        fill_color = "#FDAE61"
      )

      save_plot_pdf_png(
        p_mean,
        file.path(mean_plot_dir, paste0(target_gene_symbol, "_", ct_safe, "_mean_expression_vs_age")),
        width = 6.5,
        height = 5
      )

      save_plot_pdf_png(
        p_pct,
        file.path(pct_plot_dir, paste0(target_gene_symbol, "_", ct_safe, "_percent_positive_vs_age")),
        width = 6.5,
        height = 5
      )
    }

    # ============================================================
    # Goal 3: Positive-vs-negative comparison
    # ============================================================

    genes_for_positive_analysis <- other_genes_found[other_genes_found %in% colnames(df)]

    if (length(genes_for_positive_analysis) == 0) {
      warning("Skipping positive-vs-negative analysis for ", target_gene_symbol, ".")
      next
    }

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

    positive_plot_dir <- file.path(gene_out_dir, paste0(target_gene_symbol, "_positive_vs_negative_plots"))
    fs::dir_create(positive_plot_dir)

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
      mutate(gene_symbol = gene_label_map[gene])

    p_pos_box <- ggplot(
      plot_df,
      aes(x = target_gene_status, y = expression, fill = target_gene_status)
    ) +
      geom_boxplot(width = 0.6, outlier.size = 0.15, linewidth = 0.5) +
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
      theme_pub() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)
      )

    save_plot_pdf_png(
      p_pos_box,
      file.path(positive_plot_dir, paste0(target_gene_symbol, "_positive_vs_negative_other_genes_cell_level_boxplot")),
      width = 10,
      height = 6
    )

    # Donor-level positive-vs-negative test.
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
      filter(n_cells >= min_cells_positive_negative_all)

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

        p_value = if (
          n_distinct(target_gene_status) == 2 &&
            n_positive >= 2 &&
            n_negative >= 2
        ) {
          wilcox.test(mean_expression ~ target_gene_status)$p.value
        } else {
          NA_real_
        },

        .groups = "drop"
      ) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        p_label = format_p(p_adj, prefix = "adj. P"),
        direction = case_when(
          p_adj < 0.05 & log2FC_positive_vs_negative > 0 ~ paste0("Higher in ", target_gene_symbol, "+ cells"),
          p_adj < 0.05 & log2FC_positive_vs_negative < 0 ~ paste0("Lower in ", target_gene_symbol, "+ cells"),
          TRUE ~ "No significant difference"
        )
      ) %>%
      arrange(p_adj)

    write.csv(
      positive_stats,
      file.path(gene_out_dir, paste0(target_gene_symbol, "_positive_population_other_genes_stats.csv")),
      row.names = FALSE
    )

    label_df <- donor_pos_long %>%
      group_by(gene, gene_symbol) %>%
      summarise(
        y_max = max(mean_expression, na.rm = TRUE),
        y_min = min(mean_expression, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        y_range = ifelse(y_max - y_min == 0, y_max, y_max - y_min),
        y_range = ifelse(y_range == 0, 0.1, y_range),
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
      geom_blank(
        data = label_df,
        aes(x = group1, y = y.expand),
        inherit.aes = FALSE
      ) +
      geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.6, alpha = 0.85) +
      geom_jitter(aes(size = n_cells), width = 0.12, alpha = 0.7, color = "black") +
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
      scale_size_continuous(name = "Cells per donor", range = c(1.8, 4.8)) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0("Donor-level mean expression in ", target_gene_symbol, "+ vs ", target_gene_symbol, "- cells"),
        subtitle = "Each point represents one donor-level mean; P values are BH-adjusted Wilcoxon tests",
        x = "",
        y = "Mean expression per donor"
      ) +
      theme_pub() +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.margin = margin(16, 24, 16, 24),
        panel.spacing = unit(1.2, "lines")
      )

    save_plot_pdf_png(
      p_pos_donor,
      file.path(positive_plot_dir, paste0(target_gene_symbol, "_positive_vs_negative_other_genes_donor_level")),
      width = 10.5,
      height = 7
    )

    # ============================================================
    # Goal 3b: Positive-vs-negative comparison by cell type
    # ============================================================

    celltype_pos_plot_dir <- file.path(
      gene_out_dir,
      paste0(target_gene_symbol, "_positive_vs_negative_by_celltype_plots")
    )

    fs::dir_create(celltype_pos_plot_dir)

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
      filter(n_cells >= min_cells_positive_negative_celltype) %>%
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

    celltype_positive_stats <- celltype_pos_long %>%
      group_by(cell_type, gene, gene_symbol) %>%
      filter(
        n_distinct(target_gene_status) == 2,
        n_distinct(donor_id) >= min_donors_per_celltype
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

        p_value = if (n_positive >= 2 && n_negative >= 2) {
          wilcox.test(mean_expression ~ target_gene_status)$p.value
        } else {
          NA_real_
        },

        .groups = "drop"
      ) %>%
      group_by(cell_type) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        p_label = format_p(p_adj, prefix = "adj. P"),
        direction = case_when(
          p_adj < 0.05 & log2FC_positive_vs_negative > 0 ~ paste0("Higher in ", target_gene_symbol, "+ cells"),
          p_adj < 0.05 & log2FC_positive_vs_negative < 0 ~ paste0("Lower in ", target_gene_symbol, "+ cells"),
          TRUE ~ "No significant difference"
        )
      ) %>%
      ungroup() %>%
      arrange(cell_type, p_adj)

    write.csv(
      celltype_positive_stats,
      file.path(gene_out_dir, paste0(target_gene_symbol, "_positive_population_other_genes_by_celltype.csv")),
      row.names = FALSE
    )

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
        geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.6, alpha = 0.85) +
        geom_jitter(aes(size = n_cells), width = 0.12, alpha = 0.7, color = "black") +
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
        scale_size_continuous(name = "Cells per donor", range = c(1.8, 4.8)) +
        coord_cartesian(clip = "off") +
        labs(
          title = paste0(target_gene_symbol, "+ vs ", target_gene_symbol, "- cells in ", ct),
          subtitle = "Each point represents one donor-level mean; P values are BH-adjusted Wilcoxon tests",
          x = "",
          y = "Mean expression per donor",
          fill = "Cell status"
        ) +
        theme_pub() +
        theme(
          legend.position = "right",
          axis.text.x = element_text(angle = 30, hjust = 1)
        )

      ct_safe <- safe_filename(ct)

      save_plot_pdf_png(
        p_ct,
        file.path(
          celltype_pos_plot_dir,
          paste0(target_gene_symbol, "_positive_vs_negative_other_genes_", ct_safe)
        ),
        width = 10.5,
        height = 6.8
      )
    }

    gc()
  }

  gc()
}

# Save R session information for reproducibility.
session_info_file <- file.path(results_dir, "sessionInfo.txt")
capture.output(sessionInfo(), file = session_info_file)

message("Analysis complete.")
message("Results written to: ", normalizePath(results_dir, mustWork = FALSE))
message("Session information written to: ", normalizePath(session_info_file, mustWork = FALSE))
