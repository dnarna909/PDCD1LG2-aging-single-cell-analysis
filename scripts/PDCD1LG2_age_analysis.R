#!/usr/bin/env Rscript

# PDCD1LG2 aging single-cell analysis ============================================================

#
# Goals:
#  1. Evaluate PDCD1LG2, MMP2, MMP9, CDKN1A, CDKN2A, IFNG, and TNF
#     expression across all cells. 
#
# GitHub-ready notes:
#   - Do not hard-code local computer paths.
#   - Place input .h5ad files in the data/ directory, or set DATA_DIR.
#   - Results are written to the results/ directory, or set RESULTS_DIR.
#
# Example:
#   Rscript scripts/PDCD1LG2_age_analysis.R
#
# Optional environment variables:
#   DATA_DIR=/path/to/h5ad/files RESULTS_DIR=/path/to/results Rscript scripts/PDCD1LG2_age_analysis.R




# 1. Package setup -----------------------------


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
  "ggpubr",
  "openxlsx"
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
  library(openxlsx)
})



# 2. User configuration -----------------------------


data_dir <- Sys.getenv("DATA_DIR", unset = "data")
results_dir <- Sys.getenv("RESULTS_DIR", unset = "results")

fs::dir_create(data_dir)
fs::dir_create(results_dir)

file_list <- list(
  list(
    name = "Healthy pediatric and adult human liver tissue.h5ad",
    celltype_col = "cell_type",
    age_col = "development_stage",
    required_metadata = c("donor_id", "development_stage", "sex", "cell_type")
  )
)

analysis_gene_symbols <- c(
  "PDCD1LG2",
  "MMP2",
  "MMP9",
  "CDKN1A",
  "CDKN2A",
  "IFNG",
  "TNF"
)

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

min_cells_per_donor_celltype <- 20
min_donors_per_celltype <- 6

genes.label <- c(
  "PDCD1LG2 (PD-L2)",
  "CDKN2A (p16)",
  "CDKN1A (p21)",
  "MMP2",
  "MMP9",
  "TNF",
  "IFNG"
)

genes.order <- c(
  "PDCD1LG2",
  "CDKN2A",
  "CDKN1A",
  "MMP2",
  "MMP9",
  "TNF",
  "IFNG"
)

cell.order <- c(
  "hematopoietic stem cell"  ,
  "endothelial cell of hepatic sinusoid",
  "hepatocyte"             ,
  "macrophage"      ,
  "Kupffer cell"      ,
  "cholangiocyte"     ,
  "natural killer cell"    ,
  "T cell"         ,
  "mononuclear phagocyte"     ,
  "mature B cell"       ,
  "plasma cell"       ,   
  "gamma-delta T cell"   ,
  "unknown",
  "All cells",
  "erythrocyte"                ,
  "platelet"      ,
  "neutrophil"       ,
  "conventional dendritic cell"     ,
  "cycling myeloid cell"        ,
  "cycling plasma cell"       
)

cell.label <- c(
  "Hematopoietic stem cell"  ,
  "Endothelial cell of hepatic sinusoid",
  "Hepatocyte"             ,
  "Macrophage"      ,
  "Kupffer cell"      ,
  "Cholangiocyte"     ,
  "Natural killer cell"    ,
  "T cell"         ,
  "Mononuclear phagocyte"     ,
  "Mature B cell"       ,
  "Plasma cell"       ,   
  "Gamma-delta T cell"   ,
  "Unknown",
  "All cells",
  
  "erythrocyte"                ,
  "platelet"      ,
  "neutrophil"       ,
  "conventional dendritic cell"     ,
  "cycling myeloid cell"        ,
  "cycling plasma cell"       
)

cell.label.map <- stats::setNames(cell.label, cell.order)
gene.label.map <- stats::setNames(genes.label, genes.order)

# Reverse row order for plotting
cell.order.rev <- rev(cell.order)
cell.label.rev <- rev(cell.label)
cell.label.map.rev <- stats::setNames(cell.label.rev, cell.order.rev)


# Two model sets will be generated throughout the analysis.
model_configs <- list(
  age_only = list(
    label = "Age-only model",
    suffix = "age_only",
    covars = character(0)
  ),
  age_sex = list(
    label = "Age + sex model",
    suffix = "age_sex_adjusted",
    covars = c("sex")
  )
)


# 3. Helper functions -----------------------------


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

get_valid_covars <- function(data, covars) {
  covars[
    vapply(covars, function(x) {
      x %in% colnames(data) &&
        dplyr::n_distinct(stats::na.omit(data[[x]])) >= 2
    }, logical(1))
  ]
}

make_lm_formula <- function(y_var, valid_covars) {
  rhs <- c("age", valid_covars)
  as.formula(paste(y_var, "~", paste(rhs, collapse = " + ")))
}

extract_age_model_stats <- function(model) {
  model_summary <- summary(model)
  age_row <- broom::tidy(model) %>%
    dplyr::filter(term == "age")
  
  if (nrow(age_row) == 0) {
    return(list(
      beta_age = NA_real_,
      r2 = model_summary$r.squared,
      adj_r2 = model_summary$adj.r.squared,
      p_age = NA_real_
    ))
  }
  
  list(
    beta_age = age_row$estimate[1],
    r2 = model_summary$r.squared,
    adj_r2 = model_summary$adj.r.squared,
    p_age = age_row$p.value[1]
  )
}

get_model_label <- function(model, p_adj) {
  model_summary <- summary(model)
  r2 <- model_summary$r.squared
  adj_r2 <- model_summary$adj.r.squared
  
  paste0(
    "R² = ", round(r2, 3),
    "\nAdjusted R² = ", round(adj_r2, 3),
    "\n", format_p(p_adj, prefix = "adj. P")
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

fit_donor_models <- function(data, model_config) {
  valid_covars <- get_valid_covars(data, model_config$covars)
  
  formula_mean <- make_lm_formula("mean_target_gene", valid_covars)
  formula_pct <- make_lm_formula("pct_target_gene_pos", valid_covars)
  
  model_mean <- lm(formula_mean, data = data)
  model_pct <- lm(formula_pct, data = data)
  
  stats_mean <- extract_age_model_stats(model_mean)
  stats_pct <- extract_age_model_stats(model_pct)
  
  p_adj <- p.adjust(
    c(mean_expression = stats_mean$p_age, percent_positive = stats_pct$p_age),
    method = "BH"
  )
  
  list(
    model_mean = model_mean,
    model_pct = model_pct,
    stats_mean = stats_mean,
    stats_pct = stats_pct,
    p_adj_mean = p_adj["mean_expression"],
    p_adj_pct = p_adj["percent_positive"],
    valid_covars = valid_covars,
    formula_mean = formula_mean,
    formula_pct = formula_pct
  )
}


get_adjusted_prediction_df <- function(data, model, valid_covars, y_var, n_points = 100) {
  age_seq <- seq(
    min(data$age, na.rm = TRUE),
    max(data$age, na.rm = TRUE),
    length.out = n_points
  )
  
  pred_df <- data.frame(age = age_seq)
  
  for (covar in valid_covars) {
    covar_values <- data[[covar]]
    
    if (is.numeric(covar_values)) {
      pred_df[[covar]] <- stats::median(covar_values, na.rm = TRUE)
    } else {
      ref_value <- names(sort(table(covar_values), decreasing = TRUE))[1]
      pred_df[[covar]] <- ref_value
      
      if (is.factor(covar_values)) {
        pred_df[[covar]] <- factor(pred_df[[covar]], levels = levels(covar_values))
      }
    }
  }
  
  pred_values <- predict(model, newdata = pred_df, interval = "confidence")
  pred_df$fit <- pred_values[, "fit"]
  pred_df$lwr <- pred_values[, "lwr"]
  pred_df$upr <- pred_values[, "upr"]
  pred_df
}

make_all_cells_stats_row <- function(donor_expr, target_gene_symbol, file_name, model_name, model_config) {
  model_results <- fit_donor_models(donor_expr, model_config)
  
  data.frame(
    cell_type = "All cells",
    n_donors = dplyr::n_distinct(donor_expr$donor_id),
    beta_age_mean = model_results$stats_mean$beta_age,
    r2_mean = model_results$stats_mean$r2,
    adj_r2_mean = model_results$stats_mean$adj_r2,
    p_age_mean = model_results$stats_mean$p_age,
    beta_age_pct = model_results$stats_pct$beta_age,
    r2_pct = model_results$stats_pct$r2,
    adj_r2_pct = model_results$stats_pct$adj_r2,
    p_age_pct = model_results$stats_pct$p_age,
    p_adj_mean = model_results$p_adj_mean,
    p_adj_pct = model_results$p_adj_pct,
    direction_mean = dplyr::case_when(
      model_results$p_adj_mean < 0.05 & model_results$stats_mean$beta_age > 0 ~ "Increase with age",
      model_results$p_adj_mean < 0.05 & model_results$stats_mean$beta_age < 0 ~ "Decrease with age",
      TRUE ~ "No significant change"
    ),
    direction_pct = dplyr::case_when(
      model_results$p_adj_pct < 0.05 & model_results$stats_pct$beta_age > 0 ~ "Increase with age",
      model_results$p_adj_pct < 0.05 & model_results$stats_pct$beta_age < 0 ~ "Decrease with age",
      TRUE ~ "No significant change"
    ),
    file_name = file_name,
    target_gene_symbol = target_gene_symbol,
    model_name = model_name,
    model_label = model_config$label,
    covariates_used = paste(model_results$valid_covars, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

calculate_celltype_stats <- function(celltype_expr, model_name, model_config, file_name) {
  celltypes <- celltype_expr %>%
    dplyr::group_by(cell_type) %>%
    dplyr::filter(dplyr::n_distinct(donor_id) >= min_donors_per_celltype) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(cell_type) %>%
    dplyr::pull(cell_type)
  
  stats_list <- lapply(celltypes, function(ct) {
    dat <- celltype_expr %>% dplyr::filter(cell_type == ct)
    model_results <- fit_donor_models(dat, model_config)
    
    data.frame(
      cell_type = ct,
      n_donors = dplyr::n_distinct(dat$donor_id),
      beta_age_mean = model_results$stats_mean$beta_age,
      r2_mean = model_results$stats_mean$r2,
      adj_r2_mean = model_results$stats_mean$adj_r2,
      p_age_mean = model_results$stats_mean$p_age,
      beta_age_pct = model_results$stats_pct$beta_age,
      r2_pct = model_results$stats_pct$r2,
      adj_r2_pct = model_results$stats_pct$adj_r2,
      p_age_pct = model_results$stats_pct$p_age,
      covariates_used = paste(model_results$valid_covars, collapse = ";"),
      file_name = file_name,
      model_name = model_name,
      model_label = model_config$label,
      stringsAsFactors = FALSE
    )
  })
  
  stats_df <- dplyr::bind_rows(stats_list)
  
  if (nrow(stats_df) == 0) {
    return(stats_df)
  }
  
  stats_df %>%
    dplyr::mutate(
      p_adj_mean = p.adjust(p_age_mean, method = "BH"),
      p_adj_pct = p.adjust(p_age_pct, method = "BH"),
      direction_mean = dplyr::case_when(
        p_adj_mean < 0.05 & beta_age_mean > 0 ~ "Increase with age",
        p_adj_mean < 0.05 & beta_age_mean < 0 ~ "Decrease with age",
        TRUE ~ "No significant change"
      ),
      direction_pct = dplyr::case_when(
        p_adj_pct < 0.05 & beta_age_pct > 0 ~ "Increase with age",
        p_adj_pct < 0.05 & beta_age_pct < 0 ~ "Decrease with age",
        TRUE ~ "No significant change"
      )
    ) %>%
    dplyr::arrange(p_adj_mean)
}

make_celltype_plot <- function(data, y_var, y_label, title_text, subtitle_text,
                               line_color, fill_color, model_config, p_adj) {
  model_results <- fit_donor_models(data, model_config)
  model_to_label <- if (y_var == "mean_target_gene") {
    model_results$model_mean
  } else {
    model_results$model_pct
  }
  
  label_text <- get_model_label(model_to_label, p_adj)
  
  pred_df <- get_adjusted_prediction_df(
    data = data,
    model = model_to_label,
    valid_covars = model_results$valid_covars,
    y_var = y_var
  )
  
  ggplot(data, aes(x = age, y = .data[[y_var]])) +
    geom_point(aes(size = n_cells), alpha = 0.75, color = "black") +
    geom_ribbon(
      data = pred_df,
      aes(x = age, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      fill = fill_color,
      alpha = 0.35
    ) +
    geom_line(
      data = pred_df,
      aes(x = age, y = fit),
      inherit.aes = FALSE,
      color = line_color,
      linewidth = 1
    ) +
    add_stats_label(label_text) +
    scale_size_continuous(name = "Cells per donor", range = c(2.5, 7)) +
    labs(title = title_text, subtitle = subtitle_text, x = "Age", y = y_label) +
    theme_pub()
}

make_sig_stars <- function(p_adj) {
  dplyr::case_when(
    is.na(p_adj) ~ "",
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

make_heatmap <- function(stats_df, value_col, p_adj_col, title_text, subtitle_text,
                         output_name, heatmap_dir, shared_celltype_order, analysis_gene_symbols_order,
                         fill_type = c("signed_adjusted_r2", "signed_r2", "beta", "r2") ) {
  fill_type <- match.arg(fill_type)
  
  heatmap_df <- stats_df %>%
    dplyr::mutate(
      value = .data[[value_col]],
      p_adj = .data[[p_adj_col]],
      significance = make_sig_stars(p_adj),
      target_gene_symbol = factor(target_gene_symbol, levels = analysis_gene_symbols_order),
      cell_type = factor(cell_type, levels = shared_celltype_order)
    ) %>%
    # dplyr::filter(!is.na(value)) %>%
    mutate(
      missing_label = dplyr::if_else(
        is.na(.data[[value_col]]) | .data[[value_col]] == 0,
        "X",
        ""
      )
    )
  
  max_abs_value <- max(abs(heatmap_df[[value_col]]), na.rm = TRUE)
  
  if (!is.finite(max_abs_value) || max_abs_value == 0) {
    max_abs_value <- 1
  }
  
  n_rows <- length(unique(heatmap_df$cell_type))
  
  legend_barheight <- grid::unit(
    max(3, 0.35 * n_rows),
    "cm"
  )
  
  fill_scale <- if (fill_type == "beta") {
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      limits = c(-max_abs_value*0.25, max_abs_value*0.25),
      breaks = c(-max_abs_value*0.25, 0, max_abs_value*0.25),
      labels = scales::number_format(accuracy = 0.001),
      na.value = "grey90",
      # name = "Age beta"
      name = NULL
    )
  } else if (fill_type == "r2") {
    scale_fill_gradient(
      low = "white",
      high = "#54278F",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = scales::number_format(accuracy = 0.01),
      na.value = "grey90",
      # name = "R²",
      name = NULL
    )
  } else if (fill_type %in% c("signed_r2", "signed_adjusted_r2")) {
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      # limits = c(-max_abs_value, max_abs_value),
      # breaks = c(-max_abs_value, 0, max_abs_value),
      # labels = scales::number_format(accuracy = 0.01),
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      na.value = "grey90",
      # name = dplyr::if_else(fill_type == "signed_adjusted_r2", "Signed adjusted R²", "Signed R²")
      name = NULL
    )
  }
  
  p_heatmap <- ggplot(
    heatmap_df,
    aes(x = target_gene_symbol, y = cell_type, fill = value)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = significance),
      color = "black",
      size = 4.5,
      fontface = "bold"
    ) +
    geom_text(
      ggplot2::aes(label = missing_label),
      size = 5,
      fontface = "bold",
      color = "black",
      na.rm = TRUE
    ) +
    fill_scale +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Target gene",
      y = "Cell type"
    ) +
    theme_pub(base_size = 12) +
    theme(
      legend.title = ggplot2::element_blank(),
      legend.key.height = legend_barheight,
      legend.key.width = grid::unit(0.35, "cm"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 12),
      panel.grid = element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    ) +
    guides(
      fill = ggplot2::guide_colorbar(
        title = NULL,
        # barheight = grid::unit(1, "npc"), # too tall
        barheight = legend_barheight,
        barwidth = grid::unit(0.35, "cm"),
        frame.colour = "black",
        frame.linewidth = 0.5,
        ticks.colour = "black"
      )
    ) 
  
  save_plot_pdf_png(
    p_heatmap,
    file.path(heatmap_dir, output_name),
    width = max(6, 0.9 * length(unique(heatmap_df$target_gene_symbol))),
    height = max(5, 0.35 * length(unique(heatmap_df$cell_type)))
  )
}



# 4. Main analysis-----------------------------


for (ff in seq_along(file_list)) {
  
  file_name <- file_list[[ff]]$name
  celltype_col <- file_list[[ff]]$celltype_col
  age_col <- file_list[[ff]]$age_col
  required_metadata <- file_list[[ff]]$required_metadata
  
  # 1) Processing dataset-----------------------------
  
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
  
  celltype_values <- seurat_obj@meta.data[[celltype_col]]
  keep_cells <- !tolower(celltype_values) %in% tolower(exclude_celltypes)
  
  seurat_obj <- subset(
    seurat_obj,
    cells = colnames(seurat_obj)[keep_cells]
  )
  
  message("Remaining cells after cell type filtering: ", ncol(seurat_obj))
  
  
  if ("X_umap" %in% Reductions(seurat_obj)) {
    p_umap <- DimPlot(
      seurat_obj,
      reduction = "X_umap",
      group.by = celltype_col,
      label = TRUE,
      repel = TRUE
    ) +
      ggtitle(paste0("Cell type overview: ", tools::file_path_sans_ext(file_name)))
    
    save_plot_pdf_png(
      plot = p_umap,
      filename_base = file.path(dataset_out_dir, "celltype_umap_overview"),
      width = 12,
      height = 6
    )
  }
  
  openxlsx::write.xlsx(
    list(
      subjects = as.data.frame.matrix(table(seurat_obj$development_stage, seurat_obj$donor_id)),
      cell_type = as.data.frame.matrix(table(seurat_obj$cell_type, seurat_obj$donor_id))
    ),
    file = file.path(dataset_out_dir, "celltype_summary_results.xlsx"),
    overwrite = TRUE,
    rowNames = TRUE
  )
  
  expr_mat <- get_expression_matrix(seurat_obj)
  
  # FeaturePlot() and VlnPlot() expect an RNA "data" layer.
  # AnnData imports commonly store the expression matrix as "X".
  # Copy the same layer used in the analysis into "data" for plotting.
  assay_name <- DefaultAssay(seurat_obj)
  rna_layers <- Layers(seurat_obj[[assay_name]])
  
  if (!"data" %in% rna_layers) {
    source_layer <- if ("X" %in% rna_layers) {
      "X"
    } else if ("counts" %in% rna_layers) {
      "counts"
    } else {
      rna_layers[1]
    }
    
    message(
      "RNA 'data' layer not found; creating it from layer: ",
      source_layer
    )
    
    LayerData(
      seurat_obj,
      assay = assay_name,
      layer = "data"
    ) <- LayerData(
      seurat_obj,
      assay = assay_name,
      layer = source_layer
    )
  }
  
  feature_names <- rownames(expr_mat)
  gene_id_type <- detect_gene_id_type(feature_names)
  
  message("Detected gene identifier type: ", gene_id_type)
  
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
  
  if (celltype_col != "cell_type") {
    meta_df$cell_type <- meta_df[[celltype_col]]
  }
  
  if (age_col != "development_stage") {
    meta_df$development_stage <- meta_df[[age_col]]
    age_col <- "development_stage"
  }
  
  all_gene_celltype_stats <- list(
    age_only = list(),
    age_sex = list()
  )
  all_gene_all_cells_stats <- list(
    age_only = list(),
    age_sex = list()
  )
  
  # 2) Processing gene-----------------------------
  for (target_gene_symbol in analysis_gene_symbols) {
    message("Running target gene: ", target_gene_symbol)

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
    
    gene_out_dir <- file.path(dataset_out_dir, target_gene_symbol)
    fs::dir_create(gene_out_dir)
    
    genes_for_this_analysis <- c(target_gene, other_genes_found)
    
    expr_df <- as.data.frame(
      t(as.matrix(expr_mat[genes_for_this_analysis, colnames(seurat_obj), drop = FALSE]))
    )
    
    colnames(expr_df)[colnames(expr_df) == target_gene] <- "target_gene"
    expr_df$cell_barcode <- rownames(expr_df)
    
    meta_df <- meta_df %>%
      mutate(
        development_stage = case_when(
          donor_id == "C64" ~ "15",
          donor_id == "C70" ~ "48",
          donor_id == "C68" ~ "61",
          donor_id == "C39" ~ "65",
          donor_id == "C54" ~ "67",
          TRUE ~ development_stage  # keeps existing values; use NA if creating a new column
        )
      )
    table(meta_df$donor_id)
    table(meta_df$donor_id, meta_df$development_stage)
    
    df <- meta_df %>%
      left_join(expr_df, by = "cell_barcode") %>%
      extract_numeric_age(age_col = age_col)
    table(df$donor_id, df$development_stage)
    
    # Goal: FeaturePlot and violin plot of each target gene across all cells ---------------------------------------------------
    # These plots use the full filtered Seurat object and are saved as PNG files.
    goal3_dir <- file.path(dataset_out_dir, "goal3_feature_violin_plots")
    featureplot_dir <- file.path(goal3_dir, "feature_plots")
    violinplot_dir <- file.path(goal3_dir, "violin_plots")
    fs::dir_create(featureplot_dir)
    fs::dir_create(violinplot_dir)
    
    feature_reduction <- if ("X_umap" %in% Reductions(seurat_obj)) {
      "X_umap"
    } else if ("umap" %in% Reductions(seurat_obj)) {
      "umap"
    } else {
      NA_character_
    }
    
    summary(FetchData(seurat_obj, vars = target_gene)[,1])
    
    table(FetchData(seurat_obj, vars = target_gene)[,1] > 0)
    
    mean(FetchData(seurat_obj, vars = target_gene)[,1] > 0)
    
    if (is.na(feature_reduction)) {
      warning(
        "Skipping FeaturePlot for ", target_gene_symbol,
        " because no UMAP reduction ('X_umap' or 'umap') was found."
      )
    } else {
      p_feature <- FeaturePlot(
        object = seurat_obj,
        features = target_gene,
        keep.scale = "feature",
        reduction = feature_reduction,
        order = TRUE,
        raster = TRUE,
        # pt.size = 0.15,
        min.cutoff = NA,
       # max.cutoff = "q99",
        cols = c("grey90", "#6A00FF")
      ) +
        ggtitle(paste0(target_gene_symbol, " expression across all cells")) +
        theme_void(base_size = 12) +
        theme(
          plot.title = element_text(
            face = "bold",
            hjust = 0.5,
            size = 14
          ),
          legend.position = "right",
          plot.margin = margin(10, 10, 10, 10)
        )
      p_feature
      ggsave(
        filename = file.path(
          featureplot_dir,
          paste0(target_gene_symbol, "_featureplot_all_cells.png")
        ),
        plot = p_feature,
        width = 8,
        height = 6,
        dpi = 600
      )
    }
    
    p_violin <- VlnPlot(
      object = seurat_obj,
      features = target_gene,
      group.by = celltype_col,
      pt.size = 0
    ) +
      ggtitle(paste0(target_gene_symbol, " expression by cell type")) +
      xlab("Cell type") +
      ylab(paste0(target_gene_symbol, " expression")) +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    p_violin
    ggsave(
      filename = file.path(
        violinplot_dir,
        paste0(target_gene_symbol, "_violinplot_all_cells.png")
      ),
      plot = p_violin,
      width = 12,
      height = 6,
      dpi = 600
    )
    
    
  }
  
}

# Save R session information for reproducibility.
session_info_file <- file.path(results_dir, "sessionInfo.txt")
capture.output(sessionInfo(), file = session_info_file)

message("Analysis complete.")
message("Results written to: ", normalizePath(results_dir, mustWork = FALSE))
message("Session information written to: ", normalizePath(session_info_file, mustWork = FALSE))


message("Analysis complete.")
message("Results written to: ", normalizePath(results_dir, mustWork = FALSE))
message("Session information written to: ", normalizePath(session_info_file, mustWork = FALSE))

