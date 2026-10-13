#!/usr/bin/env Rscript

###############################################################################
# Anastatica RNA-seq post-processing, QC, and design-diagnostics script
#
# Author: Iris Sammarco
#
# Date: 25/08/2026
#
# Purpose:
#   Models post-processing QC and design-diagnostics without replacing files produced
#   by the main dream/variancePartition DEG script.
#
# This script reads existing results under degs_root and writes new summaries,
# QC outputs, metadata/design diagnostics, and plots. It does not refit DE
# models and does not delete or overwrite the primary model-result files.
#
# Main output roots:
#   degs_summary_root:
#     DEG summaries, variance-partition summaries, and merged QC notes.
#   correlation_root:
#     Metadata, confounding, canonical-correlation, MDS, and sample-structure
#     diagnostics.
#
# Statistical interpretation:
#   - FDR and log2-fold-change counts are descriptive summaries of already
#     computed tests; they do not constitute additional hypothesis tests.
#   - Omnibus F-tests are summarized without applying a logFC threshold because
#     a joint F-test does not have one single interpretable logFC.
#   - Variance-partition summaries describe distributions of per-gene fractions;
#     they are not sample-level variance estimates.
#   - z-score flags for residual variance and voom weights are exploratory QC
#     flags, not automatic exclusion criteria.
#   - Contingency tables, correlations, canCorPairs(), and MDS are diagnostics,
#     not proof of causality or independent biological effects.
###############################################################################

###############################################################################
# Package-library configuration
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(variancePartition)
})

###############################################################################
# Configuration
###############################################################################

degs_root <- "/Anastatica/degs"
degs_summary_root <- "/Anastatica/degs/degs_summary"
correlation_root <- "/Anastatica/degs/models_variable_correlations"
gse_rds <- "/Anastatica/tximeta/tximeta_gene_level.rds"

fdr_thresh <- 0.05
lfc_thresh <- log2(1.5)
z_threshold <- 3
mds_top_genes <- 500

# Process only the seven models listed below.
# Missing model outputs are logged and handled without stopping the run.
model_names <- c(
  "climate_model",
  "ECe_model",
  "distribution_range_model",
  "distribution_range_east_west_model",
  "ECe_climate_range_model",
  "ECe_climate_range_east_west_model",
  "ECe_climate_range_interaction_east_west_model"
)

model_specs <- list(
  climate_model = list(
    label = "Climate model",
    grouping_var = NULL,
    enabled = TRUE,
    climate = TRUE,
    combined = FALSE
  ),
  ECe_model = list(
    label = "ECe model",
    grouping_var = NULL,
    enabled = TRUE,
    climate = FALSE,
    combined = TRUE
  ),
  distribution_range_model = list(
    label = "Distribution range: East, Centre, West",
    grouping_var = "Site_of_distribution_range",
    enabled = TRUE,
    climate = FALSE,
    combined = FALSE
  ),
  distribution_range_east_west_model = list(
    label = "Distribution range: East, West",
    grouping_var = "Site_of_distribution_range_East_West",
    enabled = TRUE,
    climate = FALSE,
    combined = FALSE
  ),
  ECe_climate_range_model = list(
    label = "ECe plus climate plus distribution range",
    grouping_var = "Site_of_distribution_range",
    enabled = TRUE,
    climate = TRUE,
    combined = TRUE
  ),
  ECe_climate_range_east_west_model = list(
    label = "ECe plus climate plus distribution range: East, West",
    grouping_var = "Site_of_distribution_range_East_West",
    enabled = TRUE,
    climate = TRUE,
    combined = TRUE
  ),
  ECe_climate_range_interaction_east_west_model = list(
    label = "ECe plus climate plus East-West distribution-range interaction model",
    grouping_var = "Site_of_distribution_range_East_West",
    enabled = TRUE,
    climate = TRUE,
    combined = TRUE,
    interaction = TRUE
  )
)

# Create only the new output directories. Existing primary DEG files remain
# untouched; outputs created by this script have distinct names.
dir.create(degs_summary_root, recursive = TRUE, showWarnings = FALSE)
dir.create(correlation_root, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Unified run logging
###############################################################################

run_notes <- character()

add_note <- function(...) {
  txt <- paste0(...)
  run_notes <<- c(run_notes, txt)
  message(txt)
}

###############################################################################
# General helper functions
###############################################################################

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    add_note("MISSING FILE: ", path)
    return(NULL)
  }

  ans <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      add_note("READ ERROR: ", path, " -- ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(ans)) return(NULL)

  if (nrow(ans) == 0) {
    add_note("EMPTY FILE: ", path)
  } else {
    add_note("READ OK: ", path, " (rows=", nrow(ans), ", cols=", ncol(ans), ")")
  }
  ans
}

safe_write_csv <- function(x, path) {
  tryCatch({
    write.csv(x, path, row.names = FALSE, na = "")
    add_note("WROTE: ", path)
  }, error = function(e) {
    add_note("WRITE ERROR: ", path, " -- ", conditionMessage(e))
  })
}

finite_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x)]
}

safe_stat <- function(fun, x) {
  x <- finite_numeric(x)
  if (length(x) == 0) return(NA_real_)
  fun(x)
}

safe_mean <- function(x) safe_stat(mean, x)
safe_median <- function(x) safe_stat(median, x)
safe_sd <- function(x) {
  x <- finite_numeric(x)
  if (length(x) < 2) return(NA_real_)
  sd(x)
}
safe_min <- function(x) safe_stat(min, x)
safe_max <- function(x) safe_stat(max, x)

zflag <- function(x, threshold = 3) {
  x <- suppressWarnings(as.numeric(x))
  z <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  if (sum(keep) < 2) return(z)
  s <- sd(x[keep])
  if (!is.finite(s) || s == 0) return(z)
  z[keep] <- (x[keep] - mean(x[keep])) / s
  z
}

summarize_numeric <- function(x) {
  x <- finite_numeric(x)
  if (length(x) == 0) {
    return(c(n = 0, mean = NA_real_, median = NA_real_, sd = NA_real_,
             min = NA_real_, max = NA_real_))
  }
  c(n = length(x), mean = mean(x), median = median(x), sd = safe_sd(x),
    min = min(x), max = max(x))
}

order_df <- function(df, cols) {
  if (is.null(df) || nrow(df) == 0) return(df)
  cols <- intersect(cols, colnames(df))
  if (length(cols) == 0) return(df)
  ord <- do.call(order, c(df[cols], list(na.last = TRUE)))
  df[ord, , drop = FALSE]
}

###############################################################################
# DEG summaries
###############################################################################

summarize_de_table <- function(df, model_name, test_type, source_file) {
  if (is.null(df)) {
    return(data.frame(
      Model_name = model_name, Test_type = test_type, Test_label = NA_character_,
      Source_file = basename(source_file), Output_status = "missing",
      N_genes_tested = NA_integer_, N_sig_FDR05 = NA_integer_,
      N_sig_FDR05_log2FC_1.5 = NA_integer_, stringsAsFactors = FALSE
    ))
  }

  if (nrow(df) == 0) {
    return(data.frame(
      Model_name = model_name, Test_type = test_type, Test_label = NA_character_,
      Source_file = basename(source_file), Output_status = "empty",
      N_genes_tested = 0L, N_sig_FDR05 = 0L,
      N_sig_FDR05_log2FC_1.5 = NA_integer_, stringsAsFactors = FALSE
    ))
  }

  if (!"adj.P.Val" %in% colnames(df)) {
    add_note("INVALID DEG TABLE: missing adj.P.Val in ", source_file)
    return(data.frame(
      Model_name = model_name, Test_type = test_type, Test_label = NA_character_,
      Source_file = basename(source_file), Output_status = "invalid_missing_adj.P.Val",
      N_genes_tested = NA_integer_, N_sig_FDR05 = NA_integer_,
      N_sig_FDR05_log2FC_1.5 = NA_integer_, stringsAsFactors = FALSE
    ))
  }

  if (!"Test" %in% colnames(df)) {
    df$Test <- basename(source_file)
    add_note("INFO: Test column absent in ", source_file,
             "; file name used as test label.")
  }

  labels <- unique(as.character(df$Test))
  labels <- labels[!is.na(labels)]
  if (length(labels) == 0) labels <- NA_character_

  has_lfc <- "logFC" %in% colnames(df)

  out <- lapply(labels, function(label) {
    sub <- if (is.na(label)) df else df[as.character(df$Test) == label, , drop = FALSE]
    adj <- suppressWarnings(as.numeric(sub$adj.P.Val))
    tested <- is.finite(adj)
    n_fdr <- sum(tested & adj <= fdr_thresh)

    n_fdr_lfc <- if (has_lfc) {
      lfc <- suppressWarnings(as.numeric(sub$logFC))
      sum(tested & is.finite(lfc) & adj <= fdr_thresh & abs(lfc) >= lfc_thresh)
    } else {
      NA_integer_
    }

    data.frame(
      Model_name = model_name,
      Test_type = test_type,
      Test_label = label,
      Source_file = basename(source_file),
      Output_status = "processed",
      N_genes_tested = as.integer(sum(tested)),
      N_sig_FDR05 = as.integer(n_fdr),
      N_sig_FDR05_log2FC_1.5 = as.integer(n_fdr_lfc),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

run_deg_summaries <- function() {
  coeff <- list()
  contrasts <- list()
  omnibus <- list()

  for (m in model_names) {
    mdir <- file.path(degs_root, m)
    paths <- list(
      coefficient = file.path(mdir, paste0(m, "_DE_results_full.csv")),
      contrast = file.path(mdir, paste0(m, "_pairwise_contrasts_full.csv")),
      omnibus_Ftest = file.path(mdir, paste0(m, "_omnibus_Ftest_full.csv"))
    )

    coeff[[m]] <- summarize_de_table(safe_read_csv(paths$coefficient), m,
                                     "coefficient", paths$coefficient)
    contrasts[[m]] <- summarize_de_table(safe_read_csv(paths$contrast), m,
                                          "contrast", paths$contrast)
    omnibus[[m]] <- summarize_de_table(safe_read_csv(paths$omnibus_Ftest), m,
                                       "omnibus_Ftest", paths$omnibus_Ftest)
  }

  coeff_df <- order_df(do.call(rbind, coeff), c("Model_name", "Test_label"))
  contrast_df <- order_df(do.call(rbind, contrasts), c("Model_name", "Test_label"))
  omnibus_df <- order_df(do.call(rbind, omnibus), c("Model_name", "Test_label"))
  all_df <- rbind(coeff_df, contrast_df, omnibus_df)
  all_df <- order_df(all_df, c("Model_name", "Test_type", "Test_label"))

  safe_write_csv(coeff_df, file.path(degs_summary_root, "DEG_summary_coefficients.csv"))
  safe_write_csv(contrast_df, file.path(degs_summary_root, "DEG_summary_pairwise_contrasts.csv"))
  safe_write_csv(omnibus_df, file.path(degs_summary_root, "DEG_summary_omnibus_Ftests.csv"))
  safe_write_csv(all_df, file.path(degs_summary_root, "DEG_summary_all_tests.csv"))
}

###############################################################################
# Variance-partition summaries
###############################################################################

guess_feature_id_col <- function(df) {
  preferred <- c("FeatureID", "gene", "Gene", "GeneID", "gene_id", "X")
  hit <- preferred[preferred %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (ncol(df) > 0 && !is.numeric(df[[1]])) return(colnames(df)[1])
  NULL
}

summarize_varpart_model <- function(df, model_name, source_file) {
  if (is.null(df)) {
    return(list(
      by_term = data.frame(Model_name = model_name, Term = NA_character_,
                           Source_file = basename(source_file), Output_status = "missing",
                           N_genes = NA_integer_, Mean = NA_real_, Median = NA_real_,
                           SD = NA_real_, Min = NA_real_, Max = NA_real_,
                           Prop_gt_0 = NA_real_, Prop_ge_0.05 = NA_real_,
                           Prop_ge_0.10 = NA_real_, stringsAsFactors = FALSE),
      overview = data.frame(Model_name = model_name,
                            Source_file = basename(source_file), Output_status = "missing",
                            N_genes = NA_integer_, N_terms = NA_integer_,
                            Top_term_by_mean = NA_character_, Top_term_mean = NA_real_,
                            Top_term_by_median = NA_character_, Top_term_median = NA_real_,
                            stringsAsFactors = FALSE)
    ))
  }

  if (nrow(df) == 0) {
    return(list(
      by_term = data.frame(Model_name = model_name, Term = NA_character_,
                           Source_file = basename(source_file), Output_status = "empty",
                           N_genes = 0L, Mean = NA_real_, Median = NA_real_, SD = NA_real_,
                           Min = NA_real_, Max = NA_real_, Prop_gt_0 = NA_real_,
                           Prop_ge_0.05 = NA_real_, Prop_ge_0.10 = NA_real_,
                           stringsAsFactors = FALSE),
      overview = data.frame(Model_name = model_name,
                            Source_file = basename(source_file), Output_status = "empty",
                            N_genes = 0L, N_terms = 0L,
                            Top_term_by_mean = NA_character_, Top_term_mean = NA_real_,
                            Top_term_by_median = NA_character_, Top_term_median = NA_real_,
                            stringsAsFactors = FALSE)
    ))
  }

  id_col <- guess_feature_id_col(df)
  term_df <- if (is.null(id_col)) df else df[, setdiff(colnames(df), id_col), drop = FALSE]
  numeric_cols <- vapply(term_df, is.numeric, logical(1))
  term_df <- term_df[, numeric_cols, drop = FALSE]

  if (ncol(term_df) == 0) {
    add_note("INVALID VARPART TABLE: no numeric terms in ", source_file)
    return(summarize_varpart_model(data.frame(), model_name, source_file))
  }

  by_term <- lapply(colnames(term_df), function(term) {
    x <- suppressWarnings(as.numeric(term_df[[term]]))
    finite <- is.finite(x)
    xf <- x[finite]
    n <- length(xf)

    data.frame(
      Model_name = model_name,
      Term = term,
      Source_file = basename(source_file),
      Output_status = "processed",
      N_genes = as.integer(n),
      Mean = if (n > 0) mean(xf) else NA_real_,
      Median = if (n > 0) median(xf) else NA_real_,
      SD = if (n > 1) sd(xf) else NA_real_,
      Min = if (n > 0) min(xf) else NA_real_,
      Max = if (n > 0) max(xf) else NA_real_,
      Prop_gt_0 = if (n > 0) mean(xf > 0) else NA_real_,
      Prop_ge_0.05 = if (n > 0) mean(xf >= 0.05) else NA_real_,
      Prop_ge_0.10 = if (n > 0) mean(xf >= 0.10) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  by_term <- do.call(rbind, by_term)
  valid <- by_term[is.finite(by_term$Mean) | is.finite(by_term$Median), , drop = FALSE]

  if (nrow(valid) == 0) {
    overview <- data.frame(
      Model_name = model_name, Source_file = basename(source_file),
      Output_status = "processed", N_genes = nrow(df), N_terms = nrow(by_term),
      Top_term_by_mean = NA_character_, Top_term_mean = NA_real_,
      Top_term_by_median = NA_character_, Top_term_median = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    by_mean <- valid[order(-valid$Mean, valid$Term, na.last = TRUE), , drop = FALSE]
    by_median <- valid[order(-valid$Median, valid$Term, na.last = TRUE), , drop = FALSE]
    overview <- data.frame(
      Model_name = model_name, Source_file = basename(source_file),
      Output_status = "processed", N_genes = nrow(df), N_terms = nrow(by_term),
      Top_term_by_mean = by_mean$Term[1], Top_term_mean = by_mean$Mean[1],
      Top_term_by_median = by_median$Term[1], Top_term_median = by_median$Median[1],
      stringsAsFactors = FALSE
    )
  }

  list(by_term = by_term, overview = overview)
}

extract_top_terms <- function(by_term, top_n = 3) {
  if (is.null(by_term) || nrow(by_term) == 0) return(data.frame())
  by_term <- by_term[by_term$Output_status == "processed" & !is.na(by_term$Term), , drop = FALSE]
  if (nrow(by_term) == 0) return(data.frame())

  out <- lapply(unique(by_term$Model_name), function(m) {
    sub <- by_term[by_term$Model_name == m, , drop = FALSE]
    a <- head(sub[order(-sub$Mean, sub$Term, na.last = TRUE), , drop = FALSE], top_n)
    b <- head(sub[order(-sub$Median, sub$Term, na.last = TRUE), , drop = FALSE], top_n)
    rbind(
      data.frame(Model_name = m, Rank_type = "mean", Rank = seq_len(nrow(a)),
                 Term = a$Term, Value = a$Mean, stringsAsFactors = FALSE),
      data.frame(Model_name = m, Rank_type = "median", Rank = seq_len(nrow(b)),
                 Term = b$Term, Value = b$Median, stringsAsFactors = FALSE)
    )
  })
  do.call(rbind, out)
}

run_varpart_summaries <- function() {
  by_term_list <- list()
  overview_list <- list()

  for (m in model_names) {
    path <- file.path(degs_root, m, paste0(m, "_varPart_results.csv"))
    res <- summarize_varpart_model(safe_read_csv(path), m, path)
    by_term_list[[m]] <- res$by_term
    overview_list[[m]] <- res$overview
  }

  by_term <- order_df(do.call(rbind, by_term_list), c("Model_name", "Term"))
  overview <- order_df(do.call(rbind, overview_list), c("Model_name"))
  top <- order_df(extract_top_terms(by_term), c("Model_name", "Rank_type", "Rank"))

  safe_write_csv(by_term, file.path(degs_summary_root, "varPart_summary_by_term.csv"))
  safe_write_csv(top, file.path(degs_summary_root, "varPart_summary_top_terms.csv"))
  safe_write_csv(overview, file.path(degs_summary_root, "varPart_summary_model_overview.csv"))
}

###############################################################################
# Residual and voom-weight QC post-processing
###############################################################################

run_model_qc <- function(model_name) {
  outdir <- file.path(degs_root, model_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  residual_path <- file.path(outdir, paste0(model_name, "_diagnostics_residual_summary.csv"))
  weights_path <- file.path(outdir, paste0(model_name, "_diagnostics_voom_weights_summary.csv"))
  residual_df <- safe_read_csv(residual_path)
  weights_df <- safe_read_csv(weights_path)
  qc_lines <- c(paste0("QC summary for model: ", model_name), paste0("Folder: ", outdir), "")

  if (!is.null(residual_df) && "Residual_variance" %in% colnames(residual_df)) {
    residual_df$z_resid_var <- zflag(residual_df$Residual_variance, z_threshold)
    residual_df$flag_high_resid_var <- is.finite(residual_df$z_resid_var) &
      residual_df$z_resid_var >= z_threshold
    flags <- residual_df[residual_df$flag_high_resid_var, , drop = FALSE]
    flags <- flags[order(-residual_df$Residual_variance[match(rownames(flags), rownames(residual_df))]), , drop = FALSE]

    safe_write_csv(residual_df, file.path(outdir, paste0(model_name, "_diagnostics_residual_summary_with_flags.csv")))
    safe_write_csv(flags, file.path(outdir, paste0(model_name, "_diagnostics_residual_genes_flagged.csv")))

    rs <- summarize_numeric(residual_df$Residual_variance)
    qc_lines <- c(qc_lines,
      "Residual variance summary:",
      paste0("  Finite genes: ", rs["n"]),
      paste0("  Mean: ", signif(rs["mean"], 5)),
      paste0("  Median: ", signif(rs["median"], 5)),
      paste0("  SD: ", signif(rs["sd"], 5)),
      paste0("  Min: ", signif(rs["min"], 5)),
      paste0("  Max: ", signif(rs["max"], 5)),
      paste0("  Exploratory high-variance flags (z >= ", z_threshold, "): ", nrow(flags)), "")

    finite_res <- residual_df[is.finite(residual_df$Residual_variance), , drop = FALSE]
    if (nrow(finite_res) > 0) {
      p <- ggplot(finite_res, aes(Residual_variance)) +
        geom_histogram(bins = 50, fill = "grey70", color = "white") +
        theme_bw() + labs(title = paste(model_name, "- residual variance"),
                          x = "Residual variance", y = "Gene count")
      ggsave(file.path(outdir, paste0(model_name, "_diagnostics_residual_variance_histogram.pdf")), p, width = 7, height = 5)

      positive <- finite_res[finite_res$Residual_variance > 0, , drop = FALSE]
      if (nrow(positive) > 0) {
        p <- ggplot(positive, aes(log10(Residual_variance))) +
          geom_histogram(bins = 50, fill = "steelblue", color = "white") +
          theme_bw() + labs(title = paste(model_name, "- log10 residual variance"),
                            x = "log10 residual variance", y = "Gene count")
        ggsave(file.path(outdir, paste0(model_name, "_diagnostics_residual_variance_log10_histogram.pdf")), p, width = 7, height = 5)
      }
    }
  } else {
    qc_lines <- c(qc_lines, "Residual summary CSV unavailable or missing Residual_variance.", "")
  }

  if (!is.null(weights_df) && all(c("Sample", "Mean_weight", "SD_weight") %in% colnames(weights_df))) {
    weights_df$z_mean_weight <- zflag(weights_df$Mean_weight, z_threshold)
    weights_df$z_sd_weight <- zflag(weights_df$SD_weight, z_threshold)
    weights_df$flag_low_mean_weight <- is.finite(weights_df$z_mean_weight) & weights_df$z_mean_weight <= -z_threshold
    weights_df$flag_high_sd_weight <- is.finite(weights_df$z_sd_weight) & weights_df$z_sd_weight >= z_threshold
    weights_df$flag_any_weight_issue <- weights_df$flag_low_mean_weight | weights_df$flag_high_sd_weight
    flags <- weights_df[weights_df$flag_any_weight_issue, , drop = FALSE]
    flags <- flags[order(flags$Mean_weight, na.last = TRUE), , drop = FALSE]

    safe_write_csv(weights_df, file.path(outdir, paste0(model_name, "_diagnostics_voom_weights_summary_with_flags.csv")))
    safe_write_csv(flags, file.path(outdir, paste0(model_name, "_diagnostics_voom_weight_samples_flagged.csv")))

    wm <- summarize_numeric(weights_df$Mean_weight)
    ws <- summarize_numeric(weights_df$SD_weight)
    qc_lines <- c(qc_lines,
      "Voom-weight summary:",
      paste0("  Samples: ", wm["n"]),
      paste0("  Mean of mean weights: ", signif(wm["mean"], 5)),
      paste0("  Median of mean weights: ", signif(wm["median"], 5)),
      paste0("  Mean of weight SDs: ", signif(ws["mean"], 5)),
      paste0("  Median of weight SDs: ", signif(ws["median"], 5)),
      paste0("  Exploratory low-mean-weight flags: ", sum(weights_df$flag_low_mean_weight, na.rm = TRUE)),
      paste0("  Exploratory high-SD flags: ", sum(weights_df$flag_high_sd_weight, na.rm = TRUE)), "")

    p <- ggplot(weights_df, aes(reorder(Sample, Mean_weight), Mean_weight, fill = flag_any_weight_issue)) +
      geom_col() + coord_flip() + theme_bw() +
      scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
      labs(title = paste(model_name, "- mean voom weight"), x = "Sample", y = "Mean weight", fill = "Flagged")
    ggsave(file.path(outdir, paste0(model_name, "_diagnostics_voom_mean_weight_by_sample.pdf")), p, width = 9, height = 5)

    p <- ggplot(weights_df, aes(reorder(Sample, SD_weight), SD_weight, fill = flag_any_weight_issue)) +
      geom_col() + coord_flip() + theme_bw() +
      scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
      labs(title = paste(model_name, "- voom-weight SD"), x = "Sample", y = "Weight SD", fill = "Flagged")
    ggsave(file.path(outdir, paste0(model_name, "_diagnostics_voom_weight_sd_by_sample.pdf")), p, width = 9, height = 5)

    p <- ggplot(weights_df, aes(Mean_weight, SD_weight, color = flag_any_weight_issue)) +
      geom_point(alpha = 0.8) + theme_bw() +
      scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick")) +
      labs(title = paste(model_name, "- mean vs SD of voom weights"), x = "Mean weight", y = "Weight SD", color = "Flagged")
    ggsave(file.path(outdir, paste0(model_name, "_diagnostics_voom_weights_mean_vs_sd.pdf")), p, width = 6, height = 5)
  } else {
    qc_lines <- c(qc_lines, "Voom weights CSV unavailable or missing expected columns.", "")
  }

  qc_lines <- c(qc_lines,
    "Interpretation guide:",
    "- These are exploratory diagnostics, not automatic exclusion rules.",
    "- Inspect flagged genes or samples together with library size, sample identity, and biological metadata.",
    "- A right-skewed residual-variance distribution is not automatically problematic.")
  writeLines(qc_lines, file.path(outdir, paste0(model_name, "_diagnostics_QC_summary.txt")))
  add_note("QC post-processing completed for model: ", model_name)
}

###############################################################################
# Metadata, design, canonical-correlation, and MDS diagnostics
###############################################################################

write_table <- function(tab, path) {
  safe_write_csv(as.data.frame.matrix(tab), path)
}

model_required_vars <- function(spec) {
  vars <- c("Population", "Sequencing_batch")
  if (!is.null(spec$grouping_var)) vars <- c(vars, spec$grouping_var)
  if (isTRUE(spec$climate)) vars <- c(vars, "Climate_PC1", "Climate_PC2", "Climate_PC3")
  if (isTRUE(spec$combined)) vars <- c(vars, "ECe_mean")
  unique(vars)
}

run_expression_diagnostics <- function() {
  if (!file.exists(gse_rds)) {
    add_note("GSE RDS missing; expression diagnostics skipped: ", gse_rds)
    return(invisible(NULL))
  }

  gse <- tryCatch(readRDS(gse_rds), error = function(e) {
    add_note("GSE READ ERROR: ", conditionMessage(e)); NULL
  })
  if (is.null(gse)) return(invisible(NULL))

  if (!"counts" %in% assayNames(gse)) {
    add_note("GSE does not contain a counts assay; expression diagnostics skipped.")
    return(invisible(NULL))
  }

  cts <- as.matrix(assay(gse, "counts"))
  meta <- as.data.frame(colData(gse), stringsAsFactors = FALSE)
  sample_ids <- colnames(cts)
  meta_ids <- rownames(meta)

  if (is.null(sample_ids)) {
    add_note("Count matrix has no column names; expression diagnostics skipped.")
    return(invisible(NULL))
  }
  if (is.null(meta_ids)) meta_ids <- sample_ids

  if (all(sample_ids %in% meta_ids)) {
    meta <- meta[match(sample_ids, meta_ids), , drop = FALSE]
  } else if (nrow(meta) == ncol(cts)) {
    add_note("Metadata rownames did not match count columns; positional alignment used because dimensions match.")
  } else {
    add_note("Count/metadata sample IDs cannot be aligned; expression diagnostics skipped.")
    return(invisible(NULL))
  }
  rownames(meta) <- sample_ids

  factor_vars <- c("Population", "Sequencing_batch", "Site_of_distribution_range",
                   "Site_of_distribution_range_East_West")
  for (v in intersect(factor_vars, colnames(meta))) meta[[v]] <- factor(meta[[v]])

  climate_vars <- intersect(c("Climate_PC1", "Climate_PC2", "Climate_PC3"), colnames(meta))
  for (v in climate_vars) meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))

  # Add ECe_mean as a numeric covariate
  if ("ECe_mean" %in% colnames(meta)) {
    meta$ECe_mean <- suppressWarnings(as.numeric(meta$ECe_mean))
  }

  core_vars <- c("Population", "Sequencing_batch",
               "Climate_PC1", "Climate_PC2", "Climate_PC3",
               "ECe_mean")
  core_present <- intersect(core_vars, colnames(meta))
  missing_core <- setdiff(core_vars, colnames(meta))
  if (length(missing_core) > 0) add_note("Missing global metadata variables: ", paste(missing_core, collapse = ", "))

  missingness <- data.frame(
    Variable = core_present,
    N_missing = vapply(meta[core_present], function(x) sum(is.na(x)), integer(1)),
    stringsAsFactors = FALSE
  )
  safe_write_csv(missingness, file.path(correlation_root, "metadata_missingness_core.csv"))

  if (all(c("Population", "Sequencing_batch") %in% colnames(meta))) {
    write_table(table(meta$Population, meta$Sequencing_batch),
                file.path(correlation_root, "table_population_by_batch.csv"))
  }
  
  if (length(climate_vars) == 3) {
  complete <- complete.cases(meta[, climate_vars, drop = FALSE])
  cor_mat <- cor(meta[complete, climate_vars, drop = FALSE], method = "pearson")
  safe_write_csv(data.frame(Variable = rownames(cor_mat), cor_mat, check.names = FALSE),
                 file.path(correlation_root, "climate_pc_correlation.csv"))
  writeLines(capture.output(summary(meta[complete, climate_vars, drop = FALSE])),
             file.path(correlation_root, "climate_pc_summary.txt"))
  add_note("Climate correlation used ", sum(complete), " complete samples.")

  # Add ECe_mean to correlation if present
  if ("ECe_mean" %in% colnames(meta)) {
    vars_with_ece <- c(climate_vars, "ECe_mean")
    complete_ece <- complete.cases(meta[, vars_with_ece, drop = FALSE])
    cor_mat_ece <- cor(meta[complete_ece, vars_with_ece, drop = FALSE], method = "pearson")
    safe_write_csv(data.frame(Variable = rownames(cor_mat_ece), cor_mat_ece, check.names = FALSE),
                   file.path(correlation_root, "climate_pc_and_ECe_correlation.csv"))
    add_note("Climate + ECe correlation used ", sum(complete_ece), " complete samples.")
  }
}
  
  if (!all(c("Sequencing_batch") %in% colnames(meta))) {
    add_note("Sequencing_batch unavailable; MDS diagnostics skipped.")
    return(invisible(NULL))
  }

  design_filter <- model.matrix(~ Sequencing_batch, data = meta)
  keep <- filterByExpr(cts, design = design_filter)
  if (sum(keep) == 0) {
    add_note("filterByExpr retained zero genes; MDS diagnostics skipped.")
    return(invisible(NULL))
  }

  dge <- calcNormFactors(DGEList(counts = cts[keep, , drop = FALSE]))
  
  voom_formula <- if (all(c(
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3",
    "Site_of_distribution_range"
  ) %in% colnames(meta))) {

    ~ Climate_PC1 + Climate_PC2 + Climate_PC3 + Site_of_distribution_range + (1 | Sequencing_batch) + (1 | Population)

  } else if (all(c(
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3"
  ) %in% colnames(meta))) {

  ~ Climate_PC1 + Climate_PC2 + Climate_PC3 + (1 | Sequencing_batch) + (1 | Population)

  } else {
    ~ 1 + (1 | Sequencing_batch) + (1 | Population)
  }

  v <- tryCatch(
    voomWithDreamWeights(dge, formula = voom_formula, data = meta, plot = FALSE),
    error = function(e) { add_note("voomWithDreamWeights ERROR: ", conditionMessage(e)); NULL }
  )
  if (is.null(v) || is.null(v$E)) return(invisible(NULL))

  mds <- limma::plotMDS(v$E, top = min(mds_top_genes, nrow(v$E)), plot = FALSE)
  mds_df <- data.frame(Sample = colnames(v$E), Dim1 = mds$x, Dim2 = mds$y,
                       Population = meta$Population,
                       Sequencing_batch = meta$Sequencing_batch,
                       stringsAsFactors = FALSE)
  for (vname in c("Site_of_distribution_range", "Site_of_distribution_range_East_West")) {
    mds_df[[vname]] <- if (vname %in% colnames(meta)) meta[[vname]] else NA
  }
  safe_write_csv(mds_df, file.path(correlation_root, "mds_coordinates.csv"))

  for (m in model_names) {
    spec <- model_specs[[m]]
    outdir <- file.path(correlation_root, m)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    report <- c(paste0("Diagnostic report: ", m), paste0("Label: ", spec$label), "")

    if (!isTRUE(spec$enabled)) {
      report <- c(report, "Status: skipped because this model is disabled.")
      writeLines(report, file.path(outdir, "diagnostic_summary_report.txt"))
      next
    }

    req <- model_required_vars(spec)
    missing <- setdiff(req, colnames(meta))
    if (length(missing) > 0) {
      report <- c(report, "Status: skipped because required variables are missing.",
                  paste0("Missing: ", paste(missing, collapse = ", ")))
      writeLines(report, file.path(outdir, "diagnostic_summary_report.txt"))
      add_note("Model diagnostics skipped for ", m, ": missing metadata.")
      next
    }

    model_meta <- meta
    complete <- complete.cases(model_meta[, req, drop = FALSE])
    report <- c(report, paste0("Complete cases for required variables: ", sum(complete), " / ", nrow(meta)), "")
    safe_write_csv(data.frame(Variable = req,
                              N_missing = vapply(model_meta[req], function(x) sum(is.na(x)), integer(1))),
                   file.path(outdir, "metadata_missingness.csv"))

    write_table(table(model_meta$Population, model_meta$Sequencing_batch),
                file.path(outdir, "table_population_by_batch.csv"))

    if (!is.null(spec$grouping_var)) {
      g <- spec$grouping_var
      write_table(table(model_meta[[g]], model_meta$Sequencing_batch),
                  file.path(outdir, paste0("table_", g, "_by_batch.csv")))
      write_table(table(model_meta[[g]], model_meta$Population),
                  file.path(outdir, paste0("table_", g, "_by_population.csv")))
    }

    if (isTRUE(spec$climate) && length(climate_vars) == 3) {
      safe_write_csv(data.frame(Variable = rownames(cor_mat), cor_mat, check.names = FALSE),
                     file.path(outdir, "climate_pc_correlation.csv"))
      writeLines(capture.output(summary(meta[complete.cases(meta[, climate_vars, drop = FALSE]), climate_vars, drop = FALSE])),
                 file.path(outdir, "climate_pc_summary.txt"))
    }

    # canCorPairs() is wrapped because model-specific missingness or sparse
    # factor levels can make a particular diagnostic unavailable.
	  if (all(complete)) {
      terms <- if (is.null(spec$grouping_var)) character() else spec$grouping_var
      if (isTRUE(spec$climate)) terms <- c(terms, "Climate_PC1", "Climate_PC2", "Climate_PC3")
	  if (isTRUE(spec$combined)) terms <- c(terms, "ECe_mean")
      fixed <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
      form <- as.formula(paste0("~ ", fixed, " + (1 | Sequencing_batch) + (1 | Population)"))
      cc <- tryCatch(canCorPairs(form, model_meta), error = function(e) {
        report <<- c(report, paste0("canCorPairs unavailable: ", conditionMessage(e)))
        NULL
      })
      if (!is.null(cc)) {
        safe_write_csv(data.frame(Variable = rownames(cc), as.data.frame(cc), check.names = FALSE),
                       file.path(outdir, paste0("canCorPairs_", m, ".csv")))
        report <- c(report, "canCorPairs diagnostic: completed.")
      }
    } else {
      report <- c(report, "canCorPairs diagnostic: skipped because required variables contain missing values.")
    }

    p1 <- ggplot(mds_df, aes(Dim1, Dim2, color = Population, shape = Sequencing_batch)) +
      geom_point(size = 3, alpha = 0.9) + theme_bw() +
      labs(title = paste("MDS -", spec$label), color = "Population", shape = "Batch")
    ggsave(file.path(outdir, paste0("mds_", m, "_population_batch.png")), p1, width = 8, height = 6, dpi = 300)

    if (!is.null(spec$grouping_var) && spec$grouping_var %in% colnames(mds_df)) {
      p2 <- ggplot(mds_df, aes(Dim1, Dim2, color = .data[[spec$grouping_var]], shape = Sequencing_batch)) +
        geom_point(size = 3, alpha = 0.9) + theme_bw() +
        labs(title = paste("MDS -", spec$label), color = spec$grouping_var, shape = "Batch")
      ggsave(file.path(outdir, paste0("mds_", m, "_grouping_batch.png")), p2, width = 8, height = 6, dpi = 300)
    }

    report <- c(report,
      "Interpretation:",
      "- Contingency tables describe overlap and possible confounding; they do not test causality.",
      "- canCorPairs() is a design-association diagnostic, not a DE significance test.",
      "- MDS plots show expression-space structure and should be interpreted with metadata and model results.")
    writeLines(report, file.path(outdir, "diagnostic_summary_report.txt"))
    add_note("Expression diagnostics completed for model: ", m)
  }

  status <- data.frame(
    Model = model_names,
    Enabled = vapply(model_names, function(x) isTRUE(model_specs[[x]]$enabled), logical(1)),
    Required_variables = vapply(model_names, function(x) paste(model_required_vars(model_specs[[x]]), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  safe_write_csv(status, file.path(correlation_root, "diagnostic_models_summary.csv"))
  writeLines(capture.output(status), file.path(correlation_root, "diagnostic_models_summary.txt"))
}

###############################################################################
# Execute all components
###############################################################################

add_note("Starting merged Anastatica post-processing script.")
add_note("Primary DEG files will not be replaced.")

run_deg_summaries()
run_varpart_summaries()
for (m in model_names) run_model_qc(m)
run_expression_diagnostics()

notes_path <- file.path(degs_summary_root, "merged_postprocessing_run_notes.txt")
writeLines(c(
  "Merged Anastatica post-processing run notes",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("degs_root: ", degs_root),
  paste0("degs_summary_root: ", degs_summary_root),
  paste0("correlation_root: ", correlation_root),
  paste0("FDR threshold: ", fdr_thresh),
  paste0("Absolute log2FC threshold: ", signif(lfc_thresh, 8)),
  paste0("Exploratory z-score threshold: ", z_threshold),
  "",
  run_notes
), notes_path)

message("Merged post-processing completed.")
message("Summary outputs: ", degs_summary_root)
message("Correlation/QC outputs: ", correlation_root)