#!/usr/bin/env Rscript

# ==============================================================================
# Date: 28/08/2026
# Author: Iris Sammarco
# ==============================================================================
#
# Purpose:
# --------
# Differential expression analysis of Anastatica greenhouse-grown populations
# using dream (variancePartition) with mixed-effects models, explicit pairwise
# contrasts for multi-level factors, and variance partitioning.
#
# This script contains ONLY single-predictor models:
#   1) climate_model
#   2) distribution_range_model (3-level: East, Centre, West)
#   3) distribution_range_east_west_model (2-level: East, West)
#   4) ECe_model
#
# Scientific rationale:
# ---------------------
# - Uses mixed-effects models with dream for differential expression, retaining
#   Population as a random intercept to account for non-independence among
#   individuals from the same population.
# - Models Sequencing_batch as a random effect in all DE and variance
#   partitioning models.
# - Uses fitExtractVarPartModel() for variance partitioning to quantify the
#   fraction of expression variance attributable to model terms.
# - Population and Sequencing_batch are modeled as random effects, following
#   variancePartition guidance that categorical grouping variables are
#   typically treated as random for variance-fraction estimation.
# - dream() and fitExtractVarPartModel() serve different purposes:
#     * dream(): tests differential expression and estimates fixed/random effects.
#     * fitExtractVarPartModel(): decomposes variance into contributions from
#       fixed and random terms.
#
# For Dream and Mixed Models in general:
#   being a fixed effect does NOT mean it is the predictor of interest.
#
# Variance partitioning:
#   random effect -> How much variability exists among groups?
#   fixed effect  -> How much higher/lower is one level relative to reference?
#
# Models implemented in this script:
# ----------------------------------
# 1) climate_model (ACTIVE):
#    DE and Variance partitioning formulas:
#      ~ Climate_PC1 + Climate_PC2 + Climate_PC3 +
#        (1 | Sequencing_batch) + (1 | Population)
#    Tests the independent effects of the three climate PCs.
#
# 2) distribution_range_model (ACTIVE):
#    3-level distribution range model:
#      ~ Site_of_distribution_range +
#        (1 | Sequencing_batch) + (1 | Population)
#    Levels: East, Centre/Center, West
#    Pairwise contrasts: East vs Centre, West vs Centre, East vs West
#
# 3) distribution_range_east_west_model (ACTIVE):
#    2-level distribution range model:
#      ~ Site_of_distribution_range_East_West +
#        (1 | Sequencing_batch) + (1 | Population)
#    Centre populations are coded as East.
#    Tests: East vs West
#
# 4) ECe_model (ACTIVE):
#      ~ ECe_mean +
#        (1 | Sequencing_batch) + (1 | Population)
#
# From the variancePartition documentation:
# -----------------------------------------
# - Categorical variables should (almost) always be modeled as a random effect
#   for variance partitioning.
# - dream() is designed for differential expression testing with mixed models,
#   and fitExtractVarPartModel() is designed to estimate variance fractions.
#
# In the summary metrics:
# -----------------------
# For each dream model result:
# - Extracts the per-gene lmer fits from the dream object.
# - Computes per-gene AIC.
# - Saves: _AIC_per_gene.csv, _AIC_summary.csv
# - Combines all model-level summaries into one file: AIC_summary_all_models.csv
#
# ==============================================================================
# Environment setup used for this run
# ==============================================================================
# 1) Start an interactive job with enough resources:
#    qsub -I -l select=1:ncpus=2:mem=40gb -l walltime=4:00:00
#
# 2) Create the conda environment with the R version used for the analysis:
#    mamba create -n dream_degs -c conda-forge r-base=4.4 r-biocmanager r-ggplot2 r-reformulas -y
#
# 3) Install system-style libraries and build tools needed by R/Bioconductor:
#    mamba activate dream_degs
#    mamba install -c conda-forge zlib libxml2 libcurl libpng make gxx_linux-64 -y
#
# 4) Reopen the session with the conda environment active and no module R loaded:
#    mamba deactivate
#    module purge
#    mamba activate dream_degs
#    R
#
# 5) Install Bioconductor packages inside the R 4.5.3 conda session:
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager", repos = "https://cloud.r-project.org")
#    BiocManager::install(version = "3.22", ask = FALSE)
#    BiocManager::install(c("GenomeInfoDbData", "GenomeInfoDb", "SummarizedExperiment",
#      "edgeR", "limma", "variancePartition", "reformulas", "BiocParallel"), ask = FALSE)
#
# Notes:
# ------
# - GenomeInfoDbData is installed through BiocManager, not via r-genomeinfodbdata.
# - Variance partitioning is run with fitExtractVarPartModel(), modeling
#   categorical variables as random effects and continuous covariates as fixed.
# ==============================================================================

# ==============================================================================
# Library paths and loading
# ==============================================================================
message("R version: ", R.version.string)
message("variancePartition: ", as.character(packageVersion("variancePartition")))
message("lme4: ", as.character(packageVersion("lme4")))
message("lmerTest: ", as.character(packageVersion("lmerTest")))
message("Matrix: ", as.character(packageVersion("Matrix")))
message("zenith: ", if (requireNamespace("zenith", quietly = TRUE)) as.character(packageVersion("zenith")) else "not installed")

suppressPackageStartupMessages({
  library(reformulas)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(lme4)
  library(BiocParallel)
  library(ggplot2)
})


# ==============================================================================
# Configuration
# ==============================================================================
tximeta_outdir <- "/Anastatica/tximeta"
gse_rds <- file.path(tximeta_outdir, "tximeta_gene_level.rds")
degs_root <- "/Anastatica/degs"
fdr_thresh <- 0.05
lfc_thresh <- log2(1.5)
bp_param <- SerialParam()

model_names <- c("climate_model", "distribution_range_model", "distribution_range_east_west_model", "ECe_model")
model_dirs <- setNames(file.path(degs_root, model_names), model_names)
for (d in model_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# Load data
# ==============================================================================
if (!file.exists(gse_rds)) stop("Gene-level RDS not found at: ", gse_rds)
gse <- readRDS(gse_rds)
cts <- as.matrix(assay(gse, "counts"))
metadata <- as.data.frame(colData(gse), stringsAsFactors = FALSE)
match_col <- "quant_dir_name"
if (!match_col %in% colnames(metadata)) stop("metadata must contain a '", match_col, "' column.")
if (!all(colnames(cts) %in% metadata[[match_col]])) stop("The following count columns are missing from metadata: ", paste(setdiff(colnames(cts), metadata[[match_col]]), collapse = ", "))
rownames(metadata) <- metadata[[match_col]]
metadata <- metadata[colnames(cts), , drop = FALSE]


# ==============================================================================
# Metadata typing and QC
# ==============================================================================
prepare_metadata <- function(meta) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  required_vars <- c("Sample_name", "Population", "Sequencing_batch", "Site_of_distribution_range", "Site_of_distribution_range_East_West", "Climate_PC1", "Climate_PC2", "Climate_PC3", "ECe_mean")
  missing_vars <- setdiff(required_vars, colnames(meta))
  if (length(missing_vars) > 0) stop("Missing required metadata columns: ", paste(missing_vars, collapse = ", "))
  meta$Population <- factor(meta$Population)
  meta$Sequencing_batch <- factor(meta$Sequencing_batch)
  meta$Site_of_distribution_range <- factor(meta$Site_of_distribution_range)
  meta$Site_of_distribution_range_East_West <- factor(meta$Site_of_distribution_range_East_West)
  cont_vars <- c("Climate_PC1", "Climate_PC2", "Climate_PC3", "ECe_mean")
  for (v in cont_vars) meta[[v]] <- as.numeric(meta[[v]])
  if (anyNA(meta[, c("Population", "Sequencing_batch", "Site_of_distribution_range", "Site_of_distribution_range_East_West", "Climate_PC1", "Climate_PC2", "Climate_PC3", "ECe_mean")])) stop("Metadata contains NA values in variables used for modelling.")
  meta
}
metadata <- prepare_metadata(metadata)


# ==============================================================================
# Set reference levels
# ==============================================================================
message("Distribution range levels before releveling:"); print(levels(metadata$Site_of_distribution_range))
if ("Centre" %in% levels(metadata$Site_of_distribution_range)) {
  metadata$Site_of_distribution_range <- relevel(metadata$Site_of_distribution_range, ref = "Centre")
} else if ("Center" %in% levels(metadata$Site_of_distribution_range)) {
  metadata$Site_of_distribution_range <- relevel(metadata$Site_of_distribution_range, ref = "Center")
} else stop("Could not find a Centre/Center level in Site_of_distribution_range. Levels are: ", paste(levels(metadata$Site_of_distribution_range), collapse = ", "))
message("Distribution range levels after releveling:"); print(levels(metadata$Site_of_distribution_range))

if (!all(c("East", "West") %in% levels(metadata$Site_of_distribution_range_East_West)) || length(levels(metadata$Site_of_distribution_range_East_West)) != 2) stop("Site_of_distribution_range_East_West must contain exactly East and West. Found: ", paste(levels(metadata$Site_of_distribution_range_East_West), collapse = ", "))
metadata$Site_of_distribution_range_East_West <- relevel(metadata$Site_of_distribution_range_East_West, ref = "West")


# ==============================================================================
# Gene filtering and normalization
# ==============================================================================
design_filter <- model.matrix(~ Sequencing_batch, metadata)
keep <- filterByExpr(cts, design = design_filter)
cts_filt <- cts[keep, , drop = FALSE]
dge <- DGEList(counts = cts_filt)
dge <- calcNormFactors(dge)


# ==============================================================================
# Helper functions
# ==============================================================================
save_voom_summary <- function(vobj, outfile, model_name) {
  capture.output({
    cat("Voom object summary for model:", model_name, "\n\n")
    str(vobj, max.level = 1)
    cat("\nDimensions of expression matrix:", dim(vobj$E), "\n")
    cat("Dimensions of weights matrix:", dim(vobj$weights), "\n")
  }, file = outfile)
}

add_gene_id_column <- function(df) {
  df$FeatureID <- rownames(df)
  df <- df[, c("FeatureID", setdiff(colnames(df), "FeatureID")), drop = FALSE]
  rownames(df) <- NULL
  df
}

write_de_sets <- function(df, prefix, outdir, fdr_thresh, lfc_thresh, apply_lfc = TRUE) {
  df <- add_gene_id_column(df)
  write.csv(df, file = file.path(outdir, paste0(prefix, "_full.csv")), row.names = FALSE)
  df_fdr <- subset(df, adj.P.Val <= fdr_thresh)
  write.csv(df_fdr, file = file.path(outdir, paste0(prefix, "_FDR05.csv")), row.names = FALSE)
  if (apply_lfc && "logFC" %in% colnames(df)) {
    df_fdr_lfc <- subset(df, adj.P.Val <= fdr_thresh & abs(logFC) >= lfc_thresh)
    write.csv(df_fdr_lfc, file = file.path(outdir, paste0(prefix, "_FDR05_log2FC_1.5.csv")), row.names = FALSE)
  }
}

summarize_de_results <- function(df, model_name) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  required_cols <- c("Test_type", "Test", "adj.P.Val")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) stop("Result table is missing required columns: ", paste(missing_cols, collapse = ", "))
  test_groups <- unique(df[, c("Test_type", "Test"), drop = FALSE])
  summary_list <- lapply(seq_len(nrow(test_groups)), function(i) {
    test_type_i <- test_groups$Test_type[i]
    test_label_i <- test_groups$Test[i]
    test_df <- df[df$Test_type == test_type_i & df$Test == test_label_i, , drop = FALSE]
    tested <- !is.na(test_df$adj.P.Val)
    sig_lfc <- if ("logFC" %in% colnames(test_df)) sum(tested & test_df$adj.P.Val <= fdr_thresh & !is.na(test_df$logFC) & abs(test_df$logFC) >= lfc_thresh) else NA_integer_
    data.frame(Model = model_name, Test_type = as.character(test_type_i), Test_label = as.character(test_label_i), N_genes_tested = sum(tested), N_sig_FDR05 = sum(tested & test_df$adj.P.Val <= fdr_thresh), N_sig_FDR05_log2FC_1.5 = sig_lfc, stringsAsFactors = FALSE)
  })
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  summary_df
}

extract_coefficient_results <- function(fit, coefs_to_test) {
  tt_all <- lapply(coefs_to_test, function(cf) {
    tt <- variancePartition::topTable(fit, coef = cf, number = Inf, sort.by = "P")
    tt$Test <- cf
    tt$Test_type <- "coefficient"
    tt
  })
  names(tt_all) <- coefs_to_test
  res_df <- do.call(rbind, tt_all)
  res_df <- res_df[order(res_df$Test, res_df$P.Value), , drop = FALSE]
  res_df
}

extract_omnibus_results <- function(fit, coef_vector, test_name) {
  tt <- variancePartition::topTable(fit, coef = coef_vector, number = Inf, sort.by = "P")
  tt$Test <- test_name
  tt$Test_type <- "omnibus_Ftest"
  tt
}

run_model <- function(
  model_name,
  form,
  coefficient_terms = NULL,
  omnibus_coef = NULL,
  omnibus_name = NULL,
  contrast_form = NULL,
  contrast_matrix = NULL,
  contrast_names = NULL,
  dge_obj = dge,
  meta = metadata) {
  
  message("============================================")
  message("Running model: ", model_name)
  message("============================================")
  
  outdir <- model_dirs[[model_name]]
  
  write.csv(
    meta,
    file = file.path(
      outdir,
      paste0(model_name, "_sample_metadata_used.csv")
    ),
    row.names = FALSE
  )
  
  fixed_form <- reformulas::nobars(form)
  design_mat <- model.matrix(fixed_form, meta)
  
  write.csv(
    design_mat,
    file = file.path(
      outdir,
      paste0(model_name, "_design_matrix.csv")
    ),
    row.names = TRUE
  )
  
  qr_design <- qr(design_mat)
  
  if (qr_design$rank < ncol(design_mat)) {
    dropped_columns <- colnames(design_mat)[
      qr_design$pivot[
        seq.int(qr_design$rank + 1, ncol(design_mat))
      ]
    ]
    
    stop(
      "Rank-deficient fixed-effect design for model '",
      model_name,
      "'.\nRank = ",
      qr_design$rank,
      "; columns = ",
      ncol(design_mat),
      ".\nPotentially redundant columns:\n",
      paste(dropped_columns, collapse = "\n")
    )
  }
  
  # --------------------------------------------------------------------------
  # Voom
  # --------------------------------------------------------------------------
  
  vobj <- voomWithDreamWeights(
    dge_obj,
    formula = form,
    data = meta,
    plot = FALSE,
    BPPARAM = bp_param
  )
  
  weights_mat <- vobj$weights
  
  if (!is.null(weights_mat) &&
      is.null(colnames(weights_mat))) {
    colnames(weights_mat) <- colnames(dge_obj$counts)
  }
  
  save_voom_summary(
    vobj = vobj,
    outfile = file.path(
      outdir,
      paste0(model_name, "_voom_weights_summary.txt")
    ),
    model_name = model_name
  )
  
  if (is.null(vobj$E) || nrow(vobj$E) == 0) {
    stop(
      "voomWithDreamWeights returned an empty expression matrix for model: ",
      model_name
    )
  }
  
  # --------------------------------------------------------------------------
  # Base dream model
  # --------------------------------------------------------------------------
  
  old_warn <- getOption("warn")
  options(warn = 1)
  
  fit_base <- suppressWarnings(
    dream(
      vobj,
      form,
      meta,
      BPPARAM = bp_param,
      computeResiduals = TRUE,
      ddf = "Satterthwaite",
      hideErrorsInBackend = FALSE
    )
  )
  
  options(warn = old_warn)
  
  fit_base <- variancePartition::eBayes(fit_base)
  
  saveRDS(
    fit_base,
    file = file.path(
      outdir,
      paste0(model_name, "_fit_base_ebayes.rds")
    )
  )
  
  saveRDS(
    vobj,
    file = file.path(
      outdir,
      paste0(model_name, "_vobj.rds")
    )
  )
  
  message(
    "Saved base eBayes fit: ",
    file.path(
      outdir,
      paste0(model_name, "_fit_base_ebayes.rds")
    )
  )
  
  # --------------------------------------------------------------------------
  # Model coefficients
  # --------------------------------------------------------------------------
  
  coef_df <- as.data.frame(coef(fit_base))
  coef_df <- add_gene_id_column(coef_df)
  
  write.csv(
    coef_df,
    file = file.path(
      outdir,
      paste0(model_name, "_model_coefficients.csv")
    ),
    row.names = FALSE
  )
  
  available_tests <- colnames(coef(fit_base))
  
  fixed_coefs <- setdiff(
    available_tests,
    "(Intercept)"
  )
  
  if (!is.null(coefficient_terms)) {
    fixed_coefs <- intersect(
      fixed_coefs,
      coefficient_terms
    )
  }
  
  # --------------------------------------------------------------------------
  # Coefficient tests
  # --------------------------------------------------------------------------
  
  de_full <- NULL
  
  if (length(fixed_coefs) > 0) {
    
    de_full <- extract_coefficient_results(
      fit_base,
      fixed_coefs
    )
    
    write_de_sets(
      df = de_full,
      prefix = paste0(model_name, "_DE_results"),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh
    )
  }
  
  # --------------------------------------------------------------------------
  # Omnibus F-test
  # --------------------------------------------------------------------------
  
  omnibus_full <- NULL
  
  if (!is.null(omnibus_coef) &&
      length(omnibus_coef) > 0) {
    
    available_base_tests <- colnames(coef(fit_base))
    
    valid_omnibus <- intersect(
      omnibus_coef,
      available_base_tests
    )
    
    if (length(valid_omnibus) == 0) {
      stop(
        "None of the omnibus coefficients were found in the BASE fitted model: ",
        model_name
      )
    }
    
    omnibus_full <- extract_omnibus_results(
      fit = fit_base,
      coef_vector = valid_omnibus,
      test_name = omnibus_name
    )
    
    write_de_sets(
      df = omnibus_full,
      prefix = paste0(model_name, "_omnibus_Ftest"),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh,
      apply_lfc = FALSE
    )
  }
  
  # --------------------------------------------------------------------------
  # Pairwise contrasts
  # --------------------------------------------------------------------------
  
  pairwise_full <- NULL
  fit_contrast <- NULL
  
  if (!is.null(contrast_form) &&
      !is.null(contrast_matrix) &&
      ncol(contrast_matrix) > 0) {
    
    message(
      "Running pairwise contrasts for model: ",
      model_name
    )
    
    L_contrast <- contrast_matrix
    
    if (is.null(colnames(L_contrast))) {
      stop(
        "contrast_matrix must have column names for model: ",
        model_name
      )
    }
    
    message("Contrast matrix:")
    print(L_contrast)
    
    fit_contrast <- suppressWarnings(
      variancePartition::dream(
        vobj,
        contrast_form,
        meta,
        L = L_contrast,
        BPPARAM = bp_param,
        computeResiduals = FALSE,
        ddf = "Satterthwaite",
        hideErrorsInBackend = FALSE
      )
    )
    
    fit_contrast <- variancePartition::eBayes(
      fit_contrast
    )
    
    saveRDS(
      fit_contrast,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_fit_contrast_ebayes.rds"
        )
      )
    )
    
    valid_contrasts <- intersect(
      colnames(coef(fit_contrast)),
      colnames(L_contrast)
    )
    
    if (length(valid_contrasts) == 0) {
      stop(
        "No requested contrasts were found in the contrast fit for model: ",
        model_name
      )
    }
    
    pairwise_list <- lapply(
      valid_contrasts,
      function(cf) {
        
        tt <- variancePartition::topTable(
          fit_contrast,
          coef = cf,
          number = Inf,
          sort.by = "P"
        )
        
        tt$Test <- cf
        tt$Test_type <- "contrast"
        
        tt
      }
    )
    
    pairwise_full <- do.call(
      rbind,
      pairwise_list
    )
    
    pairwise_full <- pairwise_full[
      order(
        pairwise_full$Test,
        pairwise_full$P.Value
      ),
      ,
      drop = FALSE
    ]
    
    write_de_sets(
      df = pairwise_full,
      prefix = paste0(
        model_name,
        "_pairwise_contrasts"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh
    )
  }
  
  expr_log2 <- vobj$E
  form_vp <- if (model_name == "climate_model") ~ Climate_PC1 + Climate_PC2 + Climate_PC3 + (1 | Sequencing_batch) + (1 | Population) else if (model_name == "distribution_range_model") ~ (1 | Site_of_distribution_range) + (1 | Sequencing_batch) + (1 | Population) else if (model_name == "distribution_range_east_west_model") ~ (1 | Site_of_distribution_range_East_West) + (1 | Sequencing_batch) + (1 | Population) else if (model_name == "ECe_model") ~ ECe_mean + (1 | Sequencing_batch) + (1 | Population) else stop("No variance-partition formula defined for model: ", model_name)
  varPart <- fitExtractVarPartModel(expr_log2, form_vp, meta, BPPARAM = bp_param)
  write.csv(as.data.frame(varPart), file = file.path(outdir, paste0(model_name, "_varPart_results.csv")), row.names = TRUE)
  vp_summary <- apply(varPart, 2, function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)))
  capture.output({ cat("Variance fraction summary for model:", model_name, "\n\n"); print(round(vp_summary, 4)) }, file = file.path(outdir, paste0(model_name, "_varPart_summary.txt")))
  pdf(file.path(outdir, paste0(model_name, "_varPart_plot.pdf")), width = 7, height = 5); print(plotVarPart(sortCols(varPart))); dev.off()
  res_mat <- residuals(fit_base, y = vobj)
  res_var <- apply(res_mat, 1, function(x) var(x, na.rm = TRUE))
  median_gene <- names(res_var)[order(res_var)][ceiling(length(res_var) / 2)]
  expr_mat <- vobj$E
  diag_df <- data.frame(Observed = expr_mat[median_gene, ], Residuals = res_mat[median_gene, ])
  pdf(file.path(outdir, paste0(model_name, "_diagnostics_residuals_vs_observed.pdf")), width = 7, height = 5); print(ggplot(diag_df, aes(x = Observed, y = Residuals)) + geom_point(alpha = 0.6) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + theme_bw() + labs(title = paste(model_name, "- residuals vs observed (median residual variance gene)"), x = "Observed log-expression", y = "Residuals")); dev.off()
  diag_df_sqrt <- data.frame(Observed = expr_mat[median_gene, ], SqrtAbsRes = sqrt(abs(res_mat[median_gene, ])))
  pdf(file.path(outdir, paste0(model_name, "_diagnostics_scale_location.pdf")), width = 7, height = 5); print(ggplot(diag_df_sqrt, aes(x = Observed, y = SqrtAbsRes)) + geom_point(alpha = 0.6) + theme_bw() + labs(title = paste(model_name, "- scale-location (sqrt(|residuals|) vs fitted)"), x = "Observed log-expression", y = "sqrt(|Residuals|)")); dev.off()
  qq_df <- data.frame(Sample = sort(diag_df$Residuals), Theoretical = sort(qnorm(ppoints(length(diag_df$Residuals)))))
  pdf(file.path(outdir, paste0(model_name, "_diagnostics_residuals_QQ.pdf")), width = 7, height = 5); print(ggplot(qq_df, aes(x = Theoretical, y = Sample)) + geom_point(alpha = 0.6) + geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + theme_bw() + labs(title = paste(model_name, "- residuals QQ plot (median residual variance gene)"), x = "Theoretical quantiles", y = "Sample residuals")); dev.off()
  res_summary <- data.frame(Gene = names(res_var), Residual_variance = res_var)
  res_summary$Residual_sd <- sqrt(res_summary$Residual_variance)
  write.csv(res_summary, file = file.path(outdir, paste0(model_name, "_diagnostics_residual_summary.csv")), row.names = FALSE)
  if (!is.null(fit_base$random)) {
    re_example <- fit_base$random[[1]]
    re_list <- lme4::ranef(re_example)
    if ("Population" %in% names(re_list)) { pop_re <- re_list[["Population"]]; pop_re$Level <- rownames(pop_re); colnames(pop_re) <- paste0("RE_", colnames(pop_re)); write.csv(pop_re, file = file.path(outdir, paste0(model_name, "_diagnostics_random_effects_Population_example.csv")), row.names = FALSE) }
    if ("Sequencing_batch" %in% names(re_list)) { batch_re <- re_list[["Sequencing_batch"]]; batch_re$Level <- rownames(batch_re); colnames(batch_re) <- paste0("RE_", colnames(batch_re)); write.csv(batch_re, file = file.path(outdir, paste0(model_name, "_diagnostics_random_effects_Sequencing_batch_example.csv")), row.names = FALSE) }
    pop_re_var <- sapply(seq_along(fit_base$random), function(i) { re_i <- fit_base$random[[i]]; if (is.null(re_i)) return(NA_real_); re_list_i <- lme4::ranef(re_i); if (!"Population" %in% names(re_list_i)) return(NA_real_); pop_vals <- re_list_i[["Population"]][, 1]; if (length(pop_vals) < 2) return(NA_real_); var(pop_vals, na.rm = TRUE) })
    pop_re_var_df <- data.frame(Gene = names(fit_base$random), Population_RE_variance = pop_re_var)
    write.csv(pop_re_var_df, file = file.path(outdir, paste0(model_name, "_diagnostics_Population_RE_var_per_gene.csv")), row.names = FALSE)
    pdf(file.path(outdir, paste0(model_name, "_diagnostics_Population_RE_var_density.pdf")), width = 7, height = 5); print(ggplot(subset(pop_re_var_df, !is.na(Population_RE_variance)), aes(x = Population_RE_variance)) + geom_density(fill = "steelblue", alpha = 0.4) + theme_bw() + ggtitle(paste(model_name, "- Population random-effect variance across genes")) + xlab("Variance of Population random intercepts") + ylab("Density")); dev.off()
  }
  weights_mat <- vobj$weights
  
  message("Diagnostics for ", model_name, ":")
  message(" dim(vobj$E): ", if (is.null(vobj$E)) "NULL" else paste(dim(vobj$E), collapse = " x "))
  message(" dim(vobj$weights): ", if (is.null(weights_mat)) "NULL" else paste(dim(weights_mat), collapse = " x "))
  message(" nrow(dge_obj$counts): ", nrow(dge_obj$counts))
  message(" ncol(dge_obj$counts): ", ncol(dge_obj$counts))
  message(" length(colnames(weights_mat)): ", if (is.null(colnames(weights_mat))) 0 else length(colnames(weights_mat)))
  
  if (!is.null(weights_mat) && is.matrix(weights_mat) && nrow(weights_mat) > 0 && ncol(weights_mat) > 0 && ncol(weights_mat) == ncol(dge_obj$counts)) {
    mean_w <- apply(weights_mat, 2, mean, na.rm = TRUE)
    sd_w <- apply(weights_mat, 2, sd, na.rm = TRUE)
    sample_names <- colnames(weights_mat)
    if (is.null(sample_names) || length(sample_names) == 0) sample_names <- colnames(dge_obj$counts)
    if (length(sample_names) != length(mean_w)) stop("Cannot build voom weights summary: sample name length (", length(sample_names), ") does not match weights length (", length(mean_w), ").")
    weights_summary <- data.frame(Sample = sample_names, Mean_weight = as.numeric(mean_w), SD_weight = as.numeric(sd_w), row.names = NULL, stringsAsFactors = FALSE)
    write.csv(weights_summary, file = file.path(outdir, paste0(model_name, "_diagnostics_voom_weights_summary.csv")), row.names = FALSE)
    weights_long <- data.frame(Weight = as.vector(weights_mat), Sample = rep(sample_names, each = nrow(weights_mat)), row.names = NULL, stringsAsFactors = FALSE)
    pdf(file.path(outdir, paste0(model_name, "_diagnostics_voom_weights_density.pdf")), width = 7, height = 5); print(ggplot(weights_long, aes(x = Weight)) + geom_density() + theme_bw() + ggtitle(paste(model_name, "- voom weight density (all samples)"))); dev.off()
  } else {
    message("Skipping voom-weight diagnostics for ", model_name, " because weights matrix was NULL, empty, or had unexpected dimensions.")
    writeLines(c(paste("Model:", model_name), paste("dim(vobj$E):", if (is.null(vobj$E)) "NULL" else paste(dim(vobj$E), collapse = " x ")), paste("dim(vobj$weights):", if (is.null(weights_mat)) "NULL" else paste(dim(weights_mat), collapse = " x ")), paste("nrow(dge_obj$counts):", nrow(dge_obj$counts)), paste("ncol(dge_obj$counts):", ncol(dge_obj$counts)), "Voom-weight diagnostics skipped because weights matrix was NULL, empty, or malformed."), con = file.path(outdir, paste0(model_name, "_diagnostics_voom_weights_skipped.txt")))
  }
  
summary_tables <- list(
  coefficient = de_full,
  pairwise = pairwise_full,
  omnibus = omnibus_full
)

summary_tables <- summary_tables[
  !vapply(summary_tables, is.null, logical(1))
]

de_summary <- NULL

if (length(summary_tables) > 0) {
  de_summary <- do.call(
    rbind,
    lapply(
      summary_tables,
      summarize_de_results,
      model_name = model_name
    )
  )

  rownames(de_summary) <- NULL
}

write.csv(
  de_summary,
  file = file.path(
    outdir,
    paste0(model_name, "_DE_summary.csv")
  ),
  row.names = FALSE
)

# --------------------------------------------------------------------------
# Return objects needed for downstream analyses
# --------------------------------------------------------------------------

invisible(list(
  model_name = model_name,
  vobj = vobj,
  fit_base = fit_base,
  fit_contrast = fit_contrast,
  de_full = de_full,
  pairwise_full = pairwise_full,
  omnibus_full = omnibus_full,
  de_summary = de_summary
))
}

compute_aic_for_model <- function(model_name, form, vobj, meta, outdir, bpparam = bp_param) {
  message("Computing AIC metrics for model: ", model_name)
  
  aic_fxn <- function(fit) { if (inherits(fit, "lmerMod")) tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_) else if (inherits(fit, "lm")) tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_) else NA_real_ }
  fit_list <- variancePartition::fitVarPartModel(exprObj = vobj, formula = form, data = meta, REML = FALSE, useWeights = TRUE, fxn = aic_fxn, BPPARAM = bpparam)
  
  aic_vec <- unlist(fit_list, use.names = TRUE)

  metrics_df <- data.frame(
    Gene = names(aic_vec),
    AIC = as.numeric(aic_vec),
    stringsAsFactors = FALSE
  )

  metrics_df <- metrics_df[
    is.finite(metrics_df$AIC),
    ,
    drop = FALSE
  ]
  
  summary_df <- data.frame(
    Model = model_name,
    N_genes = nrow(metrics_df),
    AIC_mean = if (nrow(metrics_df) > 0) mean(metrics_df$AIC) else NA_real_,
    AIC_median = if (nrow(metrics_df) > 0) median(metrics_df$AIC) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  write.csv(metrics_df, file = file.path(outdir, paste0(model_name, "_AIC_per_gene.csv")), row.names = FALSE)
  write.csv(summary_df, file = file.path(outdir, paste0(model_name, "_AIC_summary.csv")), row.names = FALSE)
  
  invisible(list(
    fits = fit_list,
    metrics = metrics_df,
    summary = summary_df
  ))
}

# ==============================================================================
# Contrasts for Site_of_distribution_range
# ==============================================================================

range_formula_base <- ~ Site_of_distribution_range +
  (1 | Sequencing_batch) +
  (1 | Population)

range_formula_contrast <- ~ 0 + Site_of_distribution_range +
  (1 | Sequencing_batch) +
  (1 | Population)

range_center_level <- if (
  "Centre" %in% levels(metadata$Site_of_distribution_range)
) {
  "Centre"
} else if (
  "Center" %in% levels(metadata$Site_of_distribution_range)
) {
  "Center"
} else {
  stop(
    "Could not find Centre/Center level in Site_of_distribution_range."
  )
}

range_center_coef <- paste0(
  "Site_of_distribution_range",
  range_center_level
)

message(
  "Range centre coefficient for contrasts: ",
  range_center_coef
)

L_range <- variancePartition::makeContrastsDream(
  range_formula_contrast,
  metadata,
  contrasts = c(
    East_vs_Center =
      paste0(
        "Site_of_distribution_rangeEast - ",
        range_center_coef
      ),

    East_vs_West =
      "Site_of_distribution_rangeEast - Site_of_distribution_rangeWest",

    West_vs_Center =
      paste0(
        "Site_of_distribution_rangeWest - ",
        range_center_coef
      )
  )
)

message("Range contrast matrix:")
print(L_range)

# ==============================================================================
# Run single-predictor models
# ==============================================================================
climate_res <- run_model(model_name = "climate_model", form = ~ Climate_PC1 + Climate_PC2 + Climate_PC3 + (1 | Sequencing_batch) + (1 | Population), coefficient_terms = c("Climate_PC1", "Climate_PC2", "Climate_PC3"), omnibus_coef = c("Climate_PC1", "Climate_PC2", "Climate_PC3"), omnibus_name = "Climate_omnibus", dge_obj = dge, meta = metadata)

distribution_range_res <- run_model(
  model_name = "distribution_range_model",

  form = range_formula_base,

  contrast_form = range_formula_contrast,

  coefficient_terms = NULL,

  contrast_matrix = L_range,

  contrast_names = colnames(L_range),

  omnibus_coef = c(
    "Site_of_distribution_rangeEast",
    "Site_of_distribution_rangeWest"
  ),

  omnibus_name = "Site_of_distribution_range_omnibus",

  dge_obj = dge,
  meta = metadata
)

distribution_range_east_west_res <- run_model(model_name = "distribution_range_east_west_model", form = ~ Site_of_distribution_range_East_West + (1 | Sequencing_batch) + (1 | Population), coefficient_terms = "Site_of_distribution_range_East_WestEast", omnibus_coef = "Site_of_distribution_range_East_WestEast", omnibus_name = "Site_of_distribution_range_East_West_omnibus", dge_obj = dge, meta = metadata)

ECe_res <- run_model(model_name = "ECe_model", form = ~ ECe_mean + (1 | Sequencing_batch) + (1 | Population), coefficient_terms = "ECe_mean", omnibus_coef = "ECe_mean", omnibus_name = "ECe_mean_omnibus", dge_obj = dge, meta = metadata)


# ==============================================================================
# Compute AIC for each model
# ==============================================================================
aic_summaries <- list()

aic_summaries[["climate_model"]] <-
  compute_aic_for_model(
    model_name = "climate_model",
    form = ~ Climate_PC1 + Climate_PC2 + Climate_PC3 +
      (1 | Sequencing_batch) +
      (1 | Population),
    vobj = climate_res$vobj,
    meta = metadata,
    outdir = model_dirs[["climate_model"]],
    bpparam = bp_param
  )

aic_summaries[["ECe_model"]] <-
  compute_aic_for_model(
    model_name = "ECe_model",
    form = ~ ECe_mean +
      (1 | Sequencing_batch) +
      (1 | Population),
    vobj = ECe_res$vobj,
    meta = metadata,
    outdir = model_dirs[["ECe_model"]],
    bpparam = bp_param
  )

# ==============================================================================
# Common voom object for range-model AIC comparison
# ==============================================================================

# Use one common voom object for both range models.
# This is necessary because the AIC comparison should not be affected
# by different voom precision weights.

common_range_vobj <- voomWithDreamWeights(
  dge,
  formula = ~ Site_of_distribution_range +
    (1 | Sequencing_batch) +
    (1 | Population),
  data = metadata,
  plot = FALSE,
  BPPARAM = bp_param
)

# 3-level distribution range:
# East, Centre/Center, West
aic_summaries[["distribution_range_model"]] <-
  compute_aic_for_model(
    model_name = "distribution_range_model",

    form = ~ Site_of_distribution_range +
      (1 | Sequencing_batch) +
      (1 | Population),

    vobj = common_range_vobj,

    meta = metadata,

    outdir = model_dirs[["distribution_range_model"]],

    bpparam = bp_param
  )

# 2-level distribution range:
# East versus West
aic_summaries[["distribution_range_east_west_model"]] <-
  compute_aic_for_model(
    model_name = "distribution_range_east_west_model",

    form = ~ Site_of_distribution_range_East_West +
      (1 | Sequencing_batch) +
      (1 | Population),

    vobj = common_range_vobj,

    meta = metadata,

    outdir = model_dirs[["distribution_range_east_west_model"]],

    bpparam = bp_param
  )


# ==============================================================================
# Combine all model-level AIC summaries into one file
# ==============================================================================
all_model_level_summaries <- do.call(
  rbind,
  lapply(aic_summaries, `[[`, "summary")
)

rownames(all_model_level_summaries) <- NULL

aic_summary_file <- file.path(
  degs_root,
  "AIC_summary_all_models_single.csv"
)

write.csv(
  all_model_level_summaries,
  file = aic_summary_file,
  row.names = FALSE
)

message("Combined AIC summary written to: ", aic_summary_file)


# ==============================================================================
# AIC comparison: 3-level vs 2-level distribution range models
# ==============================================================================
m3_name <- "distribution_range_model"
m2_name <- "distribution_range_east_west_model"
m3_file <- file.path(model_dirs[[m3_name]], paste0(m3_name, "_AIC_per_gene.csv"))
m2_file <- file.path(model_dirs[[m2_name]], paste0(m2_name, "_AIC_per_gene.csv"))

if (file.exists(m3_file) && file.exists(m2_file)) {
  m3_df <- read.csv(m3_file, stringsAsFactors = FALSE)
  m2_df <- read.csv(m2_file, stringsAsFactors = FALSE)
  
  if (all(c("Gene", "AIC") %in% colnames(m3_df)) && all(c("Gene", "AIC") %in% colnames(m2_df))) {
    
    # Merge by Gene
    cmp_df <- merge(
	 m3_df,
	 m2_df,
	 by = "Gene",
	 suffixes = c("_3lev", "_2lev")
	)
    
    # Compute differences: 3-level minus 2-level
    cmp_df <- cmp_df[
	 is.finite(cmp_df$AIC_3lev) &
		 is.finite(cmp_df$AIC_2lev),
	 ,
	 drop = FALSE
	]

	cmp_df$delta_AIC <-
	 cmp_df$AIC_3lev -
	 cmp_df$AIC_2lev
    
    # Save combined per-gene comparison table
    cmp_file <- file.path(degs_root, "distribution_range_3lev_vs_2lev_AIC_per_gene_single.csv")
    write.csv(cmp_df, file = cmp_file, row.names = FALSE)
    message("Per-gene 3lev vs 2lev comparison written to: ", cmp_file)
    
    # Simple summary statistics
    summary_cmp <- data.frame(
      Metric = c("N_genes", "median_delta_AIC", "mean_delta_AIC", "prop_genes_delta_AIC_lt_0"),
      Value = c(
        nrow(cmp_df),
        median(cmp_df$delta_AIC, na.rm = TRUE),
        mean(cmp_df$delta_AIC, na.rm = TRUE),
        mean(cmp_df$delta_AIC < 0, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
    summary_file <- file.path(degs_root, "distribution_range_3lev_vs_2lev_AIC_summary_single.csv")
    write.csv(summary_cmp, file = summary_file, row.names = FALSE)
    message("Summary of 3lev vs 2lev comparison written to: ", summary_file)
    
    # Density plot of delta_AIC
    pdf(file.path(degs_root, "distribution_range_3lev_vs_2lev_delta_AIC_density_single.pdf"), width = 7, height = 5)
    delta_aic <- cmp_df$delta_AIC[is.finite(cmp_df$delta_AIC)]
    if (length(delta_aic) >= 2 && length(unique(delta_aic)) >= 2) {
      dens_aic <- density(delta_aic)
      plot(dens_aic, main = "AIC difference (3-level minus 2-level)", xlab = "AIC difference", ylab = "Density", col = "steelblue")
      abline(v = 0, lty = 2, col = "darkred")
    } else {
      plot.new()
      title("Insufficient finite ΔAIC values for density plot")
    }
    dev.off()
  }
}


message("Single-predictor models completed successfully.")
message("Results written to: ", degs_root)