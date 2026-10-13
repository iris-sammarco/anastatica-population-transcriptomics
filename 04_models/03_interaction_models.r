#!/usr/bin/env Rscript

# ==============================================================================
# Date: 01/10/2026
# Author: Iris Sammarco
# ==============================================================================
#
# Purpose:
# --------
# Differential expression analysis of Anastatica greenhouse-grown populations
# using dream (variancePartition) with mixed-effects models.
#
# MAIN SCIENTIFIC QUESTION:
# -------------------------
# Do climate-expression relationships differ between the eastern and western
# parts of the species' distribution range?
#
# IMPORTANT:
# ----------
#
# East and West are compared using the SAME climate PC axes calculated from
# the global dataset. Climate PCs are therefore not recalculated separately
# within East and West.
#
# Models:
# -------
# 1) climate_range_interaction_east_west_model
#
#    ~ Climate_PC1 * East_West +
#      Climate_PC2 * East_West +
#      Climate_PC3 * East_West +
#      (1 | Sequencing_batch) +
#      (1 | Population)
#
# 2) ECe_climate_range_interaction_east_west_model
#
#    ~ ECe_mean +
#      Climate_PC1 * East_West +
#      Climate_PC2 * East_West +
#      Climate_PC3 * East_West +
#      (1 | Sequencing_batch) +
#      (1 | Population)
#
# Interpretation with West as reference:
#
#    Climate_PC1
#       = West climate slope
#
#    Climate_PC1:East
#       = East - West difference in climate slope
#
#    Climate_PC1 + Climate_PC1:East
#       = East climate slope
#
# The same applies to PC2 and PC3.
#
# The script therefore explicitly tests:
#
#    1. West climate effects
#    2. East climate effects
#    3. East-West differences in climate effects
#    4. omnibus climate effects in West
#    5. omnibus climate effects in East
#    6. omnibus East-West climate interaction
#
# It also identifies:
#    - West-only climate-associated genes
#    - East-only climate-associated genes
#    - shared climate-associated genes
#    - genes not climate-associated in either region
#
# For shared genes, effect directions and magnitudes are compared between
# East and West.
#
# Population is retained as a random intercept because individuals are nested
# within populations.
#
# Sequencing_batch is retained as a random effect.
#
# ==============================================================================


# ==============================================================================
# Library paths and loading
# ==============================================================================

message("R version: ", R.version.string)
message("variancePartition: ", as.character(packageVersion("variancePartition")))
message("lme4: ", as.character(packageVersion("lme4")))
message("lmerTest: ", as.character(packageVersion("lmerTest")))
message("Matrix: ", as.character(packageVersion("Matrix")))

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

model_names <- c(
  "climate_range_interaction_east_west_model",
  "ECe_climate_range_interaction_east_west_model"
)

model_dirs <- setNames(
  file.path(degs_root, model_names),
  model_names
)

for (d in model_dirs) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# Load data
# ==============================================================================

if (!file.exists(gse_rds)) {
  stop("Gene-level RDS not found at: ", gse_rds)
}

gse <- readRDS(gse_rds)

cts <- as.matrix(
  assay(gse, "counts")
)

metadata <- as.data.frame(
  colData(gse),
  stringsAsFactors = FALSE
)

match_col <- "quant_dir_name"

if (!match_col %in% colnames(metadata)) {
  stop(
    "metadata must contain a '",
    match_col,
    "' column."
  )
}

if (!all(colnames(cts) %in% metadata[[match_col]])) {
  stop(
    "The following count columns are missing from metadata: ",
    paste(
      setdiff(
        colnames(cts),
        metadata[[match_col]]
      ),
      collapse = ", "
    )
  )
}

rownames(metadata) <- metadata[[match_col]]

metadata <- metadata[
  colnames(cts),
  ,
  drop = FALSE
]


# ==============================================================================
# Metadata preparation
# ==============================================================================

prepare_metadata <- function(meta) {

  meta <- as.data.frame(
    meta,
    stringsAsFactors = FALSE
  )

  required_vars <- c(
    "Sample_name",
    "Population",
    "Sequencing_batch",
    "Site_of_distribution_range",
    "Site_of_distribution_range_East_West",
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3",
    "ECe_mean"
  )

  missing_vars <- setdiff(
    required_vars,
    colnames(meta)
  )

  if (length(missing_vars) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  }

  meta$Population <- factor(
    meta$Population
  )

  meta$Sequencing_batch <- factor(
    meta$Sequencing_batch
  )

  meta$Site_of_distribution_range <- factor(
    meta$Site_of_distribution_range
  )

  meta$Site_of_distribution_range_East_West <- factor(
    meta$Site_of_distribution_range_East_West
  )

  cont_vars <- c(
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3",
    "ECe_mean"
  )

  for (v in cont_vars) {
    meta[[v]] <- as.numeric(
      meta[[v]]
    )
  }

  check_vars <- c(
    "Population",
    "Sequencing_batch",
    "Site_of_distribution_range_East_West",
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3",
    "ECe_mean"
  )

  if (
    any(
      is.na(
        meta[, check_vars]
      )
    )
  ) {
    stop(
      "Metadata contains NA values in variables used for modelling."
    )
  }

  meta
}

metadata <- prepare_metadata(
  metadata
)


# ==============================================================================
# Check East-West population structure
# ==============================================================================

message("\n============================================")
message("East-West sample/population structure")
message("============================================")

print(
  table(
    metadata$Site_of_distribution_range_East_West
  )
)

pop_by_region <- table(
  metadata$Population,
  metadata$Site_of_distribution_range_East_West
)

print(pop_by_region)

write.csv(
  as.data.frame(pop_by_region),
  file = file.path(
    degs_root,
    "East_West_population_structure.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# Set reference level
# ==============================================================================

if (
  !all(
    c("East", "West") %in%
      levels(
        metadata$Site_of_distribution_range_East_West
      )
  )
) {
  stop(
    "East and West must both be present. Found: ",
    paste(
      levels(
        metadata$Site_of_distribution_range_East_West
      ),
      collapse = ", "
    )
  )
}

metadata$Site_of_distribution_range_East_West <-
  droplevels(
    metadata$Site_of_distribution_range_East_West
  )

metadata$Site_of_distribution_range_East_West <-
  relevel(
    metadata$Site_of_distribution_range_East_West,
    ref = "West"
  )

message(
  "\nReference level for East-West factor: West"
)


# ==============================================================================
# Gene filtering and normalization
# ==============================================================================

design_filter <- model.matrix(
  ~ Sequencing_batch,
  metadata
)

keep <- filterByExpr(
  cts,
  design = design_filter
)

cts_filt <- cts[
  keep,
  ,
  drop = FALSE
]

message(
  "Genes before filtering: ",
  nrow(cts)
)

message(
  "Genes retained after filtering: ",
  nrow(cts_filt)
)

dge <- DGEList(
  counts = cts_filt
)

dge <- calcNormFactors(
  dge
)


# ==============================================================================
# Helper functions
# ==============================================================================

save_voom_summary <- function(
    vobj,
    outfile,
    model_name
) {

  capture.output({

    cat(
      "Voom object summary for model:",
      model_name,
      "\n\n"
    )

    str(
      vobj,
      max.level = 1
    )

    cat(
      "\nDimensions of expression matrix:",
      dim(vobj$E),
      "\n"
    )

    cat(
      "Dimensions of weights matrix:",
      dim(vobj$weights),
      "\n"
    )

  }, file = outfile)
}

add_gene_id_column <- function(df) {

  if (!"FeatureID" %in% colnames(df)) {
    df$FeatureID <- rownames(df)
  }

  df <- df[
    ,
    c("FeatureID", setdiff(colnames(df), "FeatureID")),
    drop = FALSE
  ]

  rownames(df) <- NULL
  df
}

write_de_sets <- function(
    df,
    prefix,
    outdir,
    fdr_thresh,
    lfc_thresh,
    apply_lfc = TRUE
) {

  if (is.null(df) || nrow(df) == 0) {
    return(invisible(NULL))
  }

  df <- add_gene_id_column(
    df
  )

  write.csv(
    df,
    file = file.path(
      outdir,
      paste0(
        prefix,
        "_full.csv"
      )
    ),
    row.names = FALSE
  )

  df_fdr <- subset(
    df,
    adj.P.Val <= fdr_thresh
  )

  write.csv(
    df_fdr,
    file = file.path(
      outdir,
      paste0(
        prefix,
        "_FDR05.csv"
      )
    ),
    row.names = FALSE
  )

  if (
    apply_lfc &&
    "logFC" %in% colnames(df)
  ) {

    df_fdr_lfc <- subset(
      df,
      adj.P.Val <= fdr_thresh &
        !is.na(logFC) &
        abs(logFC) >= lfc_thresh
    )

    write.csv(
      df_fdr_lfc,
      file = file.path(
        outdir,
        paste0(
          prefix,
          "_FDR05_log2FC_1.5.csv"
        )
      ),
      row.names = FALSE
    )
  }

  invisible(NULL)
}


summarize_de_results <- function(
    df,
    model_name
) {

  if (
    is.null(df) ||
    nrow(df) == 0
  ) {
    return(NULL)
  }

  required_cols <- c(
    "Test_type",
    "Test",
    "adj.P.Val"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Result table is missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  test_groups <- unique(
    df[
      ,
      c(
        "Test_type",
        "Test"
      ),
      drop = FALSE
    ]
  )

  summary_list <- lapply(
    seq_len(nrow(test_groups)),
    function(i) {

      test_type_i <-
        test_groups$Test_type[i]

      test_label_i <-
        test_groups$Test[i]

      test_df <- df[
        df$Test_type == test_type_i &
          df$Test == test_label_i,
        ,
        drop = FALSE
      ]

      tested <- !is.na(
        test_df$adj.P.Val
      )

      sig_lfc <- if (
        "logFC" %in%
        colnames(test_df)
      ) {

        sum(
          tested &
            test_df$adj.P.Val <=
              fdr_thresh &
            !is.na(test_df$logFC) &
            abs(test_df$logFC) >=
              lfc_thresh
        )

      } else {
        NA_integer_
      }

      data.frame(
        Model = model_name,
        Test_type =
          as.character(test_type_i),
        Test_label =
          as.character(test_label_i),
        N_genes_tested =
          sum(tested),
        N_sig_FDR05 =
          sum(
            tested &
              test_df$adj.P.Val <=
                fdr_thresh
          ),
        N_sig_FDR05_log2FC_1.5 =
          sig_lfc,
        stringsAsFactors = FALSE
      )
    }
  )

  summary_df <- do.call(
    rbind,
    summary_list
  )

  rownames(summary_df) <- NULL

  summary_df
}


# ==============================================================================
# Find interaction coefficient name
# ==============================================================================

get_interaction_name <- function(
    fit,
    pc
) {

  candidate_names <- c(
    paste0(
      pc,
      ":Site_of_distribution_range_East_WestEast"
    ),
    paste0(
      "Site_of_distribution_range_East_WestEast:",
      pc
    )
  )

  found <- intersect(
    candidate_names,
    colnames(
      coef(fit)
    )
  )

  if (length(found) != 1) {
    stop(
      "Could not uniquely identify interaction coefficient for ",
      pc,
      ". Candidates: ",
      paste(
        candidate_names,
        collapse = ", "
      ),
      ". Available coefficients: ",
      paste(
        colnames(coef(fit)),
        collapse = ", "
      )
    )
  }

  found
}


# ==============================================================================
# Extract coefficient results
# ==============================================================================

extract_coefficient_results <- function(
    fit,
    coefs_to_test
) {

  coefs_to_test <- intersect(
    coefs_to_test,
    colnames(
      coef(fit)
    )
  )

  if (length(coefs_to_test) == 0) {
    return(NULL)
  }

  tt_all <- lapply(
    coefs_to_test,
    function(cf) {

      tt <- variancePartition::topTable(
        fit,
        coef = cf,
        number = Inf,
        sort.by = "P"
      )

      tt$Test <- cf
      tt$Test_type <- "coefficient"

      tt$FeatureID <- rownames(tt)
      tt
    }
  )

  names(tt_all) <- coefs_to_test

  res_df <- do.call(
    rbind,
    tt_all
  )

  res_df <- res_df[
    order(
      res_df$Test,
      res_df$P.Value
    ),
    ,
    drop = FALSE
  ]

  res_df
}


# ==============================================================================
# Build East-specific climate contrasts
#
# With West as the reference:
#
#   East_PC1 = Climate_PC1 +
#              Climate_PC1:East
#
#   East_PC2 = Climate_PC2 +
#              East:Climate_PC2
#
#   East_PC3 = Climate_PC3 +
#              East:Climate_PC3
#
# The contrasts are constructed directly from the coefficient names rather
# than parsed as strings. This avoids problems caused by ":" in interaction
# coefficient names and avoids partial string replacement.
# ==============================================================================

build_east_climate_contrast_matrix <- function(fit) {

  coef_names <- colnames(
    coef(fit)
  )

  climate_pcs <- c(
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3"
  )

  interaction_names <- setNames(
    vapply(
      climate_pcs,
      function(pc) {
        get_interaction_name(
          fit,
          pc
        )
      },
      character(1)
    ),
    climate_pcs
  )

  contrast_matrix <- matrix(
    0,
    nrow = length(coef_names),
    ncol = length(climate_pcs),
    dimnames = list(
      coef_names,
      paste0(
        "East_",
        climate_pcs
      )
    )
  )

  for (pc in climate_pcs) {

    contrast_matrix[
      pc,
      paste0("East_", pc)
    ] <- 1

    contrast_matrix[
      interaction_names[[pc]],
      paste0("East_", pc)
    ] <- 1
  }

  contrast_matrix
}


# ==============================================================================
# Fit East-specific climate slopes
#
# IMPORTANT:
# For dream mixed models, contrasts are supplied to dream() itself.
#
# East slope:
#
#   Climate_PC1 + Climate_PC1:East
#   Climate_PC2 + East:Climate_PC2
#   Climate_PC3 + East:Climate_PC3
#
# A manually constructed L matrix is used because these are sums of
# regression coefficients. makeContrastsDream() requires its contrast
# weights to sum to 1.
# ==============================================================================

fit_east_climate_contrasts <- function(
    fit_raw,
    vobj,
    form,
    meta,
    BPPARAM
) {

  east_contrast_matrix <-
    build_east_climate_contrast_matrix(
      fit_raw
    )

  message(
    "\nEast climate contrast matrix:"
  )

  print(
    east_contrast_matrix[
      rowSums(
        abs(east_contrast_matrix)
      ) > 0,
      ,
      drop = FALSE
    ]
  )

  east_fit <- variancePartition::dream(
    vobj,
    form,
    meta,
    L = east_contrast_matrix,
    BPPARAM = BPPARAM,
    computeResiduals = TRUE,
    ddf = "Satterthwaite",
    hideErrorsInBackend = FALSE
  )

  east_fit <- variancePartition::eBayes(
    east_fit
  )

  east_fit
}

# ==============================================================================
# Extract omnibus F-test
# ==============================================================================

extract_omnibus_results <- function(
    fit,
    coef_vector,
    test_name
) {

  coef_vector <- intersect(
    coef_vector,
    colnames(
      coef(fit)
    )
  )

  if (length(coef_vector) == 0) {
    return(NULL)
  }

  tt <- variancePartition::topTable(
    fit,
    coef = coef_vector,
    number = Inf,
    sort.by = "P"
  )

  tt$Test <- test_name
  tt$Test_type <- "omnibus_Ftest"
  tt$FeatureID <- rownames(tt)
  tt
}


# ==============================================================================
# Run model
# ==============================================================================

run_model <- function(
    model_name,
    form,
    coefficient_terms = NULL,
    omnibus_coef = NULL,
    omnibus_name = NULL,
    dge_obj = dge,
    meta = metadata,
    do_east_west_climate_outputs = TRUE
) {

  message(
    "\n============================================"
  )

  message(
    "Running model: ",
    model_name
  )

  message(
    "============================================"
  )


  # --------------------------------------------------------------------------
  # Output directory
  # --------------------------------------------------------------------------

  outdir <- model_dirs[[model_name]]


  # --------------------------------------------------------------------------
  # Save metadata
  # --------------------------------------------------------------------------

  write.csv(
    meta,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_sample_metadata_used.csv"
      )
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------------------------
  # Design matrix
  # --------------------------------------------------------------------------

  fixed_form <- reformulas::nobars(
    form
  )

  design_mat <- model.matrix(
    fixed_form,
    meta
  )

  write.csv(
    design_mat,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_design_matrix.csv"
      )
    ),
    row.names = TRUE
  )

  qr_design <- qr(
    design_mat
  )

  if (
    qr_design$rank <
      ncol(design_mat)
  ) {

    dropped_columns <-
      colnames(design_mat)[
        qr_design$pivot[
          seq.int(
            qr_design$rank + 1,
            ncol(design_mat)
          )
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
      paste(
        dropped_columns,
        collapse = "\n"
      )
    )
  }


  # --------------------------------------------------------------------------
  # voom weights
  # --------------------------------------------------------------------------

  vobj <- voomWithDreamWeights(
    dge_obj,
    formula = form,
    data = meta,
    plot = FALSE,
    BPPARAM = bp_param
  )

  weights_mat <- vobj$weights

  if (
    !is.null(weights_mat) &&
    is.null(
      colnames(weights_mat)
    )
  ) {

    colnames(weights_mat) <-
      colnames(
        dge_obj$counts
      )
  }

  save_voom_summary(
    vobj = vobj,
    outfile = file.path(
      outdir,
      paste0(
        model_name,
        "_voom_weights_summary.txt"
      )
    ),
    model_name = model_name
  )

  if (
    is.null(vobj$E) ||
    nrow(vobj$E) == 0
  ) {

    stop(
      "voomWithDreamWeights returned an empty expression matrix for model: ",
      model_name
    )
  }


  # --------------------------------------------------------------------------
  # Fit dream model
  #
  # IMPORTANT:
  # Save the RAW dream fit before eBayes because East slopes require
  # linear contrasts of coefficients.
  # --------------------------------------------------------------------------

  old_warn <- getOption(
    "warn"
  )

  options(warn = 1)

  fit_raw <- suppressWarnings(
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


  saveRDS(
    fit_raw,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_fit_raw.rds"
      )
    )
  )


  # --------------------------------------------------------------------------
  # eBayes fit
  # --------------------------------------------------------------------------

  fit_base <- variancePartition::eBayes(
    fit_raw
  )

  saveRDS(
    fit_base,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_fit_base_ebayes.rds"
      )
    )
  )

  saveRDS(
    vobj,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_vobj.rds"
      )
    )
  )


  # --------------------------------------------------------------------------
  # Model coefficients
  # --------------------------------------------------------------------------

  coef_df <- as.data.frame(
    coef(fit_base)
  )

  coef_df <- add_gene_id_column(
    coef_df
  )

  write.csv(
    coef_df,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_model_coefficients.csv"
      )
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------------------------
  # Standard coefficient tests
  # --------------------------------------------------------------------------

  available_tests <- colnames(
    coef(fit_base)
  )

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

  de_full <- NULL

  if (
    length(fixed_coefs) > 0
  ) {

    de_full <- extract_coefficient_results(
      fit_base,
      fixed_coefs
    )

    write_de_sets(
      df = de_full,
      prefix = paste0(
        model_name,
        "_DE_results"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh
    )
  }


  # --------------------------------------------------------------------------
  # Existing omnibus test
  # --------------------------------------------------------------------------

  omnibus_full <- NULL

  if (
    !is.null(omnibus_coef) &&
    length(omnibus_coef) > 0
  ) {

    available_base_tests <-
      colnames(
        coef(fit_base)
      )

    valid_omnibus <- intersect(
      omnibus_coef,
      available_base_tests
    )

    if (
      length(valid_omnibus) == 0
    ) {

      stop(
        "None of the omnibus coefficients were found in the BASE fitted model: ",
        model_name
      )
    }

    omnibus_full <-
      extract_omnibus_results(
        fit = fit_base,
        coef_vector = valid_omnibus,
        test_name = omnibus_name
      )

    write_de_sets(
      df = omnibus_full,
      prefix = paste0(
        model_name,
        "_omnibus_Ftest"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh,
      apply_lfc = FALSE
    )
  }


  # ==========================================================================
  # EAST-WEST CLIMATE-SPECIFIC OUTPUTS
  # ==========================================================================

  ew_outputs <- NULL

  if (
    do_east_west_climate_outputs
  ) {

    message(
      "\n--------------------------------------------"
    )

    message(
      "Extracting East-West climate-specific effects"
    )

    message(
      "--------------------------------------------"
    )


    # ------------------------------------------------------------------------
    # Identify interaction coefficients
    # ------------------------------------------------------------------------

    climate_pcs <- c(
      "Climate_PC1",
      "Climate_PC2",
      "Climate_PC3"
    )

    interaction_names <- setNames(
      vapply(
        climate_pcs,
        function(pc) {
          get_interaction_name(
            fit_base,
            pc
          )
        },
        character(1)
      ),
      climate_pcs
    )

    message(
      "Interaction coefficients:"
    )

    print(
      interaction_names
    )


    # ------------------------------------------------------------------------
    # WEST slopes
    #
    # Because West is the reference:
    #
    # Climate_PC1 = West slope
    # ------------------------------------------------------------------------

    west_results_list <- lapply(
      climate_pcs,
      function(pc) {

        tt <- variancePartition::topTable(
          fit_base,
          coef = pc,
          number = Inf,
          sort.by = "P"
        )

        tt$Test <- paste0(
          "West_",
          pc
        )

        tt$Test_type <- "West_climate_slope"
		tt$FeatureID <- rownames(tt)
        tt
      }
    )

    west_results <- do.call(
      rbind,
      west_results_list
    )

    write_de_sets(
      west_results,
      prefix = paste0(
        model_name,
        "_West_climate_slopes"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh
    )


    # ------------------------------------------------------------------------
    # EAST slopes
    #
    # East slope =
    #
    #   Climate_PCx + Climate_PCx:East
    #
    # All three East-specific slopes are fitted simultaneously.
    # ------------------------------------------------------------------------

    east_climate_fit <-
      fit_east_climate_contrasts(
        fit_raw = fit_raw,
        vobj = vobj,
        form = form,
        meta = meta,
        BPPARAM = bp_param
      )

    east_results_list <- lapply(
      paste0(
        "East_",
        climate_pcs
      ),
      function(east_coef) {

        tt <- variancePartition::topTable(
          east_climate_fit,
          coef = east_coef,
          number = Inf,
          sort.by = "P"
        )

        tt$Test <- east_coef
        tt$Test_type <- "East_climate_slope"
        tt$FeatureID <- rownames(tt)
        tt
		

      }
    )

    east_results <- do.call(
      rbind,
      east_results_list
    )

    write_de_sets(
      east_results,
      prefix = paste0(
        model_name,
        "_East_climate_slopes"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh
    )

message(
  "\nChecking East-specific climate contrasts..."
)

east_test_check <- data.frame(
  Climate_PC = climate_pcs,
  N_logFC = sapply(
    paste0(
      "East_",
      climate_pcs
    ),
    function(x) {
      sum(
        !is.na(
          east_climate_fit$coefficients[, x]
        )
      )
    }
  ),
  N_P = sapply(
    paste0(
      "East_",
      climate_pcs
    ),
    function(x) {
      sum(
        !is.na(
          east_climate_fit$p.value[, x]
        )
      )
    }
  )
)

print(
  east_test_check
)

if (
  any(
    east_test_check$N_logFC == 0 |
      east_test_check$N_P == 0
  )
) {
  stop(
    "One or more East climate contrasts failed. ",
    "Check the contrast matrix and model coefficients."
  )
}
    # ------------------------------------------------------------------------
    # EAST-WEST INTERACTION effects
    # ------------------------------------------------------------------------

    interaction_results_list <- lapply(
      climate_pcs,
      function(pc) {

        interaction_name <-
          interaction_names[[pc]]

        tt <- variancePartition::topTable(
          fit_base,
          coef = interaction_name,
          number = Inf,
          sort.by = "P"
        )

        tt$Test <- paste0(
          "East_minus_West_",
          pc
        )

        tt$Test_type <- "East_West_interaction"
        tt$FeatureID <- rownames(tt)
        tt
      }
    )

    interaction_results <- do.call(
      rbind,
      interaction_results_list
    )

    write_de_sets(
      interaction_results,
      prefix = paste0(
        model_name,
        "_East_West_interactions"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh
    )


    # ------------------------------------------------------------------------
    # WEST omnibus climate test
    #
    # Tests:
    #
    # H0:
    # beta_PC1 = beta_PC2 = beta_PC3 = 0
    #
    # This asks whether West has any climate-expression association.
    # ------------------------------------------------------------------------

    west_omnibus <- extract_omnibus_results(
      fit = fit_base,
      coef_vector = climate_pcs,
      test_name = "West_climate_omnibus"
    )

    write_de_sets(
      west_omnibus,
      prefix = paste0(
        model_name,
        "_West_climate_omnibus"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh,
      apply_lfc = FALSE
    )

    # ------------------------------------------------------------------------
    # EAST omnibus climate test
    #
    # The East-specific contrast fit has already been calculated above.
    #
    # Test:
    #
    #   H0:
    #   East_PC1 = East_PC2 = East_PC3 = 0
    #
    # This is the omnibus climate test for East.
    # ------------------------------------------------------------------------

    east_omnibus <- variancePartition::topTable(
      east_climate_fit,
      coef = paste0(
        "East_",
        climate_pcs
      ),
      number = Inf,
      sort.by = "P"
    )

    east_omnibus$Test <-
      "East_climate_omnibus"

    east_omnibus$Test_type <-
      "omnibus_Ftest"
	 
	east_omnibus$FeatureID <- rownames(east_omnibus)

    write_de_sets(
      east_omnibus,
      prefix = paste0(
        model_name,
        "_East_climate_omnibus"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh,
      apply_lfc = FALSE
    )
 
    # ------------------------------------------------------------------------
    # EAST-WEST INTERACTION OMNIBUS TEST
    #
    # Tests:
    #
    # H0:
    # all three East-West climate slope differences = 0
    # ------------------------------------------------------------------------

    interaction_omnibus <- extract_omnibus_results(
      fit = fit_base,
      coef_vector =
        unname(
          interaction_names
        ),
      test_name =
        "East_West_climate_interaction_omnibus"
    )

    write_de_sets(
      interaction_omnibus,
      prefix = paste0(
        model_name,
        "_East_West_climate_interaction_omnibus"
      ),
      outdir = outdir,
      fdr_thresh = fdr_thresh,
      lfc_thresh = lfc_thresh,
      apply_lfc = FALSE
    )


    # ------------------------------------------------------------------------
    # Save separate East/West result tables
    # ------------------------------------------------------------------------

    write.csv(
      add_gene_id_column(
        west_results
      ),
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_West_climate_results.csv"
        )
      ),
      row.names = FALSE
    )

    write.csv(
      add_gene_id_column(
        east_results
      ),
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_climate_results.csv"
        )
      ),
      row.names = FALSE
    )

    write.csv(
      add_gene_id_column(
        interaction_results
      ),
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_interaction_results.csv"
        )
      ),
      row.names = FALSE
    )


    # =========================================================================
    # GENE-LEVEL EAST / WEST CLIMATE EFFECT TABLE
    # =========================================================================

    message(
      "Building gene-level East/West climate effect table..."
    )

    all_gene_ids <- rownames(
      coef(fit_base)
    )

    effects_df <- data.frame(
      FeatureID = all_gene_ids,
      stringsAsFactors = FALSE
    )

    for (pc in climate_pcs) {

      interaction_name <-
        interaction_names[[pc]]

      # ----------------------------------------------------------------------
      # West slope
      # ----------------------------------------------------------------------

      west_tt <- variancePartition::topTable(
        fit_base,
        coef = pc,
        number = Inf,
        sort.by = "none"
      )

      west_tt <- west_tt[
        match(
          all_gene_ids,
          rownames(west_tt)
        ),
        ,
        drop = FALSE
      ]

      effects_df[
        ,
        paste0(
          "West_",
          pc,
          "_logFC"
        )
      ] <- west_tt$logFC

      effects_df[
        ,
        paste0(
          "West_",
          pc,
          "_P"
        )
      ] <- west_tt$P.Value

      effects_df[
        ,
        paste0(
          "West_",
          pc,
          "_FDR"
        )
      ] <- west_tt$adj.P.Val


      # ----------------------------------------------------------------------
      # East-West difference
      # ----------------------------------------------------------------------

      interaction_tt <-
        variancePartition::topTable(
          fit_base,
          coef = interaction_name,
          number = Inf,
          sort.by = "none"
        )

      interaction_tt <-
        interaction_tt[
          match(
            all_gene_ids,
            rownames(interaction_tt)
          ),
          ,
          drop = FALSE
        ]

      effects_df[
        ,
        paste0(
          "East_minus_West_",
          pc,
          "_logFC"
        )
      ] <- interaction_tt$logFC

      effects_df[
        ,
        paste0(
          "East_minus_West_",
          pc,
          "_P"
        )
      ] <- interaction_tt$P.Value

      effects_df[
        ,
        paste0(
          "East_minus_West_",
          pc,
          "_FDR"
        )
      ] <- interaction_tt$adj.P.Val


      # ----------------------------------------------------------------------
      # East slope
      # ----------------------------------------------------------------------

      east_tt <-
        east_results[
          east_results$Test ==
            paste0(
              "East_",
              pc
            ),
          ,
          drop = FALSE
        ]

      east_tt <- east_tt[
        match(all_gene_ids, east_tt$FeatureID),
        ,
        drop = FALSE
      ]

      if (
        all(
          is.na(
            east_tt$logFC
          )
        )
      ) {
        stop(
            "East slope extraction failed for ",
          pc,
          ". All East logFC values are NA."
        )
      }

      effects_df[
        ,
        paste0(
          "East_",
          pc,
          "_logFC"
        )
      ] <- east_tt$logFC

      effects_df[
        ,
        paste0(
          "East_",
          pc,
          "_P"
        )
      ] <- east_tt$P.Value

      effects_df[
        ,
        paste0(
          "East_",
          pc,
          "_FDR"
        )
      ] <- east_tt$adj.P.Val
    }


    write.csv(
      effects_df,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_climate_effects_by_gene.csv"
        )
      ),
      row.names = FALSE
    )


    # =========================================================================
    # CLASSIFY GENES AS EAST / WEST / SHARED / NEITHER
    #
    # Classification is based on the omnibus climate test in each region.
    #
    # This is preferable to taking the minimum FDR across PC1-PC3 because the
    # omnibus test directly asks whether ANY of the three climate axes is
    # associated with expression in that region.
    # =========================================================================

    west_omni_df <- west_omnibus

    east_omni_df <- east_omnibus

    west_omni_df <- west_omni_df[
      match(
        all_gene_ids,
        rownames(west_omni_df)
      ),
      ,
      drop = FALSE
    ]

    east_omni_df <- east_omni_df[
      match(
        all_gene_ids,
        rownames(east_omni_df)
      ),
      ,
      drop = FALSE
    ]

    classification_df <- data.frame(
      FeatureID = all_gene_ids,
      stringsAsFactors = FALSE
    )

    classification_df$West_climate_FDR <-
      west_omni_df$adj.P.Val

    classification_df$East_climate_FDR <-
      east_omni_df$adj.P.Val

    classification_df$West_climate_associated <-
      classification_df$West_climate_FDR <=
      fdr_thresh

    classification_df$East_climate_associated <-
      classification_df$East_climate_FDR <=
      fdr_thresh
	
	classification_df$West_climate_associated[
      is.na(classification_df$West_climate_associated)] <- FALSE
    classification_df$East_climate_associated[
      is.na(classification_df$East_climate_associated)] <- FALSE

    classification_df$Climate_class <- "neither"

    classification_df$Climate_class[
      classification_df$West_climate_associated &
        !classification_df$East_climate_associated
    ] <- "West_only"

    classification_df$Climate_class[
      !classification_df$West_climate_associated &
        classification_df$East_climate_associated
    ] <- "East_only"

    classification_df$Climate_class[
      classification_df$West_climate_associated &
        classification_df$East_climate_associated
    ] <- "shared"

    classification_df$Model_failed <-
      is.na(classification_df$West_climate_FDR) |
      is.na(classification_df$East_climate_FDR)

    classification_df$Climate_class[
      classification_df$Model_failed
    ] <- "not_tested"

    write.csv(
      classification_df,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_climate_gene_classification.csv"
        )
      ),
      row.names = FALSE
    )


    # ------------------------------------------------------------------------
    # Class counts
    # ------------------------------------------------------------------------

    class_counts <- as.data.frame(
      table(
        classification_df$Climate_class
      )
    )

    colnames(class_counts) <- c(
      "Climate_class",
      "N_genes"
    )

    write.csv(
      class_counts,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_climate_gene_classification_counts.csv"
        )
      ),
      row.names = FALSE
    )


    # =========================================================================
    # PER-PC EAST / WEST OVERLAP
    # =========================================================================

    per_pc_overlap <- list()

    for (pc in climate_pcs) {

      west_fdr_col <- paste0(
        "West_",
        pc,
        "_FDR"
      )

      east_fdr_col <- paste0(
        "East_",
        pc,
        "_FDR"
      )

      west_sig <- effects_df[
        !is.na(
          effects_df[[west_fdr_col]]
        ) &
          effects_df[[west_fdr_col]] <=
          fdr_thresh,
        "FeatureID"
      ]

      east_sig <- effects_df[
        !is.na(
          effects_df[[east_fdr_col]]
        ) &
          effects_df[[east_fdr_col]] <=
          fdr_thresh,
        "FeatureID"
      ]

      overlap_df <- data.frame(
        PC = pc,
        West_only = length(
          setdiff(
            west_sig,
            east_sig
          )
        ),
        East_only = length(
          setdiff(
            east_sig,
            west_sig
          )
        ),
        Shared = length(
          intersect(
            west_sig,
            east_sig
          )
        ),
        West_total = length(
          west_sig
        ),
        East_total = length(
          east_sig
        ),
        stringsAsFactors = FALSE
      )

      per_pc_overlap[[pc]] <-
        overlap_df
    }

    per_pc_overlap_df <- do.call(
      rbind,
      per_pc_overlap
    )

    rownames(
      per_pc_overlap_df
    ) <- NULL

    write.csv(
      per_pc_overlap_df,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_climate_overlap_by_PC.csv"
        )
      ),
      row.names = FALSE
    )


    # =========================================================================
    # SHARED-GENE DIRECTION CONCORDANCE
    # =========================================================================

    shared_genes <- classification_df$FeatureID[
      classification_df$Climate_class ==
        "shared"
    ]

    shared_effects <- effects_df[
      effects_df$FeatureID %in%
        shared_genes,
      ,
      drop = FALSE
    ]

    if (nrow(shared_effects) > 0) {

      for (pc in climate_pcs) {

        east_col <- paste0(
          "East_",
          pc,
          "_logFC"
        )

        west_col <- paste0(
          "West_",
          pc,
          "_logFC"
        )

        concordance_col <- paste0(
          pc,
          "_direction_concordance"
        )

        shared_effects[
          ,
          concordance_col
        ] <- ifelse(
          is.na(
            shared_effects[[east_col]]
          ) |
            is.na(
              shared_effects[[west_col]]
            ),
          NA,
          ifelse(
            sign(
              shared_effects[[east_col]]
            ) ==
              sign(
                shared_effects[[west_col]]
              ),
            "concordant",
            "opposite"
          )
        )
      }

      write.csv(
        shared_effects,
        file = file.path(
          outdir,
          paste0(
            model_name,
            "_shared_climate_genes_East_West_effects.csv"
          )
        ),
        row.names = FALSE
      )
    }


    # =========================================================================
    # SUMMARY OF EAST / WEST CLIMATE RESULTS
    # =========================================================================

    ew_summary_rows <- list()

    for (pc in climate_pcs) {

      west_tt <- west_results[
        west_results$Test ==
          paste0(
            "West_",
            pc
          ),
        ,
        drop = FALSE
      ]

      east_tt <- east_results[
        east_results$Test ==
          paste0(
            "East_",
            pc
          ),
        ,
        drop = FALSE
      ]

      interaction_tt <-
        interaction_results[
          interaction_results$Test ==
            paste0(
              "East_minus_West_",
              pc
            ),
          ,
          drop = FALSE
        ]

      ew_summary_rows[[pc]] <- data.frame(
        Model = model_name,
        Climate_PC = pc,

        West_N_genes_tested =
          sum(
            !is.na(
              west_tt$adj.P.Val
            )
          ),

        West_N_sig_FDR05 =
          sum(
            !is.na(
              west_tt$adj.P.Val
            ) &
              west_tt$adj.P.Val <=
              fdr_thresh
          ),

        West_N_sig_FDR05_log2FC_1.5 =
          sum(
            !is.na(
              west_tt$adj.P.Val
            ) &
              west_tt$adj.P.Val <=
              fdr_thresh &
              !is.na(
                west_tt$logFC
              ) &
              abs(
                west_tt$logFC
              ) >=
              lfc_thresh
          ),

        East_N_genes_tested =
          sum(
            !is.na(
              east_tt$adj.P.Val
            )
          ),

        East_N_sig_FDR05 =
          sum(
            !is.na(
              east_tt$adj.P.Val
            ) &
              east_tt$adj.P.Val <=
              fdr_thresh
          ),

        East_N_sig_FDR05_log2FC_1.5 =
          sum(
            !is.na(
              east_tt$adj.P.Val
            ) &
              east_tt$adj.P.Val <=
              fdr_thresh &
              !is.na(
                east_tt$logFC
              ) &
              abs(
                east_tt$logFC
              ) >=
              lfc_thresh
          ),

        Interaction_N_genes_tested =
          sum(
            !is.na(
              interaction_tt$adj.P.Val
            )
          ),

        Interaction_N_sig_FDR05 =
          sum(
            !is.na(
              interaction_tt$adj.P.Val
            ) &
              interaction_tt$adj.P.Val <=
              fdr_thresh
          ),

        Interaction_N_sig_FDR05_log2FC_1.5 =
          sum(
            !is.na(
              interaction_tt$adj.P.Val
            ) &
              interaction_tt$adj.P.Val <=
              fdr_thresh &
              !is.na(
                interaction_tt$logFC
              ) &
              abs(
                interaction_tt$logFC
              ) >=
              lfc_thresh
          ),

        stringsAsFactors = FALSE
      )
    }

    ew_summary <- do.call(
      rbind,
      ew_summary_rows
    )

    rownames(ew_summary) <- NULL

    write.csv(
      ew_summary,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_climate_effect_summary.csv"
        )
      ),
      row.names = FALSE
    )


    # ------------------------------------------------------------------------
    # Save omnibus summaries
    # ------------------------------------------------------------------------

    omnibus_summary <- data.frame(

      Model = model_name,

      West_climate_omnibus_N_sig_FDR05 =
        sum(
          west_omnibus$adj.P.Val <=
            fdr_thresh,
          na.rm = TRUE
        ),

      East_climate_omnibus_N_sig_FDR05 =
        sum(
          east_omnibus$adj.P.Val <=
            fdr_thresh,
          na.rm = TRUE
        ),

      East_West_interaction_omnibus_N_sig_FDR05 =
        sum(
          interaction_omnibus$adj.P.Val <=
            fdr_thresh,
          na.rm = TRUE
        ),

      stringsAsFactors = FALSE
    )

    write.csv(
      omnibus_summary,
      file = file.path(
        outdir,
        paste0(
          model_name,
          "_East_West_omnibus_summary.csv"
        )
      ),
      row.names = FALSE
    )


    ew_outputs <- list(
      west = west_results,
      east = east_results,
      interaction = interaction_results,
      west_omnibus = west_omnibus,
      east_omnibus = east_omnibus,
      interaction_omnibus =
        interaction_omnibus,
      effects_by_gene =
        effects_df,
      classification =
        classification_df,
      overlap_by_PC =
        per_pc_overlap_df
    )
  }


  # ==========================================================================
  # Variance partitioning
  # ==========================================================================

  expr_log2 <- vobj$E

  if (
    model_name ==
      "climate_range_interaction_east_west_model"
  ) {

    form_vp <-
      ~ Climate_PC1 +
      Climate_PC2 +
      Climate_PC3 +
      (1 |
         Site_of_distribution_range_East_West) +
      (1 |
         Sequencing_batch) +
      (1 |
         Population)

  } else if (
    model_name ==
      "ECe_climate_range_interaction_east_west_model"
  ) {

    form_vp <-
      ~ ECe_mean +
      Climate_PC1 +
      Climate_PC2 +
      Climate_PC3 +
      (1 |
         Site_of_distribution_range_East_West) +
      (1 |
         Sequencing_batch) +
      (1 |
         Population)

  } else {

    stop(
      "No variance-partition formula defined for model: ",
      model_name
    )
  }


  varPart <- fitExtractVarPartModel(
    expr_log2,
    form_vp,
    meta,
    BPPARAM = bp_param
  )

  write.csv(
    as.data.frame(varPart),
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_varPart_results.csv"
      )
    ),
    row.names = TRUE
  )

  vp_summary <- apply(
    varPart,
    2,
    function(x) {

      c(
        mean =
          mean(
            x,
            na.rm = TRUE
          ),
        median =
          median(
            x,
            na.rm = TRUE
          ),
        sd =
          sd(
            x,
            na.rm = TRUE
          ),
        min =
          min(
            x,
            na.rm = TRUE
          ),
        max =
          max(
            x,
            na.rm = TRUE
          )
      )
    }
  )

  capture.output({

    cat(
      "Variance fraction summary for model:",
      model_name,
      "\n\n"
    )

    print(
      round(
        vp_summary,
        4
      )
    )

  },
  file = file.path(
    outdir,
    paste0(
      model_name,
      "_varPart_summary.txt"
    )
  ))

  pdf(
    file.path(
      outdir,
      paste0(
        model_name,
        "_varPart_plot.pdf"
      )
    ),
    width = 7,
    height = 5
  )

  print(
    plotVarPart(
      sortCols(varPart)
    )
  )

  dev.off()


  # ==========================================================================
  # Diagnostics
  # ==========================================================================

  res_mat <- residuals(
    fit_base,
    y = vobj
  )

  res_var <- apply(
    res_mat,
    1,
    function(x) {
      var(
        x,
        na.rm = TRUE
      )
    }
  )

  median_gene <-
    names(res_var)[
      order(res_var)[
        ceiling(
          length(res_var) / 2
        )
      ]
    ]

  expr_mat <- vobj$E

  diag_df <- data.frame(
    Observed =
      expr_mat[
        median_gene,
        ],
    Residuals =
      res_mat[
        median_gene,
        ]
  )

  pdf(
    file.path(
      outdir,
      paste0(
        model_name,
        "_diagnostics_residuals_vs_observed.pdf"
      )
    ),
    width = 7,
    height = 5
  )

  print(
    ggplot(
      diag_df,
      aes(
        x = Observed,
        y = Residuals
      )
    ) +
      geom_point(
        alpha = 0.6
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      theme_bw() +
      labs(
        title = paste(
          model_name,
          "- residuals vs observed"
        ),
        x = "Observed log-expression",
        y = "Residuals"
      )
  )

  dev.off()


  diag_df_sqrt <- data.frame(
    Observed =
      expr_mat[
        median_gene,
        ],
    SqrtAbsRes =
      sqrt(
        abs(
          res_mat[
            median_gene,
            ]
        )
      )
  )

  pdf(
    file.path(
      outdir,
      paste0(
        model_name,
        "_diagnostics_scale_location.pdf"
      )
    ),
    width = 7,
    height = 5
  )

  print(
    ggplot(
      diag_df_sqrt,
      aes(
        x = Observed,
        y = SqrtAbsRes
      )
    ) +
      geom_point(
        alpha = 0.6
      ) +
      theme_bw() +
      labs(
        title = paste(
          model_name,
          "- scale-location"
        ),
        x = "Observed log-expression",
        y = "sqrt(|Residuals|)"
      )
  )

  dev.off()


  qq_df <- data.frame(
    Sample =
      sort(
        diag_df$Residuals
      ),
    Theoretical =
      sort(
        qnorm(
          ppoints(
            length(
              diag_df$Residuals
            )
          )
        )
      )
  )

  pdf(
    file.path(
      outdir,
      paste0(
        model_name,
        "_diagnostics_residuals_QQ.pdf"
      )
    ),
    width = 7,
    height = 5
  )

  print(
    ggplot(
      qq_df,
      aes(
        x = Theoretical,
        y = Sample
      )
    ) +
      geom_point() +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      theme_bw() +
      labs(
        title = paste(
          model_name,
          "- residual QQ plot"
        ),
        x = "Theoretical quantiles",
        y = "Sample residuals"
      )
  )

  dev.off()


  res_summary <- data.frame(
    Gene = names(res_var),
    Residual_variance =
      res_var
  )

  res_summary$Residual_sd <-
    sqrt(
      res_summary$Residual_variance
    )

  write.csv(
    res_summary,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_diagnostics_residual_summary.csv"
      )
    ),
    row.names = FALSE
  )


  # ==========================================================================
  # Random-effect diagnostics
  # ==========================================================================

  if (
    !is.null(
      fit_base$random
    )
  ) {

    re_example <-
      fit_base$random[[1]]

    if (!is.null(re_example)) {

      re_list <-
        lme4::ranef(
          re_example
        )

      if (
        "Population" %in%
        names(re_list)
      ) {

        pop_re <-
		  re_list[["Population"]]

        pop_re$Level <-
          rownames(pop_re)

        colnames(pop_re) <-
          paste0(
            "RE_",
            colnames(pop_re)
          )

        write.csv(
          pop_re,
          file = file.path(
            outdir,
            paste0(
              model_name,
              "_diagnostics_random_effects_Population_example.csv"
            )
          ),
          row.names = FALSE
        )
      }

      if (
        "Sequencing_batch" %in%
        names(re_list)
      ) {

        batch_re <-
		  re_list[["Sequencing_batch"]]

        batch_re$Level <-
          rownames(batch_re)

        colnames(batch_re) <-
          paste0(
            "RE_",
            colnames(batch_re)
          )

        write.csv(
          batch_re,
          file = file.path(
            outdir,
            paste0(
              model_name,
              "_diagnostics_random_effects_Sequencing_batch_example.csv"
            )
          ),
          row.names = FALSE
        )
      }
    }
  }


  # ==========================================================================
  # Overall DE summary
  # ==========================================================================

  summary_tables <- list(
    coefficient = de_full,
    omnibus = omnibus_full
  )

  summary_tables <-
    summary_tables[
      !vapply(
        summary_tables,
        is.null,
        logical(1)
      )
    ]

  de_summary <- NULL

  if (
    length(summary_tables) > 0
  ) {

    de_summary <- do.call(
      rbind,
      lapply(
        summary_tables,
        summarize_de_results,
        model_name =
          model_name
      )
    )

    rownames(
      de_summary
    ) <- NULL
  }

  write.csv(
    de_summary,
    file = file.path(
      outdir,
      paste0(
        model_name,
        "_DE_summary.csv"
      )
    ),
    row.names = FALSE
  )


  # ==========================================================================
  # Return
  # ==========================================================================

  invisible(
    list(
      fit_raw = fit_raw,
      fit = fit_base,
      fit_base = fit_base,
      vobj = vobj,
      de_full = de_full,
      omnibus_full = omnibus_full,
      de_summary = de_summary,
      varPart = varPart,
      east_west = ew_outputs
    )
  )
}


# ==============================================================================
# Define interaction model
# ==============================================================================

climate_range_interaction_east_west_formula_base <-
  ~ Climate_PC1 *
    Site_of_distribution_range_East_West +
    Climate_PC2 *
    Site_of_distribution_range_East_West +
    Climate_PC3 *
    Site_of_distribution_range_East_West +
    (1 | Sequencing_batch) +
    (1 | Population)


# ==============================================================================
# Expected coefficient names
# ==============================================================================

climate_range_interaction_east_west_coef_names <- c(
  "Climate_PC1",
  "Climate_PC2",
  "Climate_PC3",
  "Site_of_distribution_range_East_WestEast",
  "Climate_PC1:Site_of_distribution_range_East_WestEast",
  "Site_of_distribution_range_East_WestEast:Climate_PC2",
  "Site_of_distribution_range_East_WestEast:Climate_PC3"
)


# ==============================================================================
# Inspect design matrix
# ==============================================================================

inspect_design <- function(
    formula,
    meta,
    model_name
) {

  cat(
    "\n============================================\n"
  )

  cat(
    "Design matrix:",
    model_name,
    "\n"
  )

  cat(
    "============================================\n"
  )

  fixed_form <- reformulas::nobars(
    formula
  )

  design_mat <- model.matrix(
    fixed_form,
    meta
  )

  cat(
    "\nFormula:\n"
  )

  print(
    fixed_form
  )

  cat(
    "\nDesign matrix dimensions:\n"
  )

  print(
    dim(design_mat)
  )

  cat(
    "\nDesign matrix columns:\n"
  )

  print(
    colnames(design_mat)
  )

  cat(
    "\nRank:\n"
  )

  qr_design <- qr(
    design_mat
  )

  print(
    qr_design$rank
  )

  cat(
    "\nNumber of columns:\n"
  )

  print(
    ncol(design_mat)
  )

  if (
    qr_design$rank <
      ncol(design_mat)
  ) {

    dropped_columns <-
      colnames(design_mat)[
        qr_design$pivot[
          seq.int(
            qr_design$rank + 1,
            ncol(design_mat)
          )
        ]
      ]

    cat(
      "\nWARNING: rank-deficient design.\n"
    )

    cat(
      "Potentially redundant columns:\n"
    )

    print(
      dropped_columns
    )

  } else {

    cat(
      "\nDesign is full rank.\n"
    )
  }

  invisible(
    design_mat
  )
}


# ==============================================================================
# Inspect climate-only interaction model
# ==============================================================================

inspect_design(
  climate_range_interaction_east_west_formula_base,
  metadata,
  "climate_range_interaction_east_west_model"
)


# ==============================================================================
# Inspect ECe + climate interaction model
# ==============================================================================

ECe_formula <- update(
  climate_range_interaction_east_west_formula_base,
  ~ ECe_mean + .
)

inspect_design(
  ECe_formula,
  metadata,
  "ECe_climate_range_interaction_east_west_model"
)


# ==============================================================================
# RUN MODEL 1
# ==============================================================================

climate_range_interaction_east_west_res <-
  run_model(
    model_name =
      "climate_range_interaction_east_west_model",

    form =
      climate_range_interaction_east_west_formula_base,

    coefficient_terms =
      climate_range_interaction_east_west_coef_names,

    omnibus_coef =
      c(
        "Climate_PC1:Site_of_distribution_range_East_WestEast",
        "Site_of_distribution_range_East_WestEast:Climate_PC2",
        "Site_of_distribution_range_East_WestEast:Climate_PC3"
      ),

    omnibus_name =
      "Climate_by_East_West_interaction_omnibus",

    dge_obj =
      dge,

    meta =
      metadata,

    do_east_west_climate_outputs =
      TRUE
  )


# ==============================================================================
# RUN MODEL 2
# ECe-adjusted East-West climate interaction model
# ==============================================================================

ECe_climate_range_interaction_east_west_res <-
  run_model(
    model_name =
      "ECe_climate_range_interaction_east_west_model",

    form =
      ECe_formula,

    coefficient_terms =
      c(
        "ECe_mean",
        climate_range_interaction_east_west_coef_names
      ),

    omnibus_coef =
      c(
        "Climate_PC1:Site_of_distribution_range_East_WestEast",
        "Site_of_distribution_range_East_WestEast:Climate_PC2",
        "Site_of_distribution_range_East_WestEast:Climate_PC3"
      ),

    omnibus_name =
      "ECe_climate_by_East_West_interaction_omnibus",

    dge_obj =
      dge,

    meta =
      metadata,

    do_east_west_climate_outputs =
      TRUE
  )


# ==============================================================================
# FINAL MESSAGE
# ==============================================================================

message(
  "\n============================================================"
)

message(
  "East-West climate interaction models completed successfully."
)

message(
  "Global Climate PC axes were retained."
)