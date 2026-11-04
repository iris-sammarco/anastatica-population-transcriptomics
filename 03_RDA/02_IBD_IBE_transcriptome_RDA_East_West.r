#!/usr/bin/env Rscript

# ============================================================
# AUTHOR: Iris Sammarco
# DATE: 29/10/2026
#
# REGIONAL EAST-WEST CLIMATE RDA, ADJUSTED FOR SALINITY
#
# PURPOSE
# -------
# Test whether climate is associated with population-level
# transcriptomic variation within East and West, accounting
# for soil salinity.
#
# One model is run for each of:
#   1. Full dataset
#   2. East
#   3. West
#
# Model: climate adjusted for salinity
#      Climate_PC1_z + Climate_PC2_z + Climate_PC3_z
#      + Condition(Salinity_z)
#
# The model tests the multivariate climate association
# conditional on soil salinity.
#
# IMPORTANT INTERPRETATION
# ------------------------
# Separate regional RDAs test climate associations within each
# region. They do not formally test whether climate slopes differ
# between East and West. The gene-level climate x region
# interaction model addresses that latter question.
#
# INPUTS
# ------
# Reads the population-level objects saved by:
# 01_IBD_IBE_transcriptome_RDA.R
#
# OUTPUTS
# -------
# Writes model-level results and a combined comparison table to:
# regional_East_West_RDA_analysis_salinity_adjusted/
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(readr)
  library(dplyr)
  library(tibble)
})


# ============================================================
# 1. USER SETTINGS
# ============================================================

main_output_dir <-
  "/Anastatica/RDA_transcriptome_analysis"

metadata_file <-
  "/Anastatica/metadata/population_metadata.csv"

region_column <- "Site_of_distribution_range_2-level"

regional_output_dir <-
  file.path(
    main_output_dir,
    "regional_East_West_RDA_analysis_salinity_adjusted"
  )

dir.create(
  regional_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(12345)

n_permutations <- 9999

climate_predictors <- c(
  "Climate_PC1_z",
  "Climate_PC2_z",
  "Climate_PC3_z"
)

salinity_predictor <- "Salinity_z"


# ============================================================
# 2. DEFINE THE MODEL
# ============================================================

# The model specifies:
#   predictors  = constrained predictors in the RDA
#   condition   = covariates conditioned on
#
# In this partial RDA, salinity is conditioned out before
# testing the climate predictors.

model_label <- "climate_adjusted_for_salinity"

model_spec <- list(
  predictors = climate_predictors,
  condition = salinity_predictor
)


# Require residual degrees of freedom after accounting for both
# constrained predictors and conditioning variables.
#
# This is a basic safeguard against degenerate models, not a
# guarantee of adequate statistical power.

min_extra_df <- 2


# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

safe_adj_r2 <- function(rda_object) {

  result <- vegan::RsquareAdj(rda_object)

  tibble(
    R2 = unname(result$r.squared),
    adjusted_R2 = unname(result$adj.r.squared)
  )

}


extract_anova_p <- function(x) {

  if (is.null(x)) {
    return(NA_real_)
  }

  x <- as.data.frame(x)

  if (!"Pr(>F)" %in% colnames(x)) {
    return(NA_real_)
  }

  as.numeric(x[["Pr(>F)"]][1])

}


# Run one RDA model, save its detailed outputs, and return a
# one-row summary for the combined comparison table.

run_environment_rda <- function(
    response_matrix,
    env_data,
    predictor_names,
    condition_names,
    subset_label,
    model_label,
    output_dir,
    n_permutations
) {

  n_obs <- nrow(response_matrix)
  n_pred <- length(predictor_names)
  n_condition <- length(condition_names)

  model_name <- paste(
    subset_label,
    model_label,
    sep = "_"
  )

  cat(
    "\n---",
    model_name,
    "(populations =", n_obs,
    "; constrained predictors =", n_pred,
    "; conditioning variables =", n_condition,
    ") ---\n"
  )

  # ----------------------------------------------------------
  # Check that the model has usable residual degrees of freedom
  # ----------------------------------------------------------

  minimum_n <- n_pred + n_condition + min_extra_df

  if (n_obs < minimum_n) {

    message_text <- paste0(
      "Skipped: insufficient populations. Need at least ",
      minimum_n,
      " for this model; have ",
      n_obs,
      "."
    )

    cat("  ", message_text, "\n")

    return(
      tibble(
        subset = subset_label,
        model = model_label,
        n_populations = n_obs,
        n_predictors = n_pred,
        n_conditioning_variables = n_condition,
        R2 = NA_real_,
        adjusted_R2 = NA_real_,
        P_value = NA_real_,
        status = message_text
      )
    )
  }


  # ----------------------------------------------------------
  # Build the RDA formula
  # ----------------------------------------------------------

  rhs <- paste(
    predictor_names,
    collapse = " + "
  )

  if (length(condition_names) > 0) {

    condition_term <- paste0(
      "Condition(",
      paste(condition_names, collapse = " + "),
      ")"
    )

    rhs <- paste(
      rhs,
      condition_term,
      sep = " + "
    )
  }

  model_formula <- as.formula(
    paste("response_matrix ~", rhs)
  )

  cat(
    "  Formula:",
    deparse(model_formula),
    "\n"
  )


  # ----------------------------------------------------------
  # Fit the model
  # ----------------------------------------------------------

  regional_rda <- vegan::rda(
    model_formula,
    data = env_data
  )


  # ----------------------------------------------------------
  # Overall permutation test
  #
  # This tests the climate component conditional on salinity.
  # ----------------------------------------------------------

  regional_rda_anova <- anova(
    regional_rda,
    permutations = n_permutations
  )


  # ----------------------------------------------------------
  # Marginal tests for individual constrained predictors
  #
  # Each climate PC is tested after accounting for salinity
  # and the other climate PCs.
  # ----------------------------------------------------------

  regional_rda_terms <- anova(
    regional_rda,
    by = "margin",
    permutations = n_permutations
  )


  # ----------------------------------------------------------
  # R-squared
  # ----------------------------------------------------------

  regional_rda_r2 <- safe_adj_r2(regional_rda)


  # ----------------------------------------------------------
  # Variance inflation factors
  #
  # These apply to the constrained climate predictors. Salinity
  # is a conditioning variable, not a constrained predictor.
  # ----------------------------------------------------------

  if (n_pred > 1) {

    regional_vif <- tryCatch(
      vegan::vif.cca(regional_rda),
      error = function(e) {
        message(
          "  VIF calculation failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(regional_vif)) {

      write.csv(
        data.frame(
          Predictor = names(regional_vif),
          VIF = as.numeric(regional_vif)
        ),
        file.path(
          output_dir,
          paste0(
            "RDA_",
            model_name,
            "_VIF.csv"
          )
        ),
        row.names = FALSE
      )
    }
  }


  # ----------------------------------------------------------
  # Save detailed model outputs
  # ----------------------------------------------------------

  write.csv(
    as.data.frame(regional_rda_anova),
    file.path(
      output_dir,
      paste0(
        "RDA_",
        model_name,
        "_overall_test.csv"
      )
    )
  )

  write.csv(
    as.data.frame(regional_rda_terms),
    file.path(
      output_dir,
      paste0(
        "RDA_",
        model_name,
        "_marginal_tests.csv"
      )
    )
  )

  capture.output(
    summary(regional_rda),
    file = file.path(
      output_dir,
      paste0(
        "RDA_",
        model_name,
        "_summary.txt"
      )
    )
  )


  # ----------------------------------------------------------
  # Return a summary row
  # ----------------------------------------------------------

  tibble(
    subset = subset_label,
    model = model_label,
    n_populations = n_obs,
    n_predictors = n_pred,
    n_conditioning_variables = n_condition,
    R2 = regional_rda_r2$R2,
    adjusted_R2 = regional_rda_r2$adjusted_R2,
    P_value = extract_anova_p(regional_rda_anova),
    status = "Fitted"
  )

}


# ============================================================
# 4. LOAD THE SAVED OBJECTS FROM THE MAIN ANALYSIS
# ============================================================

cat(
  "\n============================================\n",
  "Loading population-level objects\n",
  "============================================\n\n",
  sep = ""
)

main_objects_path <- file.path(
  main_output_dir,
  "RDA_complete_analysis_objects.rds"
)

if (!file.exists(main_objects_path)) {

  stop(
    paste0(
      "\nCould not find saved main-analysis objects at:\n",
      main_objects_path,
      "\nRun the main RDA_transcriptome_analysis.R script first."
    )
  )

}

main_objects <- readRDS(main_objects_path)

# Population x gene response matrix.
rda_response <- main_objects$rda_response

# Population-level standardized environmental predictors.
environment_data <- main_objects$environment_data

population_ids <- rownames(rda_response)

if (is.null(population_ids)) {
  stop("The response matrix must have population row names.")
}

if (anyDuplicated(population_ids)) {
  stop("Duplicated population IDs found in the response matrix.")
}

cat(
  "Populations loaded:",
  length(population_ids),
  "\n"
)

cat(
  "Genes in response matrix:",
  ncol(rda_response),
  "\n\n"
)


# ------------------------------------------------------------
# Check required environmental predictors
# ------------------------------------------------------------

required_predictors <- unique(
  c(climate_predictors, salinity_predictor)
)

missing_predictors <- setdiff(
  required_predictors,
  colnames(environment_data)
)

if (length(missing_predictors) > 0) {

  stop(
    paste(
      "Missing required columns in environment_data:",
      paste(missing_predictors, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# Match environment rows to response-matrix populations
# ------------------------------------------------------------

if (is.null(rownames(environment_data))) {
  stop("environment_data must have population row names.")
}

if (anyDuplicated(rownames(environment_data))) {
  stop("Duplicated population IDs found in environment_data.")
}

missing_environment <- setdiff(
  population_ids,
  rownames(environment_data)
)

if (length(missing_environment) > 0) {

  stop(
    paste(
      "Some response-matrix populations are missing from",
      "environment_data:",
      paste(missing_environment, collapse = ", ")
    )
  )
}

environment_data <- environment_data[
  population_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(environment_data),
    rownames(rda_response)
  )
)


# ------------------------------------------------------------
# Check for missing or non-finite values in model predictors
# ------------------------------------------------------------

predictor_values <- environment_data[
  ,
  required_predictors,
  drop = FALSE
]

if (!all(vapply(predictor_values, is.numeric, logical(1)))) {
  stop("All climate and salinity predictors must be numeric.")
}

if (any(!is.finite(as.matrix(predictor_values)))) {
  stop(
    "Missing or non-finite climate/salinity values detected. ",
    "Resolve these before running the RDA."
  )
}


# ============================================================
# 5. ASSIGN POPULATIONS TO EAST / WEST
# ============================================================

cat(
  "\n============================================\n",
  "Assigning populations to regions\n",
  "============================================\n\n",
  sep = ""
)

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

metadata <- as.data.frame(metadata)

if (!"Population" %in% colnames(metadata)) {
  stop("Expected a 'Population' column in the metadata file.")
}

if (!region_column %in% colnames(metadata)) {
  stop(
    paste0(
      "Expected a '",
      region_column,
      "' column in the metadata file."
    )
  )
}

metadata$Population <- as.character(metadata$Population)

region_lookup <- metadata %>%
  group_by(Population) %>%
  summarise(
    n_distinct_region = n_distinct(.data[[region_column]]),
    Region = first(.data[[region_column]]),
    .groups = "drop"
  )

if (any(region_lookup$n_distinct_region != 1)) {
  stop(
    "Some populations have missing or inconsistent region assignments."
  )
}

region_by_population <- region_lookup$Region[
  match(population_ids, region_lookup$Population)
]

if (anyNA(region_by_population)) {
  stop(
    "Could not find a region assignment for every population ",
    "in the response matrix."
  )
}

if (!all(region_by_population %in% c("East", "West"))) {
  stop(
    "The region column must contain only 'East' and 'West' ",
    "for populations used in this analysis. Check its values."
  )
}

environment_data$Region <- region_by_population

region_counts <- environment_data %>%
  count(Region, name = "n_populations")

cat("Populations per region:\n")
print(region_counts)

write_csv(
  region_counts,
  file.path(
    regional_output_dir,
    "region_population_counts.csv"
  )
)


# ============================================================
# 6. RUN THE MODEL FOR FULL DATASET, EAST AND WEST
# ============================================================

cat(
  "\n============================================\n",
  "Running regional climate RDAs adjusted for salinity\n",
  "============================================\n",
  sep = ""
)

subsets <- list(
  Full_dataset = population_ids,
  East = population_ids[region_by_population == "East"],
  West = population_ids[region_by_population == "West"]
)

all_results <- list()

for (subset_label in names(subsets)) {

  subset_ids <- subsets[[subset_label]]

  subset_response <- rda_response[
    subset_ids,
    ,
    drop = FALSE
  ]

  subset_env <- environment_data[
    subset_ids,
    ,
    drop = FALSE
  ]

  result_row <- run_environment_rda(
    response_matrix = subset_response,
    env_data = subset_env,
    predictor_names = model_spec$predictors,
    condition_names = model_spec$condition,
    subset_label = subset_label,
    model_label = model_label,
    output_dir = regional_output_dir,
    n_permutations = n_permutations
  )

  all_results[[
    paste(subset_label, model_label, sep = "_")
  ]] <- result_row

}


# ============================================================
# 7. SAVE THE COMBINED COMPARISON TABLE
# ============================================================

comparison_table <- bind_rows(all_results)

write_csv(
  comparison_table,
  file.path(
    regional_output_dir,
    "regional_climate_salinity_RDA_comparison.csv"
  )
)


# ============================================================
# 8. FINAL SUMMARY
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("REGIONAL CLIMATE RDA ADJUSTED FOR SALINITY - SUMMARY\n")
cat("============================================================\n\n")

print(comparison_table)

cat("\nModel:\n")
cat(
  "climate_adjusted_for_salinity = climate PCs + Condition(Salinity_z)\n"
)

cat("\nCombined comparison table saved to:\n")
cat(
  file.path(
    regional_output_dir,
    "regional_climate_salinity_RDA_comparison.csv"
  ),
  "\n"
)

cat("\nDetailed model outputs saved in:\n")
cat(regional_output_dir, "\n")

cat("\nInterpretation reminder:\n")
cat(
  "- The partial RDA tests climate conditional on salinity.\n",
  "- Regional RDAs do not formally test East-West slope differences.\n",
  "- Non-significant regional results may reflect limited power.\n",
  "- Use the gene-level climate x region interaction model to test\n",
  "  whether climate-response slopes differ between regions.\n",
  sep = ""
)

cat("============================================================\n")
