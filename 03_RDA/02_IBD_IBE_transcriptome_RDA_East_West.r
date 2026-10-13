#!/usr/bin/env Rscript

# ============================================================
# Anastatica hierochuntica
#
# SUPPLEMENTARY / EXPLORATORY ANALYSIS:
# REGIONAL (EAST vs WEST) ENVIRONMENT-ONLY RDA
#
# ============================================================
#
# AUTHOR: Iris Sammarco
#
# DATE: 25/09/2026
#
# ============================================================
#
# PURPOSE
# -------
#
# The primary analysis (01_IBD_IBE_transcriptome_RDA.R) asked,
# across the full set of 29 populations, whether transcriptomic
# differentiation is better explained by environment (climate +
# salinity) or by geographic/spatial structure (dbMEMs). It found
# a significant environment effect (driven by Climate_PC1 and
# Climate_PC2) and no independent spatial effect once environment
# was accounted for.
#
# This script is a SUPPLEMENTARY, EXPLORATORY robustness check. It asks:
#
#   "If populations are split into their two broad distribution
#    regions (East vs West), does the environment effect detected
#    in the full dataset still hold within each region on its
#    own, or does it look different?"
#
# It does NOT replace or override the main analysis. It is a
# sensitivity check on the main conclusion.
#
# ============================================================
#
# IMPORTANT LIMITATIONS - READ BEFORE INTERPRETING RESULTS
# ----------------------------------------------------------
#
# 1. SAMPLE SIZE. Regional sample sizes are small (East ~16,
#    West ~13 populations, see below), which gives limited
#    residual degrees of freedom for a multi-predictor RDA.
#    This analysis therefore has much lower statistical power
#    than the full-dataset analysis. A non-significant regional
#    result does NOT mean the environmental effect is absent or
#    inconsistent in that region - it may simply reflect reduced
#    power. Results here should be read as exploratory /
#    hypothesis-generating, not as a confirmatory replication.
#
# 2. THE "CENTRE" GROUP. The 3-level distribution classification
#    (East / Centre / West) has only 2 Centre populations, which
#    cannot be analysed on their own. This script therefore uses
#    the pre-existing 2-level classification in the metadata
#    ("Site_of_distribution_range_2-level"), in which the Centre
#    populations are already grouped together with East.
#
# 3. SPATIAL (dbMEM) FRAMEWORK NOT REPEATED HERE. Forward
#    selection of spatial eigenvectors was already only
#    marginally well-powered at n = 29 populations; at n = 13-16
#    it would not be informative and could easily produce a
#    misleading null result. This script therefore re-tests ONLY
#    the environment effect within each region, not the spatial
#    effect. The main analysis remains the authority on the
#    spatial question.
#
# 4. TWO PREDICTOR SETS ARE TESTED, FOR TRANSPARENCY:
#      - "full":    Climate_PC1_z + Climate_PC2_z + Climate_PC3_z
#                    + Salinity_z (matches the main model exactly)
#      - "reduced": Climate_PC1_z + Climate_PC2_z only, since
#                    these were the only two predictors significant
#                    in the marginal tests of the full-dataset
#                    model. Dropping the two non-significant terms
#                    preserves more residual degrees of freedom,
#                    which matters a lot at these small regional
#                    sample sizes.
#    Both are also re-run on the FULL dataset (all populations)
#    using the identical function, so that the full-dataset,
#    East-only and West-only results are directly comparable -
#    same code, same settings, only the population subset differs.
#
# ============================================================
#
# WHAT THIS SCRIPT DOES, STEP BY STEP
# ------------------------------------
#
#   1. Loads the population-level transcriptomic response matrix
#      and standardized environmental predictors that were saved
#      by the main analysis script (no need to redo expression
#      processing / batch correction / standardization).
#   2. Assigns each population to East or West using the 2-level
#      regional classification in the metadata file.
#   3. For each of {Full dataset, East, West} x {full predictor
#      set, reduced predictor set}, runs an environment-only RDA:
#      overall permutation test, per-term marginal tests, adjusted
#      R2, and VIFs (where more than one predictor is used).
#   4. Skips (with a clear message, not an error) any combination
#      that does not have enough populations to support the number
#      of predictors requested.
#   5. Collects everything into one combined comparison table, so
#      it is easy to see at a glance whether the direction and
#      significance of the environment effect is consistent
#      across the full dataset, East, and West.
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

# Directory where the MAIN analysis script wrote its outputs.
# This script reads saved objects from there and does not repeat
# expression processing, batch correction, or standardization.
main_output_dir <-
  "/Anastatica/IBD_IBE_transcriptome_analysis"

# Sample-level (or population-level) metadata file containing the
# regional classification column. Same file used by the main script.
metadata_file <-
  "/Anastatica/metadata/population_metadata.csv"

# Column in the metadata that gives the 2-level region assignment
# (Centre already merged into East). Kept as a variable in case the
# column gets renamed later.
region_column <- "Site_of_distribution_range_2-level"

# Where this supplementary analysis's own outputs will be written.
# Kept as a subfolder of the main output directory so everything
# related to one project stays together.
regional_output_dir <-
  file.path(
    main_output_dir,
    "regional_East_West_RDA_analysis"
  )

dir.create(
  regional_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Reproducibility - same seed as the main script
# ------------------------------------------------------------

set.seed(12345)


# ------------------------------------------------------------
# Permutations for the RDA significance tests.
# NOTE: with small regional sample sizes the number of *distinct*
# possible permutations is still far larger than 9999 (13! is an
# enormous number), so this does not need to be reduced.
# ------------------------------------------------------------

n_permutations <- 9999


# ------------------------------------------------------------
# Predictor sets to test (see limitation #4 above for rationale)
# ------------------------------------------------------------

predictor_sets <-
  list(

    full =
      c(
        "Climate_PC1_z",
        "Climate_PC2_z",
        "Climate_PC3_z",
        "Salinity_z"
      ),

    reduced =
      c(
        "Climate_PC1_z",
        "Climate_PC2_z"
      )

  )


# A model needs at least (n_predictors + 2) populations to have any
# usable residual degrees of freedom. Combinations below this are
# skipped rather than forced to run.
min_extra_df <- 2


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

# Pulls the adjusted R2 out of a vegan RDA object in a tidy form.
safe_adj_r2 <- function(rda_object) {

  result <- vegan::RsquareAdj(rda_object)

  tibble(
    R2          = unname(result$r.squared),
    adjusted_R2 = unname(result$adj.r.squared)
  )

}


# Pulls the overall model P-value out of an anova.cca() result.
extract_anova_p <- function(x) {

  if (is.null(x)) return(NA_real_)

  x <- as.data.frame(x)

  if (!"Pr(>F)" %in% colnames(x)) return(NA_real_)

  as.numeric(x$`Pr(>F)`[1])

}


# Runs one environment-only RDA on a given subset of populations
# and a given set of predictors, saves the detailed output files,
# and returns a one-row summary for the combined comparison table.
#
# Arguments:
#   response_matrix  - population x gene matrix (subset of rda_response)
#   env_data         - data.frame of standardized predictors (matching
#                       rows/order of response_matrix)
#   predictor_names  - character vector of column names in env_data
#   subset_label     - short label used in file names / summary table
#                       (e.g. "Full_dataset", "East", "West")
#   predictor_label  - "full" or "reduced" (which predictor set)
#   output_dir       - where to write the detailed result files
#   n_permutations   - number of permutations for the significance tests
run_environment_rda <- function(
    response_matrix,
    env_data,
    predictor_names,
    subset_label,
    predictor_label,
    output_dir,
    n_permutations
) {

  model_label <-
    paste(subset_label, predictor_label, sep = "_")

  n_obs <- nrow(response_matrix)
  n_pred <- length(predictor_names)

  cat(
    "\n---",
    model_label,
    "( n populations =", n_obs, ", n predictors =", n_pred, ") ---\n"
  )

  # ----------------------------------------------------------
  # Skip cleanly if there is not enough sample size to support
  # this predictor set, rather than letting rda() fail or return
  # a degenerate model.
  # ----------------------------------------------------------

  if (n_obs < (n_pred + min_extra_df)) {

    cat(
      "  Skipped: not enough populations for",
      n_pred,
      "predictors (need at least",
      n_pred + min_extra_df,
      ", have",
      n_obs,
      ").\n"
    )

    return(
      tibble(
        subset          = subset_label,
        predictor_set   = predictor_label,
        n_populations   = n_obs,
        n_predictors    = n_pred,
        R2              = NA_real_,
        adjusted_R2     = NA_real_,
        P_value         = NA_real_,
        status          = "Skipped: insufficient sample size"
      )
    )

  }

  # ----------------------------------------------------------
  # Build and fit the RDA model
  # ----------------------------------------------------------

  model_formula <-
    as.formula(
      paste(
        "response_matrix ~",
        paste(predictor_names, collapse = " + ")
      )
    )

  regional_rda <-
    vegan::rda(
      model_formula,
      data = env_data
    )

  regional_rda_anova <-
    anova(
      regional_rda,
      permutations = n_permutations
    )

  regional_rda_terms <-
    anova(
      regional_rda,
      by = "margin",
      permutations = n_permutations
    )

  regional_rda_r2 <-
    safe_adj_r2(regional_rda)

  # VIFs are only meaningful with more than one predictor.
  if (n_pred > 1) {

    regional_vif <-
      vegan::vif.cca(regional_rda)

    write.csv(
      data.frame(
        Predictor = names(regional_vif),
        VIF       = as.numeric(regional_vif)
      ),
      file.path(
        output_dir,
        paste0("RDA_", model_label, "_VIF.csv")
      ),
      row.names = FALSE
    )

  }

  # ----------------------------------------------------------
  # Save detailed outputs for this model
  # ----------------------------------------------------------

  write.csv(
    as.data.frame(regional_rda_anova),
    file.path(
      output_dir,
      paste0("RDA_", model_label, "_overall_test.csv")
    )
  )

  write.csv(
    as.data.frame(regional_rda_terms),
    file.path(
      output_dir,
      paste0("RDA_", model_label, "_marginal_tests.csv")
    )
  )

  capture.output(
    summary(regional_rda),
    file =
      file.path(
        output_dir,
        paste0("RDA_", model_label, "_summary.txt")
      )
  )

  # ----------------------------------------------------------
  # One-row summary for the combined comparison table
  # ----------------------------------------------------------

  tibble(
    subset          = subset_label,
    predictor_set   = predictor_label,
    n_populations   = n_obs,
    n_predictors    = n_pred,
    R2              = regional_rda_r2$R2,
    adjusted_R2     = regional_rda_r2$adjusted_R2,
    P_value         = extract_anova_p(regional_rda_anova),
    status          = "Fitted"
  )

}


# ============================================================
# 3. LOAD SAVED OBJECTS FROM THE MAIN ANALYSIS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Loading population-level objects from the main analysis\n"
)

cat(
  "============================================\n\n"
)

main_objects_path <-
  file.path(
    main_output_dir,
    "IBD_IBE_complete_analysis_objects.rds"
  )

if (!file.exists(main_objects_path)) {

  stop(
    paste0(
      "\nCould not find saved main-analysis objects at:\n",
      main_objects_path,
      "\nRun the main IBD_IBE_transcriptome_analysis.R script first."
    )
  )

}

main_objects <- readRDS(main_objects_path)

# Population-level transcriptomic response matrix (populations x
# genes, gene-standardized) - identical object used in the main RDAs.
rda_response <- main_objects$rda_response

# Population-level standardized environmental predictors, one row
# per population, rownames = Population.
environment_data <- main_objects$environment_data

population_ids <- rownames(rda_response)

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


# ============================================================
# 4. ASSIGN POPULATIONS TO REGIONS (EAST / WEST)
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Assigning populations to regions\n"
)

cat(
  "============================================\n\n"
)

metadata <-
  read_csv(
    metadata_file,
    show_col_types = FALSE
  )

metadata <- as.data.frame(metadata)

if (!"Population" %in% colnames(metadata)) {

  stop(
    "\nExpected a 'Population' column in the metadata file."
  )

}

if (!region_column %in% colnames(metadata)) {

  stop(
    paste0(
      "\nExpected a '",
      region_column,
      "' column in the metadata file."
    )
  )

}

metadata$Population <- as.character(metadata$Population)

# Collapse to one region value per population, and check that the
# region assignment is actually consistent within each population
# (it should be, since region is a population-level property).
region_lookup <-
  metadata %>%

  group_by(Population) %>%

  summarise(

    n_distinct_region =
      n_distinct(.data[[region_column]]),

    Region =
      first(.data[[region_column]]),

    .groups = "drop"

  )

if (any(region_lookup$n_distinct_region > 1)) {

  stop(
    "\nSome populations have more than one region value in the metadata. ",
    "Check the region column for inconsistencies before proceeding."
  )

}

# Match region assignment onto the populations actually present in
# the transcriptomic response matrix, in the same order.
region_by_population <-
  region_lookup$Region[
    match(population_ids, region_lookup$Population)
  ]

if (anyNA(region_by_population)) {

  stop(
    "\nCould not find a region assignment for one or more populations ",
    "in the transcriptomic response matrix."
  )

}

# Attach region to the environment data (needed so that subsetting
# by region also subsets environment_data consistently).
environment_data$Region <- region_by_population

region_counts <-
  environment_data %>%
  count(Region, name = "n_populations")

cat(
  "Populations per region:\n"
)

print(region_counts)

write_csv(
  region_counts,
  file.path(
    regional_output_dir,
    "region_population_counts.csv"
  )
)


# ============================================================
# 5. RUN ENVIRONMENT-ONLY RDA FOR EACH SUBSET x PREDICTOR SET
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Running environment-only RDA models\n"
)

cat(
  "============================================\n"
)

# Named list of population subsets to loop over. "Full_dataset" is
# included so that it is generated with EXACTLY the same function
# and settings as the regional subsets, making the comparison
# table below a fair like-for-like comparison rather than mixing
# numbers computed in different ways.
subsets <-
  list(

    Full_dataset =
      population_ids,

    East =
      population_ids[region_by_population == "East"],

    West =
      population_ids[region_by_population == "West"]

  )

# Collect one summary row per (subset x predictor set) combination.
all_results <- list()

for (subset_label in names(subsets)) {

  subset_ids <- subsets[[subset_label]]

  subset_response <-
    rda_response[subset_ids, , drop = FALSE]

  subset_env <-
    environment_data[subset_ids, , drop = FALSE]

  for (predictor_label in names(predictor_sets)) {

    predictor_names <-
      predictor_sets[[predictor_label]]

    result_row <-
      run_environment_rda(
        response_matrix = subset_response,
        env_data         = subset_env,
        predictor_names  = predictor_names,
        subset_label     = subset_label,
        predictor_label  = predictor_label,
        output_dir       = regional_output_dir,
        n_permutations   = n_permutations
      )

    all_results[[
      paste(subset_label, predictor_label, sep = "_")
    ]] <- result_row

  }

}

comparison_table <-
  bind_rows(all_results)

write_csv(
  comparison_table,
  file.path(
    regional_output_dir,
    "regional_environment_RDA_comparison.csv"
  )
)


# ============================================================
# 6. FINAL CONSOLE SUMMARY
# ============================================================

cat(
  "\n\n"
)

cat(
  "============================================================\n"
)

cat(
  "REGIONAL (EAST vs WEST) ENVIRONMENT RDA - SUMMARY\n"
)

cat(
  "============================================================\n\n"
)

print(comparison_table)

cat(
  "\nFull results table saved to:\n",
  file.path(
    regional_output_dir,
    "regional_environment_RDA_comparison.csv"
  ),
  "\n\n"
)

cat(
  "REMINDER: regional sample sizes are small (see counts above).\n"
)

cat(
  "Treat any regional non-significant result as a power limitation,\n"
)

cat(
  "not as evidence that the environmental effect is absent in that\n"
)

cat(
  "region. This script is a supplementary robustness check; the\n"
)

cat(
  "main analysis remains the primary result.\n\n"
)

cat(
  "============================================================\n"
)