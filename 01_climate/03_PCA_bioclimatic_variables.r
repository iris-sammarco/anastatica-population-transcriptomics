#!/usr/bin/env Rscript

# ============================================================
# Author: Iris Sammarco
# Date: 06/07/2026
# PCA on selected WorldClim bioclim variables (12 uncorrelated Bioclims)
# Output: scores, loadings, eigenvalues, variance explained
# ============================================================

library(readr)
library(dplyr)
library(tibble)
library(stats)
library(tidyr)

# ----------------------------
# User settings
# ----------------------------
input_file <- "/Anastatica/climate/extract_coordinates/Anastatica_all_worldclim_bioclim_30s.csv"
output_dir  <- "/Anastatica/climate/pca_output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Read data
# ----------------------------
dat <- read_csv(input_file, show_col_types = FALSE)

# ----------------------------
# Select variables
# ----------------------------
vars <- c(
  "bio1",  "bio2",  "bio3",  "bio7",
  "bio8",  "bio9",  "bio10", "bio11",
  "bio15", "bio17", "bio18", "bio19"
)

missing_vars <- setdiff(vars, names(dat))
if (length(missing_vars) > 0) {
  stop(paste("These variables are missing from the input data:", paste(missing_vars, collapse = ", ")))
}

# Keep only rows with complete climatic data for PCA
pca_dat <- dat %>%
  select(all_of(vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  drop_na()

# Preserve row identity for complete cases
complete_idx <- which(complete.cases(dat[, vars]))
id_dat <- dat[complete_idx, , drop = FALSE]

# ----------------------------
# Run PCA
# ----------------------------
pca_res <- prcomp(pca_dat, center = TRUE, scale. = TRUE)

# ----------------------------
# Extract scores
# ----------------------------
scores <- as.data.frame(pca_res$x) %>%
  rownames_to_column("pca_row_id")

scores <- bind_cols(id_dat, scores)

# ----------------------------
# Extract loadings
# ----------------------------
loadings <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column("variable")

# ----------------------------
# Variance explained
# ----------------------------
eigvals <- pca_res$sdev^2
var_exp <- eigvals / sum(eigvals)
cum_var_exp <- cumsum(var_exp)

variance_table <- tibble(
  PC = paste0("PC", seq_along(eigvals)),
  eigenvalue = eigvals,
  proportion_variance = var_exp,
  cumulative_variance = cum_var_exp
)

# ----------------------------
# Save outputs
# ----------------------------
write_csv(scores, file.path(output_dir, "bioclim_pca_scores.csv"))
write_csv(loadings, file.path(output_dir, "bioclim_pca_loadings.csv"))
write_csv(variance_table, file.path(output_dir, "bioclim_pca_variance_explained.csv"))

saveRDS(pca_res, file.path(output_dir, "bioclim_pca_object.rds"))

# ----------------------------
# Optional: merge back scores to original data
# ----------------------------
# If you want the full table with PCs for modelling:
# keep only rows used in PCA, then merge as needed
model_dat <- scores %>%
  select(locality_code, country_code, site_of_distribution_range,
         latitude, longitude, starts_with("PC"))

write_csv(model_dat, file.path(output_dir, "bioclim_pca_model_table.csv"))

# ----------------------------
# Print a compact summary
# ----------------------------
cat("\nPCA completed successfully.\n")
cat("Rows used in PCA:", nrow(pca_dat), "\n")
cat("Variables used:", paste(vars, collapse = ", "), "\n")
print(variance_table)