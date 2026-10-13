#!/usr/bin/env Rscript

# ============================================================
# Author: Iris Sammarco
# Date: 20/07/2026
# Description:
#   Assess correlations among WorldClim bioclim variables (30s resolution)
#   for Anastatica localities, identify highly correlated variable pairs,
#   and evaluate the 12-variable subset selected a priori to be used for PCA and downstream
#   DREAM models (DEG detection). Save correlation outputs for
#   reproducibility and later inspection.
# ============================================================

# ------------------------------
# 1. Load data and select bioclims
# ------------------------------

# Load the data (bioclim variables extracted for Anastatica localities)
dat <- read.csv("/Anastatica/climate/extract_coordinates/Anastatica_all_worldclim_bioclim_30s.csv")

# Select only columns whose names start with "bio" (WorldClim bioclim variables)
bioclims <- dat[, grep("^bio", names(dat))]

# Define output directory for correlation results
out_dir <- "/Anastatica/climate/bioclims_30s_correlations"

# Create output directory if it does not already exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ------------------------------
# 2. Full correlation matrix for all bioclims
# ------------------------------

# Compute Pearson correlation matrix for all bioclim variables,
# using pairwise deletion for missing values
cor_mat <- cor(bioclims, use = "pairwise.complete.obs", method = "pearson")

# Inspect the full correlation matrix (rounded to 2 decimals) in the console
print(round(cor_mat, 2))

# ------------------------------
# 3. Identify highly correlated pairs (|r| >= 0.7)
# ------------------------------

# Set absolute correlation threshold for "high" correlation
thr <- 0.7

# Find index pairs in the correlation matrix that meet the threshold
# (excluding self-correlations of 1)
high_pairs <- which(abs(cor_mat) >= thr & abs(cor_mat) < 1, arr.ind = TRUE)

# Build a data.frame listing variable pairs and their correlation coefficient
pairs_df <- data.frame(
  var1 = rownames(cor_mat)[high_pairs[, 1]],
  var2 = colnames(cor_mat)[high_pairs[, 2]],
  r    = cor_mat[high_pairs]
)

# Print highly correlated pairs to the console
print(pairs_df)

# ------------------------------
# 4. Correlations in the selected 12-variable subset
# ------------------------------

# Define the 12-variable subset selected for PCA and downstream models.
# This set captures mean temperatures, temperature range/seasonality,
# growing-season conditions, and several aspects of dry/wet quarters.
selected <- bioclims[, c("bio1","bio2","bio3","bio7","bio8","bio9",
                         "bio10","bio11","bio15","bio17","bio18","bio19")]

# Compute the correlation matrix for the selected 12-variable subset
cor_sel <- cor(selected, use = "pairwise.complete.obs")

# Print the correlation matrix for the selected subset (rounded to 2 decimals)
print(round(cor_sel, 2))

# ------------------------------
# 5. Save outputs for reproducibility
# ------------------------------

# Save the full bioclim correlation matrix as CSV
write.csv(
  cor_mat,
  file = file.path(out_dir, "Anastatica_all_worldclim_bioclim_30s_cor_mat.csv"),
  row.names = TRUE
)

# Save the highly correlated pairs table as CSV
write.csv(
  pairs_df,
  file = file.path(out_dir, "Anastatica_all_worldclim_bioclim_30s_pairs_df_thr0.7.csv"),
  row.names = FALSE
)

# Save the correlation matrix of the selected 12-variable subset as CSV
write.csv(
  cor_sel,
  file = file.path(out_dir, "Anastatica_all_worldclim_bioclim_30s_cor_sel_12vars.csv"),
  row.names = TRUE
)