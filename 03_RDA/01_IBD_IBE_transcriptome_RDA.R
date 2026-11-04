#!/usr/bin/env Rscript

# ============================================================
# AUTHOR: Iris Sammarco
# DATE: 29/10/2026
#
# Anastatica hierochuntica
# Population-level transcriptomic association with climate,
# salinity, and spatial structure.
#
# Main models:
#   1. Climate-only RDA: climate PCs 1 + 2 + 3
#   2. Salinity-only RDA: ECe_mean
#   3. Spatial-only RDA: selected overall positive dbMEMs
#   4. Partial RDA: climate | salinity + spatial
#   5. Partial RDA: salinity | climate + spatial
#   6. Partial RDA: spatial | climate + salinity
#
# IMPORTANT:
# - The independent sampling units for these population-level
#   RDAs are populations, not individual RNA-seq samples.
# - Permutation tests assume exchangeability of populations.
#   If populations have a hierarchical or restricted sampling
#   design, adapt the permutation scheme accordingly.
# - Environmental variables are standardized before RDA.
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(vegan)
  library(adespatial)
  library(ade4)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(matrixStats)
})

# ============================================================
# 1. USER SETTINGS
# ============================================================

metadata_file <-
  "/Anastatica/metadata/Anastatica_RNA_metadata_tximeta.csv"

tximeta_input <-
  "/Anastatica/tximeta/tximeta_gene_level.rds"

output_dir <-
  "/Anastatica/RDA_transcriptome_analysis"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(12345)

n_permutations <- 9999
use_logCPM <- TRUE
scale_response_for_rda <- TRUE

count_assay_name <- NULL
count_assay_candidates <- c("counts", "count", "raw_counts", "NumReads")

climate_variables <- c("Climate_PC1", "Climate_PC2", "Climate_PC3")
salinity_variable <- "ECe_mean"

latitude_variable <- "Latitude"
longitude_variable <- "Longitude"

sample_variable <- "Sample_name"
population_variable <- "Population"
batch_variable <- "Sequencing_batch"

spatial_mem_autocorrelation <- "positive"
spatial_selection_method <- "FWD"
spatial_selection_alpha <- 0.05
spatial_selection_nperm <- 999
spatial_global_nperm <- 9999
max_spatial_predictors <- 6

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

check_columns <- function(data, columns, data_name = "data") {
  missing_columns <- setdiff(columns, colnames(data))
  if (length(missing_columns) > 0) {
    stop(
      "\nMissing columns in ", data_name, ":\n",
      paste(missing_columns, collapse = ", "), "\n"
    )
  }
}

check_numeric_variable <- function(x, variable_name) {
  if (!is.numeric(x)) {
    stop("\nVariable `", variable_name, "` must be numeric.")
  }
  if (all(is.na(x))) {
    stop("\nVariable `", variable_name, "` contains only NA values.")
  }
}

z_standardize <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("Cannot standardize a variable containing missing/non-finite values.")
  }
  s <- sd(x)
  if (!is.finite(s) || s == 0) {
    stop("Cannot standardize variable with zero or undefined variance.")
  }
  as.numeric(scale(x))
}

safe_adj_r2 <- function(rda_object) {
  result <- vegan::RsquareAdj(rda_object)
  tibble(
    R2 = unname(result$r.squared),
    adjusted_R2 = unname(result$adj.r.squared)
  )
}

extract_anova_p <- function(x) {
  if (is.null(x)) return(NA_real_)
  x <- as.data.frame(x)
  if (!"Pr(>F)" %in% colnames(x)) return(NA_real_)
  as.numeric(x$`Pr(>F)`[1])
}

extract_mem_p <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (!is.null(x$pvalue)) return(as.numeric(x$pvalue[1]))
  if (!is.null(x$Pvalue)) return(as.numeric(x$Pvalue[1]))
  NA_real_
}

cap_mem_matrix <- function(mem_matrix, max_predictors) {
  if (is.null(mem_matrix)) return(NULL)
  mem_matrix <- as.matrix(mem_matrix)
  if (ncol(mem_matrix) <= max_predictors) return(mem_matrix)
  mem_matrix[, seq_len(max_predictors), drop = FALSE]
}

write_rda_outputs <- function(model, prefix, output_dir,
                              n_permutations = 9999,
                              marginal = TRUE) {
  overall_test <- anova(model, permutations = n_permutations)
  model_r2 <- safe_adj_r2(model)

  write.csv(
    as.data.frame(overall_test),
    file.path(output_dir, paste0(prefix, "_overall_test.csv"))
  )
  write_csv(
    model_r2,
    file.path(output_dir, paste0(prefix, "_R2.csv"))
  )
  capture.output(
    summary(model),
    file = file.path(output_dir, paste0(prefix, "_summary.txt"))
  )

  term_test <- NULL
  if (marginal) {
    term_test <- anova(
      model, by = "margin", permutations = n_permutations
    )
    write.csv(
      as.data.frame(term_test),
      file.path(output_dir, paste0(prefix, "_marginal_tests.csv"))
    )
  }

  vif_values <- tryCatch(
    vegan::vif.cca(model),
    error = function(e) NULL
  )
  if (!is.null(vif_values)) {
    write.csv(
      data.frame(Predictor = names(vif_values), VIF = as.numeric(vif_values)),
      file.path(output_dir, paste0(prefix, "_VIF.csv")),
      row.names = FALSE
    )
  }

  list(
    overall_test = overall_test,
    term_test = term_test,
    r2 = model_r2,
    p = extract_anova_p(overall_test)
  )
}

# ============================================================
# 3. SOFTWARE VERSIONS
# ============================================================

package_names <- c(
  "SummarizedExperiment", "limma", "vegan", "adespatial", "ade4",
  "ggplot2", "readr", "dplyr", "tibble", "matrixStats"
)

package_versions <- tibble(
  package = c("R", package_names),
  version = c(
    as.character(getRversion()),
    vapply(package_names, function(x) as.character(packageVersion(x)),
           character(1))
  )
)

write_csv(package_versions, file.path(output_dir, "software_versions.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

# ============================================================
# 4. LOAD AND VALIDATE EXTERNAL METADATA
# ============================================================

metadata <- read_csv(metadata_file, show_col_types = FALSE) |> as.data.frame()

required_columns <- c(
  sample_variable, population_variable, batch_variable,
  latitude_variable, longitude_variable,
  climate_variables, salinity_variable
)
check_columns(metadata, required_columns, "metadata")

metadata <- metadata |>
  rename(
    Sample = all_of(sample_variable),
    Population = all_of(population_variable),
    Batch = all_of(batch_variable),
    Latitude = all_of(latitude_variable),
    Longitude = all_of(longitude_variable),
    ECe_mean = all_of(salinity_variable)
  )

metadata$Sample <- as.character(metadata$Sample)
metadata$Population <- as.character(metadata$Population)
metadata$Batch <- as.factor(metadata$Batch)

if (anyDuplicated(metadata$Sample) > 0) {
  stop("\nDuplicate sample IDs detected in external metadata.")
}

for (v in c("Latitude", "Longitude", climate_variables, "ECe_mean")) {
  check_numeric_variable(metadata[[v]], v)
}

if (any(!complete.cases(metadata[, c(
  "Sample", "Population", "Latitude", "Longitude",
  climate_variables, "ECe_mean", "Batch"
)]))) {
  stop(
    "\nMissing values detected in required metadata fields. ",
    "Resolve these explicitly before fitting the models."
  )
}

population_summary <- metadata |>
  group_by(Population) |>
  summarise(
    n_samples = n(),
    Batch_count = n_distinct(Batch),
    Latitude_mean = mean(Latitude),
    Longitude_mean = mean(Longitude),
    .groups = "drop"
  )

write_csv(
  population_summary,
  file.path(output_dir, "population_metadata_summary.csv")
)

# ============================================================
# 5. LOAD SUMMARIZED EXPERIMENT AND IDENTIFY COUNT ASSAY
# ============================================================

if (!file.exists(tximeta_input)) {
  stop("Could not find tximeta input: ", tximeta_input)
}

se <- readRDS(tximeta_input)
if (!inherits(se, "SummarizedExperiment")) {
  stop("\nRDS object is not a SummarizedExperiment.")
}

available_assays <- assayNames(se)
cat("Available assays:\n")
print(available_assays)

if (!is.null(count_assay_name)) {
  if (!count_assay_name %in% available_assays) {
    stop("Specified assay not found: ", count_assay_name)
  }
  selected_assay <- count_assay_name
} else {
  candidate_matches <- intersect(count_assay_candidates, available_assays)
  if (length(candidate_matches) == 0) {
    stop("\nNo raw-count assay identified. Set count_assay_name explicitly.")
  }
  selected_assay <- candidate_matches[1]
}

cat("Selected expression assay:", selected_assay, "\n")
expression_input <- as.matrix(assay(se, selected_assay))

if (!is.numeric(expression_input) ||
    anyNA(expression_input) ||
    any(!is.finite(expression_input)) ||
    any(expression_input < 0)) {
  stop("\nExpression assay must contain finite, non-negative numeric counts.")
}

# ============================================================
# 6. MATCH SAMPLES TO METADATA
# ============================================================

sample_meta <- as.data.frame(SummarizedExperiment::colData(se))
check_columns(
  sample_meta,
  c(sample_variable, population_variable, batch_variable),
  "tximeta sample metadata"
)

if (is.null(rownames(sample_meta)) ||
    !identical(rownames(sample_meta), colnames(expression_input))) {
  stop("\ntximeta colData row names do not match expression columns.")
}

sample_meta[[sample_variable]] <- as.character(sample_meta[[sample_variable]])
sample_meta[[population_variable]] <- as.character(sample_meta[[population_variable]])
sample_meta[[batch_variable]] <- as.character(sample_meta[[batch_variable]])

tximeta_sample_names <- sample_meta[[sample_variable]]
missing_external_metadata <- setdiff(tximeta_sample_names, metadata$Sample)

if (length(missing_external_metadata) > 0) {
  stop(
    "\nSamples missing from external metadata:\n",
    paste(missing_external_metadata, collapse = ", ")
  )
}

metadata_expression <- metadata[
  match(tximeta_sample_names, metadata$Sample), , drop = FALSE
]

if (!identical(metadata_expression$Sample, tximeta_sample_names)) {
  stop("\nFailed to synchronize metadata.")
}

# Ensure the population labels agree between tximeta and external metadata.
if (!identical(
  as.character(sample_meta[[population_variable]]),
  as.character(metadata_expression$Population)
)) {
  stop(
    "\nPopulation labels in tximeta colData and external metadata disagree."
  )
}

# ============================================================
# 7. FILTER ZERO-VARIANCE GENES AND TRANSFORM COUNTS
# ============================================================

gene_variance <- matrixStats::rowVars(expression_input)
keep_genes <- is.finite(gene_variance) & gene_variance > 0
expression_input <- expression_input[keep_genes, , drop = FALSE]

if (nrow(expression_input) == 0) stop("No variable genes remain after filtering.")

if (use_logCPM) {
  library_sizes <- colSums(expression_input)
  if (any(library_sizes <= 0)) stop("One or more samples have zero library size.")
  cpm_matrix <- sweep(expression_input, 2, library_sizes, "/") * 1e6
  expression_matrix <- log2(cpm_matrix + 1)
  transformation_description <- "log2(CPM + 1)"
} else {
  expression_matrix <- expression_input
  transformation_description <- "No transformation"
}

write.csv(
  expression_matrix,
  file.path(output_dir, "sample_level_expression_transformed.csv")
)

# ============================================================
# 8. BATCH ASSESSMENT / CORRECTION
# ============================================================

batch <- factor(metadata_expression$Batch)
population_factor <- factor(metadata_expression$Population)

batch_population_table <- table(metadata_expression$Population, batch)
write.csv(
  as.data.frame(batch_population_table),
  file.path(output_dir, "population_by_batch_table.csv"),
  row.names = FALSE
)

if (nlevels(batch) > 1) {
  design_preserve_population <- model.matrix(~ population_factor)
  expression_batch_corrected <- limma::removeBatchEffect(
    expression_matrix,
    batch = batch,
    design = design_preserve_population
  )
  cat("Batch correction performed while preserving population effects.\n")
} else {
  expression_batch_corrected <- expression_matrix
  cat("Only one batch detected; batch correction skipped.\n")
}

write.csv(
  expression_batch_corrected,
  file.path(output_dir, "sample_level_expression_batch_corrected.csv")
)

# ============================================================
# 9. POPULATION-LEVEL TRANSCRIPTOMIC PROFILES
# ============================================================

population_ids <- unique(metadata_expression$Population)

population_expression <- sapply(population_ids, function(pop) {
  sample_indices <- which(metadata_expression$Population == pop)
  rowMeans(expression_batch_corrected[, sample_indices, drop = FALSE])
})
population_expression <- as.matrix(population_expression)
colnames(population_expression) <- population_ids

write.csv(
  population_expression,
  file.path(output_dir, "population_level_expression.csv")
)

# ============================================================
# 10. TRANSCRIPTOMIC PCA (DESCRIPTIVE)
# ============================================================

transcriptomic_pca <- prcomp(
  t(population_expression),
  center = TRUE,
  scale. = TRUE
)

pca_variance <- transcriptomic_pca$sdev^2 / sum(transcriptomic_pca$sdev^2)
pca_scores <- as.data.frame(transcriptomic_pca$x)
pca_scores$Population <- rownames(pca_scores)
pca_scores <- left_join(pca_scores, population_summary, by = "Population")

write_csv(pca_scores, file.path(output_dir, "transcriptomic_PCA_scores.csv"))

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    title = "Population-level transcriptomic PCA",
    x = paste0("PC1 (", round(100 * pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_variance[2], 1), "%)")
  )

ggsave(
  file.path(output_dir, "transcriptomic_PCA.pdf"),
  pca_plot, width = 7, height = 6
)

# ============================================================
# 11. RESPONSE MATRIX FOR RDA
# ============================================================

if (scale_response_for_rda) {
  rda_response <- scale(t(population_expression), center = TRUE, scale = TRUE)
  rda_response_description <- "Population expression profiles standardized gene-wise"
} else {
  rda_response <- t(population_expression)
  rda_response_description <- "Population expression profiles not gene-standardized"
}

rda_response <- as.matrix(rda_response)
rownames(rda_response) <- population_ids

if (any(!is.finite(rda_response))) {
  stop("\nRDA response contains non-finite values.")
}

write.csv(
  rda_response,
  file.path(output_dir, "RDA_transcriptomic_response_matrix.csv")
)

# ============================================================
# 12. POPULATION-LEVEL ENVIRONMENTAL DATA
# ============================================================

population_environment <- metadata_expression |>
  group_by(Population) |>
  summarise(
    across(all_of(climate_variables), mean),
    ECe_mean = mean(ECe_mean),
    .groups = "drop"
  )

population_environment <- population_environment[
  match(population_ids, population_environment$Population), , drop = FALSE
]

if (!identical(population_environment$Population, population_ids)) {
  stop("Population environmental data are not aligned to response matrix.")
}

write_csv(
  population_environment,
  file.path(output_dir, "population_environment.csv")
)

environment_data <- population_environment |>
  mutate(
    across(all_of(climate_variables), z_standardize, .names = "{.col}_z"),
    Salinity_z = z_standardize(ECe_mean)
  ) |>
  as.data.frame()

rownames(environment_data) <- environment_data$Population

climate_predictor_names <- paste0(climate_variables, "_z")
salinity_predictor_names <- "Salinity_z"
environment_predictor_names <- c(climate_predictor_names, salinity_predictor_names)

climate_matrix <- as.matrix(environment_data[, climate_predictor_names, drop = FALSE])
salinity_matrix <- as.matrix(environment_data[, salinity_predictor_names, drop = FALSE])
environment_matrix <- as.matrix(environment_data[, environment_predictor_names, drop = FALSE])

write.csv(
  climate_matrix,
  file.path(output_dir, "climate_predictors_standardized.csv")
)
write.csv(
  salinity_matrix,
  file.path(output_dir, "salinity_predictor_standardized.csv")
)
write.csv(
  environment_matrix,
  file.path(output_dir, "climate_salinity_predictors_standardized.csv")
)

environment_correlations <- cor(environment_matrix)
write.csv(
  environment_correlations,
  file.path(output_dir, "environment_predictor_correlations.csv")
)

# Check whether climate PCs and salinity have sufficient independent variation.
if (qr(climate_matrix)$rank < ncol(climate_matrix)) {
  warning("Climate predictor matrix is rank deficient; inspect climate PCs.")
}
if (qr(cbind(climate_matrix, salinity_matrix))$rank <
    ncol(cbind(climate_matrix, salinity_matrix))) {
  warning("Climate + salinity predictor matrix is rank deficient.")
}

# ============================================================
# 13. GEOGRAPHIC DISTANCE AND dbMEMs
# ============================================================

haversine_distance <- function(lat1, lon1, lat2, lon2) {
  earth_radius_km <- 6371.0088
  lat1 <- lat1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon1 <- lon1 * pi / 180
  lon2 <- lon2 * pi / 180
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 +
    cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * earth_radius_km * atan2(sqrt(a), sqrt(1 - a))
}

calculate_geographic_distance <- function(coords) {
  n <- nrow(coords)
  distance_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(distance_matrix) <- coords$Population
  colnames(distance_matrix) <- coords$Population

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      distance_matrix[i, j] <- haversine_distance(
        coords$Latitude[i], coords$Longitude[i],
        coords$Latitude[j], coords$Longitude[j]
      )
    }
  }
  as.dist(distance_matrix)
}

population_coordinates <- metadata_expression |>
  group_by(Population) |>
  summarise(
    Latitude = mean(Latitude),
    Longitude = mean(Longitude),
    .groups = "drop"
  )

population_coordinates <- population_coordinates[
  match(population_ids, population_coordinates$Population), , drop = FALSE
]

if (any(!is.finite(population_coordinates$Latitude)) ||
    any(!is.finite(population_coordinates$Longitude))) {
  stop("Missing/non-finite geographic coordinates.")
}

geographic_distance <- calculate_geographic_distance(population_coordinates)
if (!identical(rownames(as.matrix(geographic_distance)), population_ids)) {
  stop("Geographic distance population order does not match response matrix.")
}

write.csv(
  as.matrix(geographic_distance),
  file.path(output_dir, "geographic_distance_matrix_km.csv")
)

geographic_distance_is_euclidean <- ade4::is.euclid(geographic_distance)

if (geographic_distance_is_euclidean) {
  geographic_distance_dbmem <- geographic_distance
  dbmem_distance_transformation <- "None; raw great-circle distance used"
} else {
  geographic_distance_dbmem <- ade4::cailliez(geographic_distance)
  dbmem_distance_transformation <- "Cailliez transformation applied for dbMEM construction only"
}

attr(geographic_distance_dbmem, "Labels") <- population_ids
write.csv(
  as.matrix(geographic_distance_dbmem),
  file.path(output_dir, "geographic_distance_matrix_for_dbMEM.csv")
)
writeLines(
  dbmem_distance_transformation,
  file.path(output_dir, "dbMEM_distance_transformation.txt")
)

dbmem_object <- adespatial::dbmem(
  geographic_distance_dbmem,
  thresh = NULL,
  MEM.autocor = spatial_mem_autocorrelation,
  store.listw = TRUE,
  silent = FALSE
)

dbmem_matrix <- as.matrix(as.data.frame(dbmem_object))
rownames(dbmem_matrix) <- population_ids

if (ncol(dbmem_matrix) == 0) stop("No positive dbMEMs were generated.")

write.csv(dbmem_matrix, file.path(output_dir, "all_positive_dbMEMs.csv"))

dbmem_eigenvalues <- attr(dbmem_object, "values")
write.csv(
  data.frame(
    dbMEM = colnames(dbmem_matrix),
    Eigenvalue = dbmem_eigenvalues[seq_len(ncol(dbmem_matrix))]
  ),
  file.path(output_dir, "dbMEM_eigenvalues.csv"),
  row.names = FALSE
)

# ============================================================
# 14. OVERALL SPATIAL MEM SELECTION
# ============================================================

mem_selection_overall <- adespatial::mem.select(
  x = rda_response,
  listw = attr(dbmem_object, "listw"),
  MEM.autocor = spatial_mem_autocorrelation,
  method = spatial_selection_method,
  MEM.all = FALSE,
  nperm = spatial_selection_nperm,
  nperm.global = spatial_global_nperm,
  alpha = spatial_selection_alpha,
  verbose = TRUE
)

overall_spatial_p <- extract_mem_p(mem_selection_overall$global.test)
overall_spatial_detected <- is.finite(overall_spatial_p) &&
  overall_spatial_p < spatial_selection_alpha

capture.output(
  print(mem_selection_overall$global.test),
  file = file.path(output_dir, "dbMEM_overall_global_test.txt")
)

if (overall_spatial_detected && !is.null(mem_selection_overall$MEM.select)) {
  overall_selected_mem_matrix <- as.matrix(mem_selection_overall$MEM.select)
  rownames(overall_selected_mem_matrix) <- population_ids
} else {
  overall_selected_mem_matrix <- matrix(
    nrow = nrow(rda_response), ncol = 0,
    dimnames = list(population_ids, character(0))
  )
}

overall_selected_mem_matrix <- cap_mem_matrix(
  overall_selected_mem_matrix, max_spatial_predictors
)

overall_selected_mem_names <- colnames(overall_selected_mem_matrix)

# Handle the case where no overall MEMs were selected.
if (is.null(overall_selected_mem_names)) {
  overall_selected_mem_names <- character(0)
}

write.csv(
  overall_selected_mem_matrix,
  file.path(output_dir, "overall_selected_dbMEMs.csv")
)
writeLines(
  overall_selected_mem_names,
  file.path(output_dir, "overall_selected_dbMEM_names.txt")
)
if (!is.null(mem_selection_overall$summary)) {
  write.csv(
    as.data.frame(mem_selection_overall$summary),
    file.path(output_dir, "overall_dbMEM_selection_summary.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 15. CLIMATE- AND SALINITY-CONDITIONED SPATIAL MEM SELECTION
# ============================================================

# This diagnostic asks whether positive spatial structure remains
# after accounting for climate PCs + salinity jointly.
environment_rda_for_residuals <- vegan::rda(
  rda_response ~ Climate_PC1_z + Climate_PC2_z +
    Climate_PC3_z + Salinity_z,
  data = environment_data
)

environment_residuals <- as.matrix(residuals(environment_rda_for_residuals))
rownames(environment_residuals) <- population_ids

mem_selection_conditional <- adespatial::mem.select(
  x = environment_residuals,
  listw = attr(dbmem_object, "listw"),
  MEM.autocor = spatial_mem_autocorrelation,
  method = spatial_selection_method,
  MEM.all = FALSE,
  nperm = spatial_selection_nperm,
  nperm.global = spatial_global_nperm,
  alpha = spatial_selection_alpha,
  verbose = TRUE
)

conditional_spatial_p <- extract_mem_p(mem_selection_conditional$global.test)
conditional_spatial_detected <- is.finite(conditional_spatial_p) &&
  conditional_spatial_p < spatial_selection_alpha

capture.output(
  print(mem_selection_conditional$global.test),
  file = file.path(output_dir, "dbMEM_conditional_global_test.txt")
)

if (conditional_spatial_detected &&
    !is.null(mem_selection_conditional$MEM.select)) {
  conditional_selected_mem_matrix <- as.matrix(mem_selection_conditional$MEM.select)
  rownames(conditional_selected_mem_matrix) <- population_ids
} else {
  conditional_selected_mem_matrix <- matrix(
    nrow = nrow(rda_response), ncol = 0,
    dimnames = list(population_ids, character(0))
  )
}


conditional_selected_mem_matrix <- cap_mem_matrix(
  conditional_selected_mem_matrix, max_spatial_predictors
)

conditional_selected_mem_names <- colnames(conditional_selected_mem_matrix)

# Handle the case where no conditional MEMs were selected.
if (is.null(conditional_selected_mem_names)) {
  conditional_selected_mem_names <- character(0)
}


write.csv(
  conditional_selected_mem_matrix,
  file.path(output_dir, "conditional_selected_dbMEMs.csv")
)
writeLines(
  conditional_selected_mem_names,
  file.path(output_dir, "conditional_selected_dbMEM_names.txt")
)
if (!is.null(mem_selection_conditional$summary)) {
  write.csv(
    as.data.frame(mem_selection_conditional$summary),
    file.path(output_dir, "conditional_dbMEM_selection_summary.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 16. FIT THE PRIMARY SINGLE-PREDICTOR-SET RDAs
# ============================================================

# Climate-only RDA: Climate_PC1 + Climate_PC2 + Climate_PC3.
climate_rda <- vegan::rda(
  rda_response ~ Climate_PC1_z + Climate_PC2_z + Climate_PC3_z,
  data = environment_data
)
climate_results <- write_rda_outputs(
  climate_rda, "RDA_climate_only", output_dir, n_permutations
)

# Salinity-only RDA.
salinity_rda <- vegan::rda(
  rda_response ~ Salinity_z,
  data = environment_data
)
salinity_results <- write_rda_outputs(
  salinity_rda, "RDA_salinity_only", output_dir, n_permutations
)

# Spatial-only RDA using selected overall positive dbMEMs.
if (length(overall_selected_mem_names) > 0) {
  overall_spatial_data <- as.data.frame(overall_selected_mem_matrix)
  spatial_rda_overall <- vegan::rda(rda_response ~ ., data = overall_spatial_data)
  spatial_results <- write_rda_outputs(
    spatial_rda_overall, "RDA_spatial_overall", output_dir, n_permutations
  )
  spatial_rda_overall_r2 <- spatial_results$r2
  spatial_rda_overall_anova <- spatial_results$overall_test
} else {
  overall_spatial_data <- data.frame(row.names = population_ids)
  spatial_rda_overall <- NULL
  spatial_rda_overall_r2 <- tibble(R2 = NA_real_, adjusted_R2 = NA_real_)
  spatial_rda_overall_anova <- NULL
}

# Data frame combining environmental predictors and selected overall
# spatial MEMs, used by the partial RDAs below.
model_data <- cbind(
  environment_data,
  as.data.frame(overall_selected_mem_matrix)
)

# ============================================================
# 17. PARTIAL RDAs
# ============================================================

# Climate conditional on salinity + spatial structure.
if (length(overall_selected_mem_names) > 0) {
  climate_given_salinity_space_formula <- as.formula(
    paste(
      "rda_response ~",
      paste(climate_predictor_names, collapse = " + "),
      "+ Condition(",
      paste(c(salinity_predictor_names, overall_selected_mem_names),
            collapse = " + "),
      ")"
    )
  )
} else {
  climate_given_salinity_space_formula <- as.formula(
    paste(
      "rda_response ~",
      paste(climate_predictor_names, collapse = " + "),
      "+ Condition(Salinity_z)"
    )
  )
}

climate_given_salinity_space_rda <- vegan::rda(
  climate_given_salinity_space_formula, data = model_data
)
climate_given_salinity_space_results <- write_rda_outputs(
  climate_given_salinity_space_rda,
  "partial_RDA_climate_given_salinity_and_space",
  output_dir, n_permutations
)

# Salinity conditional on climate + spatial structure.
if (length(overall_selected_mem_names) > 0) {
  salinity_given_climate_space_formula <- as.formula(
    paste(
      "rda_response ~ Salinity_z + Condition(",
      paste(c(climate_predictor_names, overall_selected_mem_names),
            collapse = " + "),
      ")"
    )
  )
} else {
  salinity_given_climate_space_formula <- as.formula(
    paste(
      "rda_response ~ Salinity_z + Condition(",
      paste(climate_predictor_names, collapse = " + "),
      ")"
    )
  )
}

salinity_given_climate_space_rda <- vegan::rda(
  salinity_given_climate_space_formula, data = model_data
)
salinity_given_climate_space_results <- write_rda_outputs(
  salinity_given_climate_space_rda,
  "partial_RDA_salinity_given_climate_and_space",
  output_dir, n_permutations
)

# Spatial conditional on climate + salinity.
if (length(overall_selected_mem_names) > 0) {
  spatial_given_environment_formula <- as.formula(
    paste(
      "rda_response ~",
      paste(overall_selected_mem_names, collapse = " + "),
      "+ Condition(",
      paste(environment_predictor_names, collapse = " + "),
      ")"
    )
  )

  spatial_given_environment_rda <- vegan::rda(
    spatial_given_environment_formula, data = model_data
  )
  spatial_given_environment_results <- write_rda_outputs(
    spatial_given_environment_rda,
    "partial_RDA_space_given_climate_salinity",
    output_dir, n_permutations
  )
} else {
  spatial_given_environment_rda <- NULL
  spatial_given_environment_results <- list(
    overall_test = NULL,
    term_test = NULL,
    r2 = tibble(R2 = NA_real_, adjusted_R2 = NA_real_),
    p = NA_real_
  )
}

# ============================================================
# 18. VARIATION PARTITIONING
# ============================================================

# Three explanatory tables: climate, salinity, and spatial MEMs.
# This separates unique/shared fractions among these three groups.
if (length(overall_selected_mem_names) > 0) {
  variation_partitioning <- vegan::varpart(
    rda_response,
    climate_matrix,
    salinity_matrix,
    overall_selected_mem_matrix,
    scale = FALSE
  )

  capture.output(
    print(variation_partitioning),
    file = file.path(output_dir, "variation_partitioning_climate_salinity_space.txt")
  )
  capture.output(
    summary(variation_partitioning),
    file = file.path(output_dir, "variation_partitioning_climate_salinity_space_summary.txt")
  )

  write.csv(
    variation_partitioning$part$indfract,
    file.path(output_dir, "variation_partitioning_climate_salinity_space_individual_fractions.csv"),
    row.names = TRUE
  )
  write.csv(
    variation_partitioning$part$fract,
    file.path(output_dir, "variation_partitioning_climate_salinity_space_fractions.csv"),
    row.names = TRUE
  )

  pdf(
    file.path(output_dir, "variation_partitioning_climate_salinity_space.pdf"),
    width = 9, height = 7
  )
  plot(variation_partitioning)
  dev.off()
} else {
  variation_partitioning <- NULL
  writeLines(
    "Variation partitioning skipped because no overall positive spatial MEMs were selected.",
    file.path(output_dir, "variation_partitioning_skipped.txt")
  )
}

# ============================================================
# 19. SUMMARY TABLES
# ============================================================

primary_results <- tibble(
  analysis = c(
    "Climate only (PC1 + PC2 + PC3)",
    "Salinity only",
    "Spatial only (selected overall dbMEMs)",
    "Climate | salinity + spatial",
    "Salinity | climate + spatial",
    "Spatial | climate + salinity"
  ),
  R2 = c(
    climate_results$r2$R2,
    salinity_results$r2$R2,
    spatial_rda_overall_r2$R2,
    climate_given_salinity_space_results$r2$R2,
    salinity_given_climate_space_results$r2$R2,
    spatial_given_environment_results$r2$R2
  ),
  adjusted_R2 = c(
    climate_results$r2$adjusted_R2,
    salinity_results$r2$adjusted_R2,
    spatial_rda_overall_r2$adjusted_R2,
    climate_given_salinity_space_results$r2$adjusted_R2,
    salinity_given_climate_space_results$r2$adjusted_R2,
    spatial_given_environment_results$r2$adjusted_R2
  ),
  P_value = c(
    climate_results$p,
    salinity_results$p,
    extract_anova_p(spatial_rda_overall_anova),
    climate_given_salinity_space_results$p,
    salinity_given_climate_space_results$p,
    spatial_given_environment_results$p
  )
)

write_csv(primary_results, file.path(output_dir, "primary_RDA_results.csv"))

spatial_interpretation <- tibble(
  overall_spatial_structure_detected = overall_spatial_detected,
  overall_spatial_global_P = overall_spatial_p,
  overall_MEMs_selected = length(overall_selected_mem_names),
  climate_salinity_conditioned_spatial_structure_detected = conditional_spatial_detected,
  climate_salinity_conditioned_spatial_global_P = conditional_spatial_p,
  conditional_MEMs_selected = length(conditional_selected_mem_names),
  interpretation = if (
    overall_spatial_detected && !conditional_spatial_detected
  ) {
    "Overall positive spatial structure detected, but no significant positive spatial structure was detected in residuals after climate + salinity were fitted."
  } else if (
    overall_spatial_detected && conditional_spatial_detected
  ) {
    "Positive spatial structure detected both overall and in residuals after climate + salinity were fitted."
  } else {
    "No significant overall positive spatial structure detected, or inspect saved global test outputs."
  }
)

write_csv(
  spatial_interpretation,
  file.path(output_dir, "spatial_interpretation_summary.csv")
)

analysis_summary <- tibble(
  n_samples = nrow(metadata_expression),
  n_populations = length(population_ids),
  n_genes = nrow(population_expression),
  expression_transformation = transformation_description,
  response_standardized_gene_wise = scale_response_for_rda,
  n_dbMEMs_generated = ncol(dbmem_matrix),
  n_overall_MEMs_selected = length(overall_selected_mem_names),
  n_conditional_MEMs_selected = length(conditional_selected_mem_names),
  overall_spatial_global_P = overall_spatial_p,
  conditional_spatial_global_P = conditional_spatial_p,
  climate_only_P = climate_results$p,
  climate_only_adjusted_R2 = climate_results$r2$adjusted_R2,
  salinity_only_P = salinity_results$p,
  salinity_only_adjusted_R2 = salinity_results$r2$adjusted_R2,
  spatial_only_P = extract_anova_p(spatial_rda_overall_anova),
  spatial_only_adjusted_R2 = spatial_rda_overall_r2$adjusted_R2,
  climate_given_salinity_space_P = climate_given_salinity_space_results$p,
  climate_given_salinity_space_adjusted_R2 = climate_given_salinity_space_results$r2$adjusted_R2,
  salinity_given_climate_space_P = salinity_given_climate_space_results$p,
  salinity_given_climate_space_adjusted_R2 = salinity_given_climate_space_results$r2$adjusted_R2
)

write_csv(analysis_summary, file.path(output_dir, "analysis_summary.csv"))

# ============================================================
# 20. SAVE COMPLETE ANALYSIS OBJECTS
# ============================================================

saveRDS(
  list(
    population_expression = population_expression,
    rda_response = rda_response,
    population_environment = population_environment,
    environment_data = environment_data,
    climate_matrix = climate_matrix,
    salinity_matrix = salinity_matrix,
    environment_matrix = environment_matrix,
    geographic_distance = geographic_distance,
    dbmem_object = dbmem_object,
    dbmem_matrix = dbmem_matrix,
    mem_selection_overall = mem_selection_overall,
    overall_selected_mem_matrix = overall_selected_mem_matrix,
    mem_selection_conditional = mem_selection_conditional,
    conditional_selected_mem_matrix = conditional_selected_mem_matrix,
    climate_rda = climate_rda,
    salinity_rda = salinity_rda,
    spatial_rda_overall = spatial_rda_overall,
    climate_given_salinity_space_rda = climate_given_salinity_space_rda,
    salinity_given_climate_space_rda = salinity_given_climate_space_rda,
    spatial_given_environment_rda = spatial_given_environment_rda,
    variation_partitioning = variation_partitioning
  ),
  file.path(output_dir, "RDA_complete_analysis_objects.rds")
)

# ============================================================
# 21. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("RDA TRANSCRIPTOME ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Output directory: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
cat("Samples: ", nrow(metadata_expression), "\n", sep = "")
cat("Populations: ", length(population_ids), "\n", sep = "")
cat("Genes: ", nrow(population_expression), "\n", sep = "")
cat("Climate-only adjusted R2: ",
    round(climate_results$r2$adjusted_R2, 4), "\n", sep = "")
cat("Climate-only P: ",
    format.pval(climate_results$p, digits = 4), "\n", sep = "")
cat("Salinity-only adjusted R2: ",
    round(salinity_results$r2$adjusted_R2, 4), "\n", sep = "")
cat("Salinity-only P: ",
    format.pval(salinity_results$p, digits = 4), "\n", sep = "")
cat("Overall spatial MEMs selected: ",
    length(overall_selected_mem_names), "\n", sep = "")
cat("Conditional spatial MEMs selected: ",
    length(conditional_selected_mem_names), "\n", sep = "")
cat("============================================================\n")
