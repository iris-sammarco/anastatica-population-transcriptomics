#!/usr/bin/env Rscript

# ============================================================
# Anastatica hierochuntica
#
# POPULATION-LEVEL TRANSCRIPTOMIC ASSOCIATION WITH
# GEOGRAPHY AND ENVIRONMENT
#
# ============================================================
#
# AUTHOR: Iris Sammarco
#
# DATE: 25/09/2026
#
# ============================================================
#
# BIOLOGICAL QUESTION
# -------------------
#
# Are populations that are geographically close more
# transcriptomically similar, and/or are populations from
# environmentally similar locations more transcriptomically
# similar?
#
# ============================================================
#
# PRIMARY ANALYTICAL FRAMEWORK
# ----------------------------
#
# The primary inferential framework is:
#
#   1. Population-level RDA
#   2. Spatial dbMEMs
#   3. Forward selection of spatial MEMs
#   4. Partial RDA
#   5. Variation partitioning
#
# ============================================================
#
# SPATIAL ANALYSIS
# ----------------
#
# Two spatial analyses are performed:
#
# A. OVERALL SPATIAL STRUCTURE
#
#    MEM selection is performed on the transcriptomic response
#    matrix itself.
#
#    This asks:
#
#       "Is transcriptomic variation spatially structured at all?"
#
#
# B. ENVIRONMENT-CONDITIONED SPATIAL STRUCTURE
#
#    MEM selection is performed on transcriptomic residuals
#    after climate + salinity are removed.
#
#    This asks:
#
#       "Is there spatial structure remaining after accounting
#        for measured environmental predictors?"
#
# This distinction is important because environment itself can
# be spatially structured.
#
# ============================================================
#
# IMPORTANT INTERPRETATION
# ------------------------
#
# If:
#
#   overall spatial structure = significant
#
# but:
#
#   environment-conditioned spatial structure = not significant
#
# then the data are consistent with transcriptomic spatial
# structure being explained largely by environmental structure,
# rather than independent residual geography.
#
# If both are significant, spatial structure remains after
# accounting for measured environment.
#
# If neither is significant, there is no detectable positive
# spatial structure under this model.
#
# ============================================================


# ============================================================
# 0. PACKAGES
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
  "/Anastatica/IBD_IBE_transcriptome_analysis"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

set.seed(12345)


# ------------------------------------------------------------
# Permutations
# ------------------------------------------------------------

n_permutations <- 9999


# ------------------------------------------------------------
# Expression transformation
# ------------------------------------------------------------

use_logCPM <- TRUE


# ------------------------------------------------------------
# Transcriptomic distance
# ------------------------------------------------------------

transcriptomic_distance_method <- "euclidean"


# ------------------------------------------------------------
# Standardize genes before RDA
# ------------------------------------------------------------

scale_response_for_rda <- TRUE


# ------------------------------------------------------------
# Expression assay
# ------------------------------------------------------------

count_assay_name <- NULL

count_assay_candidates <- c(
  "counts",
  "count",
  "raw_counts",
  "NumReads"
)


# ------------------------------------------------------------
# Environmental predictors
# ------------------------------------------------------------

climate_variables <- c(
  "Climate_PC1",
  "Climate_PC2",
  "Climate_PC3"
)

salinity_variable <- "ECe_mean"


# ------------------------------------------------------------
# Geographic coordinates
# ------------------------------------------------------------

latitude_variable <- "Latitude"
longitude_variable <- "Longitude"


# ------------------------------------------------------------
# Metadata variables
# ------------------------------------------------------------

sample_variable <- "Sample_name"
population_variable <- "Population"
batch_variable <- "Sequencing_batch"


# ------------------------------------------------------------
# Spatial settings
# ------------------------------------------------------------

spatial_mem_autocorrelation <- "positive"

spatial_selection_method <- "FWD"

spatial_selection_alpha <- 0.05

spatial_selection_nperm <- 999

spatial_global_nperm <- 9999


# ------------------------------------------------------------
# Maximum number of spatial predictors
# ------------------------------------------------------------

max_spatial_predictors <- 6


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================


check_columns <- function(
    data,
    columns,
    data_name = "data"
) {

  missing_columns <-
    setdiff(
      columns,
      colnames(data)
    )

  if (length(missing_columns) > 0) {

    stop(
      paste0(
        "\nMissing columns in ",
        data_name,
        ":\n",
        paste(
          missing_columns,
          collapse = ", "
        ),
        "\n"
      )
    )

  }

}


check_numeric_variable <- function(
    x,
    variable_name
) {

  if (!is.numeric(x)) {

    stop(
      paste0(
        "\nVariable `",
        variable_name,
        "` must be numeric."
      )
    )

  }

  if (all(is.na(x))) {

    stop(
      paste0(
        "\nVariable `",
        variable_name,
        "` contains only NA values."
      )
    )

  }

}


z_standardize <- function(x) {

  x <- as.numeric(x)

  s <- sd(
    x,
    na.rm = TRUE
  )

  if (!is.finite(s) || s == 0) {

    stop(
      "\nCannot standardize variable with zero or undefined variance."
    )

  }

  as.numeric(
    scale(x)
  )

}

distance_pair_table <- function(
    response_distance,
    predictor_distance,
    response_name,
    predictor_name
) {

  response_matrix <-
    as.matrix(
      response_distance
    )

  predictor_matrix <-
    as.matrix(
      predictor_distance
    )

  if (!identical(
    rownames(response_matrix),
    rownames(predictor_matrix)
  )) {

    stop(
      "\nDistance matrices do not have identical population order."
    )

  }

  idx <-
    upper.tri(
      response_matrix,
      diag = FALSE
    )

  tibble(

    !!response_name :=
      response_matrix[idx],

    !!predictor_name :=
      predictor_matrix[idx]

  )

}


haversine_distance <- function(
    lat1,
    lon1,
    lat2,
    lon2
) {

  earth_radius_km <- 6371.0088

  lat1 <- lat1 * pi / 180
  lat2 <- lat2 * pi / 180

  lon1 <- lon1 * pi / 180
  lon2 <- lon2 * pi / 180

  dlat <- lat2 - lat1
  dlon <- lon2 - lon1

  a <-
    sin(dlat / 2)^2 +
    cos(lat1) *
    cos(lat2) *
    sin(dlon / 2)^2

  c <- 2 * atan2(
    sqrt(a),
    sqrt(1 - a)
  )

  earth_radius_km * c

}


calculate_geographic_distance <- function(
    coords
) {

  n <- nrow(coords)

  distance_matrix <-
    matrix(
      0,
      nrow = n,
      ncol = n
    )

  rownames(distance_matrix) <-
    coords$Population

  colnames(distance_matrix) <-
    coords$Population

  for (i in seq_len(n)) {

    for (j in seq_len(n)) {

      distance_matrix[i, j] <-
        haversine_distance(
          coords$Latitude[i],
          coords$Longitude[i],
          coords$Latitude[j],
          coords$Longitude[j]
        )

    }

  }

  as.dist(
    distance_matrix
  )

}


validate_distance <- function(
    d,
    expected_names
) {

  if (!inherits(
    d,
    "dist"
  )) {

    stop(
      "\nObject is not a dist object."
    )

  }

  m <- as.matrix(d)

  if (!identical(
    rownames(m),
    expected_names
  )) {

    stop(
      "\nDistance matrix population names/order do not match."
    )

  }

  if (any(!is.finite(m))) {

    stop(
      "\nDistance matrix contains non-finite values."
    )

  }

}


safe_adj_r2 <- function(
    rda_object
) {

  result <-
    vegan::RsquareAdj(
      rda_object
    )

  tibble(

    R2 =
      unname(
        result$r.squared
      ),

    adjusted_R2 =
      unname(
        result$adj.r.squared
      )

  )

}


extract_anova_p <- function(
    x
) {

  if (is.null(x)) {

    return(
      NA_real_
    )

  }

  x <-
    as.data.frame(x)

  if (!"Pr(>F)" %in% colnames(x)) {

    return(
      NA_real_
    )

  }

  as.numeric(
    x$`Pr(>F)`[1]
  )

}


extract_mem_p <- function(
    x
) {

  if (is.null(x)) {

    return(
      NA_real_
    )

  }

  if (!is.null(x$pvalue)) {

    return(
      as.numeric(
        x$pvalue[1]
      )
    )

  }

  if (!is.null(x$Pvalue)) {

    return(
      as.numeric(
        x$Pvalue[1]
      )
    )

  }

  NA_real_

}


cap_mem_matrix <- function(
    mem_matrix,
    max_predictors
) {

  if (is.null(mem_matrix)) {

    return(
      NULL
    )

  }

  mem_matrix <-
    as.matrix(
      mem_matrix
    )

  if (ncol(mem_matrix) == 0) {

    return(
      mem_matrix
    )

  }

  if (ncol(mem_matrix) <= max_predictors) {

    return(
      mem_matrix
    )

  }

  mem_matrix[
    ,
    seq_len(max_predictors),
    drop = FALSE
  ]

}


# ============================================================
# 3. SOFTWARE VERSIONS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Saving software versions\n"
)

cat(
  "============================================\n\n"
)


package_names <- c(

  "SummarizedExperiment",
  "limma",
  "vegan",
  "adespatial",
  "ade4",
  "ggplot2",
  "readr",
  "dplyr",
  "tibble",
  "matrixStats"

)


package_versions <-
  tibble(

    package =
      c(
        "R",
        package_names
      ),

    version =
      c(

        as.character(
          getRversion()
        ),

        sapply(
          package_names,
          function(x) {

            as.character(
              packageVersion(x)
            )

          }
        )

      )

  )


write_csv(
  package_versions,
  file.path(
    output_dir,
    "software_versions.csv"
  )
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      output_dir,
      "sessionInfo.txt"
    )
)


# ============================================================
# 4. LOAD EXTERNAL METADATA
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Loading metadata\n"
)

cat(
  "============================================\n\n"
)


metadata <-
  read_csv(
    metadata_file,
    show_col_types = FALSE
  )


metadata <-
  as.data.frame(
    metadata
  )


required_columns <- c(

  sample_variable,
  population_variable,
  batch_variable,

  latitude_variable,
  longitude_variable,

  climate_variables,

  salinity_variable

)


check_columns(
  metadata,
  required_columns,
  "metadata"
)


metadata <-
  metadata %>%

  rename(

    Sample =
      all_of(sample_variable),

    Population =
      all_of(population_variable),

    Batch =
      all_of(batch_variable),

    Latitude =
      all_of(latitude_variable),

    Longitude =
      all_of(longitude_variable),

    ECe_mean =
      all_of(salinity_variable)

  )


metadata$Sample <-
  as.character(
    metadata$Sample
  )


metadata$Population <-
  as.character(
    metadata$Population
  )


metadata$Batch <-
  as.factor(
    metadata$Batch
  )


if (anyDuplicated(
  metadata$Sample
) > 0) {

  stop(
    "\nDuplicate sample IDs detected."
  )

}


check_numeric_variable(
  metadata$Latitude,
  "Latitude"
)

check_numeric_variable(
  metadata$Longitude,
  "Longitude"
)

for (v in climate_variables) {

  check_numeric_variable(
    metadata[[v]],
    v
  )

}

check_numeric_variable(
  metadata$ECe_mean,
  "ECe_mean"
)


population_summary <-
  metadata %>%

  group_by(
    Population
  ) %>%

  summarise(

    n_samples =
      n(),

    Batch_count =
      n_distinct(
        Batch
      ),

    Latitude_mean =
      mean(
        Latitude,
        na.rm = TRUE
      ),

    Longitude_mean =
      mean(
        Longitude,
        na.rm = TRUE
      ),

    .groups = "drop"

  )


write_csv(
  population_summary,
  file.path(
    output_dir,
    "population_metadata_summary.csv"
  )
)


# ============================================================
# 5. LOAD SUMMARIZED EXPERIMENT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Loading transcriptomic data\n"
)

cat(
  "============================================\n\n"
)


if (!file.exists(
  tximeta_input
)) {

  stop(
    paste0(
      "\nCould not find:\n",
      tximeta_input
    )
  )

}


se <-
  readRDS(
    tximeta_input
  )


if (!inherits(
  se,
  "SummarizedExperiment"
)) {

  stop(
    "\nRDS object is not a SummarizedExperiment."
  )

}


cat(
  "Genes/features:",
  nrow(se),
  "\n"
)

cat(
  "Samples:",
  ncol(se),
  "\n\n"
)


# ============================================================
# 6. IDENTIFY EXPRESSION ASSAY
# ============================================================

available_assays <-
  assayNames(
    se
  )


cat(
  "Available assays:\n"
)

print(
  available_assays
)

cat(
  "\n"
)


if (!is.null(
  count_assay_name
)) {

  if (!count_assay_name %in%
      available_assays) {

    stop(
      paste0(
        "\nSpecified assay not found: ",
        count_assay_name
      )
    )

  }

  selected_assay <-
    count_assay_name

} else {

  candidate_matches <-
    intersect(
      count_assay_candidates,
      available_assays
    )

  if (length(
    candidate_matches
  ) == 0) {

    stop(
      "\nNo raw-count assay identified."
    )

  }

  selected_assay <-
    candidate_matches[1]

}


cat(
  "Selected expression assay:",
  selected_assay,
  "\n\n"
)


expression_input <-
  assay(
    se,
    selected_assay
  )


expression_input <-
  as.matrix(
    expression_input
  )


if (!is.numeric(
  expression_input
)) {

  stop(
    "\nExpression assay is not numeric."
  )

}


if (anyNA(
  expression_input
)) {

  stop(
    "\nExpression assay contains NA values."
  )

}


if (any(
  expression_input < 0
)) {

  stop(
    "\nExpression assay contains negative values."
  )

}


# ============================================================
# 7. MATCH TXIMETA SAMPLES TO EXTERNAL METADATA
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Matching tximeta samples to external metadata\n"
)

cat(
  "============================================\n\n"
)


sample_meta <-
  as.data.frame(
    SummarizedExperiment::colData(
      se
    )
  )


check_columns(
  sample_meta,
  c(
    sample_variable,
    population_variable,
    batch_variable
  ),
  "tximeta sample metadata"
)


if (is.null(
  rownames(sample_meta)
)) {

  stop(
    "\ntximeta colData has no row names."
  )

}


if (!identical(
  rownames(sample_meta),
  colnames(expression_input)
)) {

  stop(
    "\ntximeta colData row names do not match expression columns."
  )

}


sample_meta[[sample_variable]] <-
  as.character(
    sample_meta[[sample_variable]]
  )


sample_meta[[population_variable]] <-
  as.character(
    sample_meta[[population_variable]]
  )


sample_meta[[batch_variable]] <-
  as.character(
    sample_meta[[batch_variable]]
  )


tximeta_sample_names <-
  sample_meta[[sample_variable]]


missing_external_metadata <-
  setdiff(
    tximeta_sample_names,
    metadata$Sample
  )


if (length(
  missing_external_metadata
) > 0) {

  stop(
    paste0(
      "\nSamples missing from external metadata:\n",
      paste(
        missing_external_metadata,
        collapse = ", "
      )
    )
  )

}


metadata_expression <-
  metadata[
    match(
      tximeta_sample_names,
      metadata$Sample
    ),
    ,
    drop = FALSE
  ]


if (!identical(
  metadata_expression$Sample,
  tximeta_sample_names
)) {

  stop(
    "\nFailed to synchronize metadata."
  )

}


# ============================================================
# 8. FILTER ZERO-VARIANCE GENES
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Filtering genes\n"
)

cat(
  "============================================\n\n"
)


gene_variance <-
  matrixStats::rowVars(
    expression_input
  )


keep_genes <-
  is.finite(
    gene_variance
  ) &
  gene_variance > 0


expression_input <-
  expression_input[
    keep_genes,
    ,
    drop = FALSE
  ]


cat(
  "Genes retained:",
  nrow(expression_input),
  "\n"
)


# ============================================================
# 9. EXPRESSION TRANSFORMATION
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Expression transformation\n"
)

cat(
  "============================================\n\n"
)


if (use_logCPM) {

  library_sizes <-
    colSums(
      expression_input
    )

  if (any(
    library_sizes <= 0
  )) {

    stop(
      "\nOne or more samples have zero library size."
    )

  }

  cpm_matrix <-
    sweep(
      expression_input,
      MARGIN = 2,
      STATS = library_sizes,
      FUN = "/"
    ) * 1e6

  expression_matrix <-
    log2(
      cpm_matrix + 1
    )

  transformation_description <-
    "log2(CPM + 1)"

} else {

  expression_matrix <-
    expression_input

  transformation_description <-
    "No transformation"

}


cat(
  "Expression transformation:",
  transformation_description,
  "\n"
)


write.csv(
  expression_matrix,
  file.path(
    output_dir,
    "sample_level_expression_transformed.csv"
  )
)


# ============================================================
# 10. BATCH ASSESSMENT / CORRECTION
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Batch assessment / correction\n"
)

cat(
  "============================================\n\n"
)


batch <-
  factor(
    metadata_expression$Batch
  )


population_factor <-
  factor(
    metadata_expression$Population
  )


batch_population_table <-
  table(
    metadata_expression$Population,
    batch
  )


write.csv(
  as.data.frame(
    batch_population_table
  ),
  file.path(
    output_dir,
    "population_by_batch_table.csv"
  ),
  row.names = FALSE
)


if (nlevels(batch) > 1) {

  design_preserve_population <-
    model.matrix(
      ~ population_factor
    )

  expression_batch_corrected <-
    limma::removeBatchEffect(
      expression_matrix,
      batch = batch,
      design = design_preserve_population
    )

  cat(
    "Batch correction performed while preserving population effects.\n"
  )

} else {

  expression_batch_corrected <-
    expression_matrix

  cat(
    "Only one batch detected; batch correction skipped.\n"
  )

}


write.csv(
  expression_batch_corrected,
  file.path(
    output_dir,
    "sample_level_expression_batch_corrected.csv"
  )
)


# ============================================================
# 11. POPULATION-LEVEL TRANSCRIPTOMIC PROFILES
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Creating population-level transcriptomic profiles\n"
)

cat(
  "============================================\n\n"
)


population_ids <-
  unique(
    metadata_expression$Population
  )


population_expression <-
  sapply(

    population_ids,

    function(pop) {

      sample_indices <-
        which(
          metadata_expression$Population ==
            pop
        )

      rowMeans(
        expression_batch_corrected[
          ,
          sample_indices,
          drop = FALSE
        ]
      )

    }

  )


population_expression <-
  as.matrix(
    population_expression
  )


colnames(
  population_expression
) <-
  population_ids


cat(
  "Number of populations:",
  ncol(population_expression),
  "\n"
)


write.csv(
  population_expression,
  file.path(
    output_dir,
    "population_level_expression.csv"
  )
)


# ============================================================
# 12. TRANSCRIPTOMIC PCA
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Transcriptomic PCA\n"
)

cat(
  "============================================\n\n"
)


transcriptomic_pca <-
  prcomp(
    t(
      population_expression
    ),
    center = TRUE,
    scale. = TRUE
  )


pca_variance <-
  transcriptomic_pca$sdev^2 /
  sum(
    transcriptomic_pca$sdev^2
  )


pca_scores <-
  as.data.frame(
    transcriptomic_pca$x
  )


pca_scores$Population <-
  rownames(
    pca_scores
  )


pca_scores <-
  left_join(
    pca_scores,
    population_summary,
    by = "Population"
  )


write_csv(
  pca_scores,
  file.path(
    output_dir,
    "transcriptomic_PCA_scores.csv"
  )
)


pca_plot <-
  ggplot(
    pca_scores,
    aes(
      x = PC1,
      y = PC2
    )
  ) +

  geom_point(
    size = 3
  ) +

  theme_classic() +

  labs(

    title =
      "Population-level transcriptomic PCA",

    x =
      paste0(
        "PC1 (",
        round(
          100 * pca_variance[1],
          1
        ),
        "%)"
      ),

    y =
      paste0(
        "PC2 (",
        round(
          100 * pca_variance[2],
          1
        ),
        "%)"
      )

  )


ggsave(
  file.path(
    output_dir,
    "transcriptomic_PCA.pdf"
  ),
  pca_plot,
  width = 7,
  height = 6
)


# ============================================================
# 13. RESPONSE MATRIX FOR RDA
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Preparing transcriptomic response matrix for RDA\n"
)

cat(
  "============================================\n\n"
)


if (scale_response_for_rda) {

  rda_response <-
    scale(
      t(
        population_expression
      ),
      center = TRUE,
      scale = TRUE
    )

  rda_response_description <-
    "Population expression profiles standardized gene-wise"

} else {

  rda_response <-
    t(
      population_expression
    )

  rda_response_description <-
    "Population expression profiles not gene-standardized"

}


rda_response <-
  as.matrix(
    rda_response
  )


rownames(
  rda_response
) <-
  population_ids


cat(
  rda_response_description,
  "\n"
)


write.csv(
  rda_response,
  file.path(
    output_dir,
    "RDA_transcriptomic_response_matrix.csv"
  )
)


# ============================================================
# 14. TRANSCRIPTOMIC DISTANCE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Calculating transcriptomic distances\n"
)

cat(
  "============================================\n\n"
)


transcriptomic_distance <-
  dist(
    rda_response,
    method =
      transcriptomic_distance_method
  )


attr(
  transcriptomic_distance,
  "Labels"
) <-
  population_ids


validate_distance(
  transcriptomic_distance,
  population_ids
)


write.csv(
  as.matrix(
    transcriptomic_distance
  ),
  file.path(
    output_dir,
    "transcriptomic_distance_matrix.csv"
  )
)


# ============================================================
# 15. GEOGRAPHIC DISTANCE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Calculating geographic distances\n"
)

cat(
  "============================================\n\n"
)


population_coordinates <-
  metadata_expression %>%

  group_by(
    Population
  ) %>%

  summarise(

    Latitude =
      mean(
        Latitude,
        na.rm = TRUE
      ),

    Longitude =
      mean(
        Longitude,
        na.rm = TRUE
      ),

    .groups = "drop"

  )


population_coordinates <-
  population_coordinates[
    match(
      population_ids,
      population_coordinates$Population
    ),
    ,
    drop = FALSE
  ]


if (anyNA(
  population_coordinates$Latitude
) ||
anyNA(
  population_coordinates$Longitude
)) {

  stop(
    "\nMissing geographic coordinates."
  )

}


geographic_distance <-
  calculate_geographic_distance(
    population_coordinates
  )


validate_distance(
  geographic_distance,
  population_ids
)


write.csv(
  as.matrix(
    geographic_distance
  ),
  file.path(
    output_dir,
    "geographic_distance_matrix_km.csv"
  )
)


# ============================================================
# 16. POPULATION-LEVEL ENVIRONMENT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Creating population-level environmental data\n"
)

cat(
  "============================================\n\n"
)


population_environment <-
  metadata_expression %>%

  group_by(
    Population
  ) %>%

  summarise(

    across(
      all_of(
        climate_variables
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),

    ECe_mean =
      mean(
        ECe_mean,
        na.rm = TRUE
      ),

    .groups = "drop"

  )


population_environment <-
  population_environment[
    match(
      population_ids,
      population_environment$Population
    ),
    ,
    drop = FALSE
  ]


write_csv(
  population_environment,
  file.path(
    output_dir,
    "population_environment.csv"
  )
)


# ============================================================
# 17. STANDARDIZE ENVIRONMENT
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Standardizing environmental predictors\n"
)

cat(
  "============================================\n\n"
)


# Convert to ordinary data.frame first.
#
# This avoids:
#
#   "Setting row names on a tibble is deprecated."
#

environment_data <-
  population_environment %>%

  mutate(

    across(
      all_of(
        climate_variables
      ),
      z_standardize,
      .names = "{.col}_z"
    ),

    Salinity_z =
      z_standardize(
        ECe_mean
      )

  )


environment_data <-
  as.data.frame(
    environment_data
  )


rownames(
  environment_data
) <-
  environment_data$Population


environment_predictor_names <-
  c(

    paste0(
      climate_variables,
      "_z"
    ),

    "Salinity_z"

  )


environment_matrix <-
  environment_data[
    ,
    environment_predictor_names,
    drop = FALSE
  ]


write.csv(
  environment_matrix,
  file.path(
    output_dir,
    "environment_predictors_standardized.csv"
  )
)


# ============================================================
# 18. ENVIRONMENTAL DISTANCES
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Calculating environmental distances\n"
)

cat(
  "============================================\n\n"
)


climate_matrix <-
  environment_data[
    ,
    paste0(
      climate_variables,
      "_z"
    ),
    drop = FALSE
  ]


climate_distance <-
  dist(
    climate_matrix,
    method = "euclidean"
  )


attr(
  climate_distance,
  "Labels"
) <-
  population_ids


salinity_distance <-
  dist(
    environment_data$Salinity_z
  )


attr(
  salinity_distance,
  "Labels"
) <-
  population_ids


environment_distance <-
  dist(
    environment_matrix,
    method = "euclidean"
  )


attr(
  environment_distance,
  "Labels"
) <-
  population_ids


write.csv(
  as.matrix(
    climate_distance
  ),
  file.path(
    output_dir,
    "climate_distance_matrix.csv"
  )
)


write.csv(
  as.matrix(
    salinity_distance
  ),
  file.path(
    output_dir,
    "salinity_distance_matrix.csv"
  )
)


write.csv(
  as.matrix(
    environment_distance
  ),
  file.path(
    output_dir,
    "combined_environment_distance_matrix.csv"
  )
)


# ============================================================
# 19. ENVIRONMENTAL COLLINEARITY
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Checking environmental collinearity\n"
)

cat(
  "============================================\n\n"
)


environment_correlations <-
  cor(
    environment_matrix
  )


write.csv(
  environment_correlations,
  file.path(
    output_dir,
    "environment_predictor_correlations.csv"
  )
)


# ============================================================
# 20. ENVIRONMENT-ONLY RDA
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Environment-only RDA\n"
)

cat(
  "============================================\n\n"
)


environment_rda <-
  vegan::rda(

    rda_response ~
      Climate_PC1_z +
      Climate_PC2_z +
      Climate_PC3_z +
      Salinity_z,

    data =
      environment_data

  )


environment_rda_anova <-
  anova(
    environment_rda,
    permutations =
      n_permutations
  )


environment_rda_terms <-
  anova(
    environment_rda,
    by = "margin",
    permutations =
      n_permutations
  )


environment_rda_r2 <-
  safe_adj_r2(
    environment_rda
  )


environment_vif <-
  vegan::vif.cca(
    environment_rda
  )


write.csv(
  as.data.frame(
    environment_rda_anova
  ),
  file.path(
    output_dir,
    "RDA_environment_overall_test.csv"
  )
)


write.csv(
  as.data.frame(
    environment_rda_terms
  ),
  file.path(
    output_dir,
    "RDA_environment_marginal_tests.csv"
  )
)


write.csv(
  data.frame(

    Predictor =
      names(
        environment_vif
      ),

    VIF =
      as.numeric(
        environment_vif
      )

  ),
  file.path(
    output_dir,
    "RDA_environment_VIF.csv"
  ),
  row.names = FALSE
)


write_csv(
  environment_rda_r2,
  file.path(
    output_dir,
    "RDA_environment_R2.csv"
  )
)


capture.output(
  summary(
    environment_rda
  ),
  file =
    file.path(
      output_dir,
      "RDA_environment_summary.txt"
    )
)


# ============================================================
# 21. dbMEMs
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Generating spatial dbMEMs\n"
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# Check whether the geographic distance is Euclidean.
#
# dbmem() accepts a dist object, but internally performs PCoA.
# A non-Euclidean geographic distance can therefore generate
# a warning from dudi.pco().
#
# The great-circle distance is retained for geographic
# distance analyses.
#
# For dbMEM construction only, a Cailliez-transformed Euclidean
# version is used if necessary.
# ------------------------------------------------------------


geographic_distance_is_euclidean <-
  ade4::is.euclid(
    geographic_distance
  )


cat(
  "Is raw geographic distance Euclidean?:",
  geographic_distance_is_euclidean,
  "\n"
)


if (geographic_distance_is_euclidean) {

  geographic_distance_dbmem <-
    geographic_distance

  dbmem_distance_transformation <-
    "None; raw great-circle distance used"

} else {

  cat(
    "\nRaw geographic distance is non-Euclidean.\n"
  )

  cat(
    "Applying Cailliez transformation for dbMEM construction only.\n\n"
  )

  geographic_distance_dbmem <-
    ade4::cailliez(
      geographic_distance
    )

  dbmem_distance_transformation <-
    "Cailliez transformation applied for dbMEM only"

}


attr(
  geographic_distance_dbmem,
  "Labels"
) <-
  population_ids


write.csv(
  as.matrix(
    geographic_distance_dbmem
  ),
  file.path(
    output_dir,
    "geographic_distance_matrix_for_dbMEM.csv"
  )
)


writeLines(
  dbmem_distance_transformation,
  file.path(
    output_dir,
    "dbMEM_distance_transformation.txt"
  )
)


# ------------------------------------------------------------
# Generate dbMEMs
# ------------------------------------------------------------

dbmem_object <-
  adespatial::dbmem(

    geographic_distance_dbmem,

    thresh =
      NULL,

    MEM.autocor =
      spatial_mem_autocorrelation,

    store.listw =
      TRUE,

    silent =
      FALSE

  )


dbmem_matrix <-
  as.data.frame(
    dbmem_object
  )


# IMPORTANT:
#
# Keep the native MEM names from adespatial.
#
# Do NOT rename them manually, because mem.select()
# returns MEM names that must match these columns.


dbmem_matrix <-
  as.matrix(
    dbmem_matrix
  )


rownames(
  dbmem_matrix
) <-
  population_ids


if (ncol(
  dbmem_matrix
) == 0) {

  stop(
    "\nNo positive dbMEMs were generated."
  )

}


write.csv(
  dbmem_matrix,
  file.path(
    output_dir,
    "all_positive_dbMEMs.csv"
  )
)


dbmem_eigenvalues <-
  attr(
    dbmem_object,
    "values"
  )


write.csv(
  data.frame(

    dbMEM =
      colnames(
        dbmem_matrix
      ),

    Eigenvalue =
      dbmem_eigenvalues[
        seq_len(
          ncol(
            dbmem_matrix
          )
        )
      ]

  ),
  file.path(
    output_dir,
    "dbMEM_eigenvalues.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. OVERALL SPATIAL MEM SELECTION
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Selecting spatial dbMEMs from transcriptomic response\n"
)

cat(
  "============================================\n\n"
)


mem_selection_overall <-
  adespatial::mem.select(

    x =
      rda_response,

    listw =
      attr(
        dbmem_object,
        "listw"
      ),

    MEM.autocor =
      spatial_mem_autocorrelation,

    method =
      spatial_selection_method,

    MEM.all =
      FALSE,

    nperm =
      spatial_selection_nperm,

    nperm.global =
      spatial_global_nperm,

    alpha =
      spatial_selection_alpha,

    verbose =
      TRUE

  )


overall_spatial_p <-
  extract_mem_p(
    mem_selection_overall$global.test
  )


overall_spatial_detected <-
  is.finite(
    overall_spatial_p
  ) &&
  overall_spatial_p <
  spatial_selection_alpha


capture.output(

  print(
    mem_selection_overall$global.test
  ),

  file =
    file.path(
      output_dir,
      "dbMEM_overall_global_test.txt"
    )

)


if (
  overall_spatial_detected &&
  !is.null(
    mem_selection_overall$MEM.select
  )
) {

  overall_selected_mem_matrix <-
    as.matrix(
      mem_selection_overall$MEM.select
    )

  rownames(
    overall_selected_mem_matrix
  ) <-
    population_ids

} else {

  overall_selected_mem_matrix <-
    matrix(
      nrow =
        nrow(rda_response),
      ncol =
        0
    )

  rownames(
    overall_selected_mem_matrix
  ) <-
    population_ids

}


overall_selected_mem_matrix <-
  cap_mem_matrix(
    overall_selected_mem_matrix,
    max_spatial_predictors
  )


overall_selected_mem_names <-
  colnames(
    overall_selected_mem_matrix
  )


if (length(
  overall_selected_mem_names
) == 0) {

  cat(
    "\nNo significant overall positive spatial structure detected.\n"
  )

} else {

  cat(
    "\nOverall spatial MEMs selected:",
    length(
      overall_selected_mem_names
    ),
    "\n"
  )

  writeLines(
    overall_selected_mem_names,
    file.path(
      output_dir,
      "overall_selected_dbMEM_names.txt"
    )
  )

}


write.csv(
  overall_selected_mem_matrix,
  file.path(
    output_dir,
    "overall_selected_dbMEMs.csv"
  )
)


if (!is.null(
  mem_selection_overall$summary
)) {

  write.csv(
    as.data.frame(
      mem_selection_overall$summary
    ),
    file.path(
      output_dir,
      "overall_dbMEM_selection_summary.csv"
    ),
    row.names = FALSE
  )

}


# ============================================================
# 23. SPATIAL-ONLY RDA:
#     OVERALL GEOGRAPHIC STRUCTURE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Spatial-only RDA: overall geographic structure\n"
)

cat(
  "============================================\n\n"
)


if (length(
  overall_selected_mem_names
) > 0) {

  overall_spatial_data <-
    as.data.frame(
      overall_selected_mem_matrix
    )


  spatial_rda_overall <-
    vegan::rda(

      rda_response ~ .,

      data =
        overall_spatial_data

    )


  spatial_rda_overall_anova <-
    anova(
      spatial_rda_overall,
      permutations =
        n_permutations
    )


  spatial_rda_overall_terms <-
    anova(
      spatial_rda_overall,
      by = "margin",
      permutations =
        n_permutations
    )


  spatial_rda_overall_r2 <-
    safe_adj_r2(
      spatial_rda_overall
    )

  # --------------------------------------------------------
  # VIF for the spatial-only model
  #
  # Only meaningful/computable with 2+ predictors, which is
  # guaranteed here since this block only runs when
  # overall_selected_mem_names has length > 0. With exactly
  # one selected MEM, vif.cca() would fail or be trivial, so
  # we guard on ncol > 1.
  # --------------------------------------------------------

  if (ncol(overall_spatial_data) > 1) {

    spatial_vif <-
      vegan::vif.cca(
        spatial_rda_overall
      )

    write.csv(
      data.frame(

        Predictor =
          names(
            spatial_vif
          ),

        VIF =
          as.numeric(
            spatial_vif
          )

      ),
      file.path(
        output_dir,
        "RDA_spatial_overall_VIF.csv"
      ),
      row.names = FALSE
    )

  } else {

    cat(
      "\nOnly one spatial predictor selected; VIF not computed (not meaningful with a single predictor).\n"
    )

  }

  write.csv(
    as.data.frame(
      spatial_rda_overall_anova
    ),
    file.path(
      output_dir,
      "RDA_spatial_overall_structure_test.csv"
    )
  )


  write.csv(
    as.data.frame(
      spatial_rda_overall_terms
    ),
    file.path(
      output_dir,
      "RDA_spatial_overall_structure_marginal.csv"
    )
  )


  write_csv(
    spatial_rda_overall_r2,
    file.path(
      output_dir,
      "RDA_spatial_overall_structure_R2.csv"
    )
  )


  capture.output(
    summary(
      spatial_rda_overall
    ),
    file =
      file.path(
        output_dir,
        "RDA_spatial_overall_structure_summary.txt"
      )
  )

} else {

  spatial_rda_overall <-
    NULL

  spatial_rda_overall_anova <-
    NULL

  spatial_rda_overall_terms <-
    NULL

  spatial_rda_overall_r2 <-
    tibble(
      R2 = NA_real_,
      adjusted_R2 = NA_real_
    )

}


# ============================================================
# 24. ENVIRONMENT-CONDITIONED SPATIAL MEM SELECTION
# This answer: is there additional undiscovered spatial signal?
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Selecting spatial dbMEMs after accounting for environment\n"
)

cat(
  "============================================\n\n"
)


environment_residuals <-
  residuals(
    environment_rda
  )


environment_residuals <-
  as.matrix(
    environment_residuals
  )


rownames(
  environment_residuals
) <-
  population_ids


mem_selection_conditional <-
  adespatial::mem.select(

    x =
      environment_residuals,

    listw =
      attr(
        dbmem_object,
        "listw"
      ),

    MEM.autocor =
      spatial_mem_autocorrelation,

    method =
      spatial_selection_method,

    MEM.all =
      FALSE,

    nperm =
      spatial_selection_nperm,

    nperm.global =
      spatial_global_nperm,

    alpha =
      spatial_selection_alpha,

    verbose =
      TRUE

  )


conditional_spatial_p <-
  extract_mem_p(
    mem_selection_conditional$global.test
  )


conditional_spatial_detected <-
  is.finite(
    conditional_spatial_p
  ) &&
  conditional_spatial_p <
  spatial_selection_alpha


capture.output(

  print(
    mem_selection_conditional$global.test
  ),

  file =
    file.path(
      output_dir,
      "dbMEM_conditional_global_test.txt"
    )

)


if (
  conditional_spatial_detected &&
  !is.null(
    mem_selection_conditional$MEM.select
  )
) {

  conditional_selected_mem_matrix <-
    as.matrix(
      mem_selection_conditional$MEM.select
    )

  rownames(
    conditional_selected_mem_matrix
  ) <-
    population_ids

} else {

  conditional_selected_mem_matrix <-
    matrix(
      nrow =
        nrow(rda_response),
      ncol =
        0
    )

  rownames(
    conditional_selected_mem_matrix
  ) <-
    population_ids

}


conditional_selected_mem_matrix <-
  cap_mem_matrix(
    conditional_selected_mem_matrix,
    max_spatial_predictors
  )


conditional_selected_mem_names <-
  colnames(
    conditional_selected_mem_matrix
  )


if (length(
  conditional_selected_mem_names
) == 0) {

  cat(
    "\nNo significant positive spatial structure remains after accounting for environment.\n"
  )

} else {

  cat(
    "\nEnvironment-conditioned spatial MEMs selected:",
    length(
      conditional_selected_mem_names
    ),
    "\n"
  )

  writeLines(
    conditional_selected_mem_names,
    file.path(
      output_dir,
      "conditional_selected_dbMEM_names.txt"
    )
  )

}


write.csv(
  conditional_selected_mem_matrix,
  file.path(
    output_dir,
    "conditional_selected_dbMEMs.csv"
  )
)


if (!is.null(
  mem_selection_conditional$summary
)) {

  write.csv(
    as.data.frame(
      mem_selection_conditional$summary
    ),
    file.path(
      output_dir,
      "conditional_dbMEM_selection_summary.csv"
    ),
    row.names = FALSE
  )

}


# ============================================================
# 25. FULL ENVIRONMENT + SPATIAL RDA
#
#     Uses overall (unconditioned) spatial MEMs
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Full environment + spatial RDA\n"
)

cat(
  "============================================\n\n"
)

if (length(
  overall_selected_mem_names
) > 0) {

  model_data <-
    cbind(

      environment_data,

      as.data.frame(
        overall_selected_mem_matrix
      )

    )

  full_predictor_names <-
    c(
      environment_predictor_names,
      overall_selected_mem_names
    )

} else {

  model_data <-
    environment_data

  full_predictor_names <-
    environment_predictor_names

}

full_formula <-
  as.formula(
    paste(
      "rda_response ~",
      paste(
        full_predictor_names,
        collapse = " + "
      )
    )
  )


full_rda <-
  vegan::rda(
    full_formula,
    data =
      model_data
  )


full_rda_anova <-
  anova(
    full_rda,
    permutations =
      n_permutations
  )


full_rda_terms <-
  anova(
    full_rda,
    by = "margin",
    permutations =
      n_permutations
  )


full_rda_r2 <-
  safe_adj_r2(
    full_rda
  )


full_vif <-
  vegan::vif.cca(
    full_rda
  )


write.csv(
  as.data.frame(
    full_rda_anova
  ),
  file.path(
    output_dir,
    "RDA_full_overall_test.csv"
  )
)


write.csv(
  as.data.frame(
    full_rda_terms
  ),
  file.path(
    output_dir,
    "RDA_full_marginal_tests.csv"
  )
)


write.csv(
  data.frame(

    Predictor =
      names(
        full_vif
      ),

    VIF =
      as.numeric(
        full_vif
      )

  ),
  file.path(
    output_dir,
    "RDA_full_VIF.csv"
  ),
  row.names = FALSE
)


write_csv(
  full_rda_r2,
  file.path(
    output_dir,
    "RDA_full_R2.csv"
  )
)


capture.output(
  summary(
    full_rda
  ),
  file =
    file.path(
      output_dir,
      "RDA_full_summary.txt"
    )
)


# ============================================================
# 26. PARTIAL RDA:
#     ENVIRONMENT | SPATIAL STRUCTURE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Partial RDA: environment conditioned on spatial structure\n"
)

cat(
  "============================================\n\n"
)

if (length(
  overall_selected_mem_names
) > 0) {

  environment_given_space_formula <-
    as.formula(

      paste(

        "rda_response ~",

        paste(
          environment_predictor_names,
          collapse = " + "
        ),

        "+ Condition(",

        paste(
          overall_selected_mem_names,
          collapse = " + "
        ),

        ")"

      )

    )

} else {

  environment_given_space_formula <-
    as.formula(

      paste(

        "rda_response ~",

        paste(
          environment_predictor_names,
          collapse = " + "
        )

      )

    )

}

environment_given_space_rda <-
  vegan::rda(
    environment_given_space_formula,
    data =
      model_data
  )


environment_given_space_anova <-
  anova(
    environment_given_space_rda,
    permutations =
      n_permutations
  )


environment_given_space_terms <-
  anova(
    environment_given_space_rda,
    by = "margin",
    permutations =
      n_permutations
  )


environment_given_space_r2 <-
  safe_adj_r2(
    environment_given_space_rda
  )


write.csv(
  as.data.frame(
    environment_given_space_anova
  ),
  file.path(
    output_dir,
    "partial_RDA_environment_given_space_overall.csv"
  )
)


write.csv(
  as.data.frame(
    environment_given_space_terms
  ),
  file.path(
    output_dir,
    "partial_RDA_environment_given_space_marginal.csv"
  )
)


write_csv(
  environment_given_space_r2,
  file.path(
    output_dir,
    "partial_RDA_environment_given_space_R2.csv"
  )
)


capture.output(
  summary(
    environment_given_space_rda
  ),
  file =
    file.path(
      output_dir,
      "partial_RDA_environment_given_space_summary.txt"
    )
)


# ============================================================
# 27. PARTIAL RDA: SPACE | ENVIRONMENT
# This answers: do the known spatial axes matter independent of environment?
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Partial RDA: spatial structure conditioned on environment\n"
)

cat(
  "============================================\n\n"
)

if (length(
  overall_selected_mem_names
) > 0) {

  spatial_given_environment_formula <-
    as.formula(

      paste(

        "rda_response ~",

        paste(
          overall_selected_mem_names,
          collapse = " + "
        ),

        "+ Condition(",

        paste(
          environment_predictor_names,
          collapse = " + "
        ),

        ")"

      )

    )

  spatial_given_environment_rda <-
    vegan::rda(
      spatial_given_environment_formula,
      data =
        model_data
    )

    spatial_given_environment_anova <-
    anova(
      spatial_given_environment_rda,
      permutations =
        n_permutations
    )


  spatial_given_environment_terms <-
    anova(
      spatial_given_environment_rda,
      by = "margin",
      permutations =
        n_permutations
    )


  spatial_given_environment_r2 <-
    safe_adj_r2(
      spatial_given_environment_rda
    )


  write.csv(
    as.data.frame(
      spatial_given_environment_anova
    ),
    file.path(
      output_dir,
      "partial_RDA_space_given_environment_overall.csv"
    )
  )


  write.csv(
    as.data.frame(
      spatial_given_environment_terms
    ),
    file.path(
      output_dir,
      "partial_RDA_space_given_environment_marginal.csv"
    )
  )


  write_csv(
    spatial_given_environment_r2,
    file.path(
      output_dir,
      "partial_RDA_space_given_environment_R2.csv"
    )
  )


  capture.output(
    summary(
      spatial_given_environment_rda
    ),
    file =
      file.path(
        output_dir,
        "partial_RDA_space_given_environment_summary.txt"
      )
  )

} else {

  spatial_given_environment_rda <- NULL
  spatial_given_environment_anova <- NULL
  spatial_given_environment_terms <- NULL
  spatial_given_environment_r2 <-
    tibble(R2 = NA_real_, adjusted_R2 = NA_real_)

}

# ============================================================
# 28. VARIATION PARTITIONING
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Variation partitioning\n"
)

cat(
  "============================================\n\n"
)

if (length(
  overall_selected_mem_names
) > 0) {

  variation_partitioning <-
    vegan::varpart(

      rda_response,

      environment_matrix,

      overall_selected_mem_matrix,

      scale = FALSE

    )
	
	  capture.output(
    print(
      variation_partitioning
    ),
    file =
      file.path(
        output_dir,
        "variation_partitioning.txt"
      )
  )


  capture.output(
    summary(
      variation_partitioning
    ),
    file =
      file.path(
        output_dir,
        "variation_partitioning_summary.txt"
      )
  )


  write.csv(
    variation_partitioning$part$indfract,
    file.path(
      output_dir,
      "variation_partitioning_individual_fractions.csv"
    ),
    row.names = TRUE
  )


  write.csv(
    variation_partitioning$part$fract,
    file.path(
      output_dir,
      "variation_partitioning_fractions.csv"
    ),
    row.names = TRUE
  )

} else {

  variation_partitioning <- NULL

  cat(
    "Variation partitioning skipped because no overall positive spatial structure was selected.\n"
  )

}

# ============================================================
# 29. PAIRWISE DATA
# ============================================================

ibd_plot_data <-
  distance_pair_table(

    transcriptomic_distance,

    geographic_distance,

    "Transcriptomic_distance",

    "Geographic_distance_km"

  )


ibe_climate_plot_data <-
  distance_pair_table(

    transcriptomic_distance,

    climate_distance,

    "Transcriptomic_distance",

    "Climate_distance"

  )


ibe_environment_plot_data <-
  distance_pair_table(

    transcriptomic_distance,

    environment_distance,

    "Transcriptomic_distance",

    "Environmental_distance"

  )


write_csv(
  ibd_plot_data,
  file.path(
    output_dir,
    "IBD_pairwise_distances.csv"
  )
)


write_csv(
  ibe_climate_plot_data,
  file.path(
    output_dir,
    "IBE_climate_pairwise_distances.csv"
  )
)


write_csv(
  ibe_environment_plot_data,
  file.path(
    output_dir,
    "IBE_environment_pairwise_distances.csv"
  )
)


# ============================================================
# 30. DISTANCE SUMMARY
# ============================================================

distance_summary <-
  tibble(

    distance_type =
      c(
        "Transcriptomic",
        "Geographic",
        "Climate",
        "Salinity",
        "Climate + salinity"
      ),

    mean_distance =
      c(

        mean(
          transcriptomic_distance
        ),

        mean(
          geographic_distance
        ),

        mean(
          climate_distance
        ),

        mean(
          salinity_distance
        ),

        mean(
          environment_distance
        )

      ),

    median_distance =
      c(

        median(
          transcriptomic_distance
        ),

        median(
          geographic_distance
        ),

        median(
          climate_distance
        ),

        median(
          salinity_distance
        ),

        median(
          environment_distance
        )

      )

  )


write_csv(
  distance_summary,
  file.path(
    output_dir,
    "distance_matrix_summary.csv"
  )
)


# ============================================================
# 31. RESULTS
# ============================================================

primary_results <-
  tibble(

    analysis =
      c(

        "Environment only",

        "Overall spatial structure",

        "Environment | spatial structure",

        "Spatial structure | environment",

        "Full environment + overall spatial model"

      ),

    R2 =
      c(

        environment_rda_r2$R2,

        spatial_rda_overall_r2$R2,

        environment_given_space_r2$R2,

        if (
          !is.null(
            spatial_given_environment_rda
          )
        )
          spatial_given_environment_r2$R2
        else
          NA_real_,

        full_rda_r2$R2

      ),

    adjusted_R2 =
      c(

        environment_rda_r2$adjusted_R2,

        spatial_rda_overall_r2$adjusted_R2,

        environment_given_space_r2$adjusted_R2,

        if (
          !is.null(
            spatial_given_environment_rda
          )
        )
          spatial_given_environment_r2$adjusted_R2
        else
          NA_real_,

        full_rda_r2$adjusted_R2

      ),

    P_value =
      c(

        extract_anova_p(
          environment_rda_anova
        ),

        extract_anova_p(
          spatial_rda_overall_anova
        ),

        extract_anova_p(
          environment_given_space_anova
        ),

        extract_anova_p(
          spatial_given_environment_anova
        ),

        extract_anova_p(
          full_rda_anova
        )

      )

  )


write_csv(
  primary_results,
  file.path(
    output_dir,
    "primary_RDA_results.csv"
  )
)


# ============================================================
# 32. SPATIAL INTERPRETATION SUMMARY
# ============================================================

spatial_interpretation <-
  tibble(

    overall_spatial_structure_detected =
      overall_spatial_detected,

    overall_spatial_global_P =
      overall_spatial_p,

    overall_MEMs_selected =
      length(
        overall_selected_mem_names
      ),

    environment_conditioned_spatial_structure_detected =
      conditional_spatial_detected,

    environment_conditioned_spatial_global_P =
      conditional_spatial_p,

    environment_conditioned_MEMs_selected =
      length(
        conditional_selected_mem_names
      ),

    interpretation =
      if (
        overall_spatial_detected &&
        !conditional_spatial_detected
      ) {

        "Overall spatial structure detected, but no significant positive spatial structure remained after accounting for measured environment."

      } else if (
        overall_spatial_detected &&
        conditional_spatial_detected
      ) {

        "Spatial structure detected both overall and after accounting for measured environment."

      } else if (
        !overall_spatial_detected
      ) {

        "No significant overall positive spatial structure detected."

      } else {

        "Spatial result requires inspection of the saved model outputs."

      }

  )


write_csv(
  spatial_interpretation,
  file.path(
    output_dir,
    "spatial_interpretation_summary.csv"
  )
)


# ============================================================
# 36. ANALYSIS SUMMARY
# ============================================================

analysis_summary <-
  tibble(

    n_samples =
      nrow(
        metadata_expression
      ),

    n_populations =
      length(
        population_ids
      ),

    n_genes =
      nrow(
        population_expression
      ),

    n_dbMEMs_generated =
      ncol(
        dbmem_matrix
      ),

    n_overall_MEMs_selected =
      length(
        overall_selected_mem_names
      ),

    n_conditional_MEMs_selected =
      length(
        conditional_selected_mem_names
      ),

    overall_spatial_global_P =
      overall_spatial_p,

    conditional_spatial_global_P =
      conditional_spatial_p,

    environment_overall_P =
      extract_anova_p(
        environment_rda_anova
      ),

    environment_overall_adjusted_R2 =
      environment_rda_r2$adjusted_R2,

    overall_spatial_adjusted_R2 =
      spatial_rda_overall_r2$adjusted_R2,

    environment_given_space_P =
      extract_anova_p(
        environment_given_space_anova
      ),

    environment_given_space_adjusted_R2 =
      environment_given_space_r2$adjusted_R2,

    spatial_given_environment_P =
      extract_anova_p(
        spatial_given_environment_anova
      ),

    spatial_given_environment_adjusted_R2 =
      if (
        !is.null(
          spatial_given_environment_rda
        )
      )
        spatial_given_environment_r2$adjusted_R2
      else
        NA_real_,

    full_model_P =
      extract_anova_p(
        full_rda_anova
      ),

    full_model_adjusted_R2 =
      full_rda_r2$adjusted_R2

  )


write_csv(
  analysis_summary,
  file.path(
    output_dir,
    "analysis_summary.csv"
  )
)


# ============================================================
# 34. SAVE COMPLETE ANALYSIS OBJECTS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Saving analysis objects\n"
)

cat(
  "============================================\n\n"
)


saveRDS(

  list(

    population_expression =
      population_expression,

    rda_response =
      rda_response,

    transcriptomic_distance =
      transcriptomic_distance,

    geographic_distance =
      geographic_distance,

    geographic_distance_dbmem =
      geographic_distance_dbmem,

    climate_distance =
      climate_distance,

    salinity_distance =
      salinity_distance,

    environment_distance =
      environment_distance,

    environment_data =
      environment_data,

    dbmem_object =
      dbmem_object,

    dbmem_matrix =
      dbmem_matrix,

    mem_selection_overall =
      mem_selection_overall,

    overall_selected_mem_matrix =
      overall_selected_mem_matrix,

    mem_selection_conditional =
      mem_selection_conditional,

    conditional_selected_mem_matrix =
      conditional_selected_mem_matrix,

    environment_rda =
      environment_rda,

    spatial_rda_overall =
      spatial_rda_overall,

    full_rda =
      full_rda,

    environment_given_space_rda =
      environment_given_space_rda,

    spatial_given_environment_rda =
      spatial_given_environment_rda,

    variation_partitioning =
      variation_partitioning
  ),

  file.path(
    output_dir,
    "IBD_IBE_complete_analysis_objects.rds"
  )

)


# ============================================================
# 38. FINAL CONSOLE SUMMARY
# ============================================================

cat(
  "\n\n"
)

cat(
  "============================================================\n"
)

cat(
  "POPULATION-LEVEL TRANSCRIPTOMIC GEOGRAPHY / ENVIRONMENT\n"
)

cat(
  "ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Number of populations:",
  length(
    population_ids
  ),
  "\n"
)


cat(
  "Number of genes:",
  nrow(
    population_expression
  ),
  "\n\n"
)


cat(
  "------------------------------------------------------------\n"
)

cat(
  "ENVIRONMENT\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "Environment-only adjusted R2:",
  round(
    environment_rda_r2$adjusted_R2,
    4
  ),
  "\n"
)


cat(
  "Environment-only P:",
  format.pval(
    extract_anova_p(
      environment_rda_anova
    ),
    digits = 4
  ),
  "\n\n"
)


cat(
  "------------------------------------------------------------\n"
)

cat(
  "OVERALL SPATIAL STRUCTURE\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "Overall spatial global P:",
  format.pval(
    overall_spatial_p,
    digits = 4
  ),
  "\n"
)


cat(
  "Overall spatial MEMs selected:",
  length(
    overall_selected_mem_names
  ),
  "\n"
)


cat(
  "Overall spatial RDA adjusted R2:",
  round(
    spatial_rda_overall_r2$adjusted_R2,
    4
  ),
  "\n\n"
)


cat(
  "------------------------------------------------------------\n"
)

cat(
  "SPATIAL STRUCTURE AFTER ENVIRONMENT\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "Conditional spatial global P:",
  format.pval(
    conditional_spatial_p,
    digits = 4
  ),
  "\n"
)


cat(
  "Conditional spatial MEMs selected:",
  length(
    conditional_selected_mem_names
  ),
  "\n\n"
)


if (
  !conditional_spatial_detected
) {

  cat(
    "Interpretation:\n"
  )

  cat(
    "No significant positive spatial structure remained after\n"
  )

  cat(
    "accounting for climate + salinity.\n\n"
  )

}


cat(
  "------------------------------------------------------------\n"
)

cat(
  "PARTIAL RDA\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "Environment | space adjusted R2:",
  round(
    environment_given_space_r2$adjusted_R2,
    4
  ),
  "\n"
)


cat(
  "Environment | space P:",
  format.pval(
    extract_anova_p(
      environment_given_space_anova
    ),
    digits = 4
  ),
  "\n\n"
)


if (!is.null(
  spatial_given_environment_rda
)) {

  cat(
    "Space | environment adjusted R2:",
    round(
      spatial_given_environment_r2$adjusted_R2,
      4
    ),
    "\n"
  )

  cat(
    "Space | environment P:",
    format.pval(
      extract_anova_p(
        spatial_given_environment_anova
      ),
      digits = 4
    ),
    "\n\n"
  )

} else {

  cat(
    "Space | environment: not fitted because no significant\n"
  )

  cat(
    "positive residual spatial structure was detected.\n\n"
  )

}


cat(
  "------------------------------------------------------------\n"
)

cat(
  "FULL MODEL\n"
)

cat(
  "------------------------------------------------------------\n"
)


cat(
  "Full model adjusted R2:",
  round(
    full_rda_r2$adjusted_R2,
    4
  ),
  "\n"
)


cat(
  "Full model P:",
  format.pval(
    extract_anova_p(
      full_rda_anova
    ),
    digits = 4
  ),
  "\n\n"
)

cat(
  "Results saved to:\n",
  normalizePath(
    output_dir,
    winslash = "/"
  ),
  "\n\n"
)


cat(
  "Software versions:\n",
  file.path(
    output_dir,
    "software_versions.csv"
  ),
  "\n"
)


cat(
  "Session information:\n",
  file.path(
    output_dir,
    "sessionInfo.txt"
  ),
  "\n\n"
)


cat(
  "============================================================\n"
)