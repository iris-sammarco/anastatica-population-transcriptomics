#!/usr/bin/env Rscript

# ==============================================================================
# Gene-set analysis with zenith for the FULL salinity + climate + distribution
# range dream model using multiple functional annotation collections
#
# Species  : Anastatica hierochuntica
# Author   : Iris Sammarco
# Date     : 09/09/2026
#
# ==============================================================================
# PURPOSE
# ==============================================================================
#
# This script performs gene-set analyses for the common-garden transcriptome
# using the pre-fitted FULL dream/eBayes model:
#
#   expression ~ ECe_mean
#              + Climate_PC1 + Climate_PC2 + Climate_PC3
#              + Site_of_distribution_range
#              + (1 | Sequencing_batch)
#              + (1 | Population)
#
# The saved dream/eBayes fit is used directly; the model is NOT re-fitted.
#
# The functional analyses are performed with zenith_gsa(), which implements
# competitive gene-set testing for dream results while accounting for
# correlation among genes within gene sets.
#
#
# ==============================================================================
# PRIMARY AND SECONDARY FUNCTIONAL ANALYSES
# ==============================================================================
#
# The following annotation collections are analyzed separately:
#
#   1. GO Biological Process (GO BP)
#      PRIMARY functional analysis.
#
#   2. GO Molecular Function (GO MF)
#      PRIMARY functional analysis.
#
#   3. GO Cellular Component (GO CC)
#      PRIMARY functional analysis.
#
#   4. InterPro
#      SECONDARY functional/domain analysis.
#
#   5. PFAM
#      SECONDARY protein-domain analysis.
#
#   6. COG
#      SECONDARY broad functional-category analysis.
#
#   7. CAZy
#      EXPLORATORY / TARGETED carbohydrate-active-enzyme analysis.
#
# Each annotation collection is analyzed independently and receives its own
# set of zenith tests and within-test FDR correction.
#
# The annotation collections are deliberately NOT pooled into one statistical
# family because they represent different annotation systems with different
# biological meanings, coverage, redundancy, and numbers of gene sets.
#
#
# ==============================================================================
# BIOLOGICAL QUESTIONS
# ==============================================================================
#
# 1. Which functional gene sets are associated with each climatic axis after
#    accounting for salinity, distribution range, sequencing batch, and
#    population?
#
#    Tests:
#      Climate_PC1
#      Climate_PC2
#      Climate_PC3
#
# 2. Which functional gene sets differ among distribution-range classes after
#    accounting for salinity, climate PCs, sequencing batch, and population?
#
#    Tests:
#      East_vs_Center
#      West_vs_Center
#      East_vs_West
#
#
# ==============================================================================
# STATISTICAL PRINCIPLES
# ==============================================================================
#
# * The THREE-LEVEL East/Centre/West range factor is retained as the primary
#   biological definition.
#
# * At the gene level, the omnibus F-test for the three-level range factor is
#   the primary test of "any range effect", while pairwise contrasts describe
#   WHERE the range differences occur.
#
# * At the gene-set level, zenith_gsa() is run on individual estimable
#   coefficients/contrasts. zenith_gsa() does not take an arbitrary three-level
#   omnibus F-test as its normal input.
#
# * Therefore, the range omnibus F-test remains the primary gene-level
#   hypothesis, while the three pairwise zenith analyses provide directional
#   functional interpretation.
#
# * For climate, PC1, PC2 and PC3 are tested separately because they represent
#   distinct environmental axes. A separate climate omnibus F-test should be
#   considered at the gene level if available, but is NOT substituted for the
#   three directional zenith analyses.
#
# * inter.gene.cor = NA is used so that zenith estimates the average
#   within-gene-set correlation from the dream residuals rather than imposing
#   a fixed correlation across all sets.
#
# * zenith FDR values are controlled across gene sets WITHIN EACH
#   coefficient/contrast AND annotation collection.
#
# * A secondary global BH adjustment is calculated separately within each
#   annotation collection across all six tests. This is a sensitivity analysis.
#
# * The salinity term (ECe_mean) is retained in the fitted model as a covariate.
#   It is not removed simply because the corresponding gene-level DE analysis
#   yielded no significant genes.
#
#
# ==============================================================================
# IMPORTANT INTERPRETATION
# ==============================================================================
#
# A significant climate zenith result means that a functional gene set shows
# coordinated expression association with that climate PC CONDITIONAL ON the
# other terms in the full model.
#
# A significant range-contrast zenith result means that a functional gene set
# shows coordinated expression differentiation for that specific range
# contrast, CONDITIONAL ON salinity, climate PCs, sequencing batch, and
# population.
#
# Overlap between significant terms from different annotation collections is
# descriptive evidence of shared functional structure. It is NOT a formal
# interaction test and does not represent independent statistical evidence
# because annotation collections can overlap substantially.
#
#
# ==============================================================================
# ANNOTATION COLLECTIONS
# ==============================================================================
#
# The supplied annotation file contains:
#
#   GeneID
#   TranscriptID
#   Feature
#   Contig
#   Start
#   Stop
#   Strand
#   Name
#   Product
#   Alias/Synonyms
#   EC_number
#   BUSCO
#   PFAM
#   InterPro
#   EggNog
#   COG
#   GO Terms
#   Secreted
#   Membrane
#   Protease
#   CAZyme
#   Notes
#   gDNA
#   mRNA
#   CDS-transcript
#   Translation
#
# Parsing:
#
#   GO       -> GO:1234567 identifiers extracted from GO Terms
#   PFAM     -> PF12345 identifiers extracted from PFAM
#   InterPro -> IPR123456 identifiers extracted from InterPro
#   COG      -> one-letter COG functional categories extracted from COG
#   CAZy     -> comma-separated CAZyme identifiers
#
# GO terms are additionally filtered to Biological Process using GO.db.
#
#
# ==============================================================================
# REQUIRED INPUTS
# ==============================================================================
#
# 1. Saved full dream/eBayes fit:
#
#   /Anastatica/degs/ECe_climate_range_model/ECe_climate_range_model_fit_contrast_ebayes.rds
#
# 2. Anastatica annotation:
#
#   /Anastatica/genome_annotation/Anastatica_hierochuntica_ANAS.annotations.xlsx
#
#
# ==============================================================================
# OUTPUT STRUCTURE
# ==============================================================================
#
#   /Anastatica/zenith/ECe_climate_range_model/
#
#       GO_BP/
#           GO_gene_sets_summary.csv
#           GO_BP_gene_set_sizes.csv
#           zenith_Climate_PC1.csv
#           zenith_Climate_PC2.csv
#           zenith_Climate_PC3.csv
#           zenith_East_vs_Center.csv
#           zenith_West_vs_Center.csv
#           zenith_East_vs_West.csv
#           zenith_all_tests_combined.csv
#           zenith_all_tests_global_BH.csv
#           zenith_diagnostic_summary.csv
#           zenith_significant_FDR05.csv
#           zenith_significant_overlap_matrix.csv
#           zenith_climate_vs_range_overlap_summary.csv
#           zenith_shared_climate_range_terms.csv
#
#       GO_MF/
#           GO_MF_gene_sets_summary.csv
#           GO_MF_gene_set_sizes.csv
#           zenith_Climate_PC1.csv
#           ...
#
#       GO_CC/
#           GO_CC_gene_sets_summary.csv
#           GO_CC_gene_set_sizes.csv
#           zenith_Climate_PC1.csv
#           ...
#
#       InterPro/
#           InterPro_gene_sets_summary.csv
#           InterPro_gene_set_sizes.csv
#           ...
#
#       PFAM/
#           PFAM_gene_sets_summary.csv
#           PFAM_gene_set_sizes.csv
#           ...
#
#       COG/
#           COG_gene_sets_summary.csv
#           COG_gene_set_sizes.csv
#           ...
#
#       CAZy/
#           CAZy_gene_sets_summary.csv
#           CAZy_gene_set_sizes.csv
#           ...
#
#       sessionInfo.txt
#       analysis_parameters.txt
#
#
# ==============================================================================


###############################################################################
# 1. Libraries and reproducibility
###############################################################################

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(zenith)
  library(GSEABase)
  library(GO.db)
  library(AnnotationDbi)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
})


###############################################################################
# 2. Print software versions
###############################################################################

message("============================================================")
message("R / package versions")
message("============================================================")
message("R: ", R.version.string)
message("variancePartition: ", packageVersion("variancePartition"))
message("zenith: ", packageVersion("zenith"))
message("limma: ", packageVersion("limma"))
message("edgeR: ", packageVersion("edgeR"))
message("GSEABase: ", packageVersion("GSEABase"))
message("GO.db: ", packageVersion("GO.db"))
message("============================================================")


###############################################################################
# 3. User configuration
###############################################################################

# ---------------------------------------------------------------------------
# Main input locations
# ---------------------------------------------------------------------------

fit_dir <- "/Anastatica/degs/ECe_climate_range_model"

annotation_dir <- "/Anastatica/genome_annotation"

output_root <- "/Anastatica/zenith/ECe_climate_range_model"


# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------

fit_rds <- file.path(
  fit_dir,
  "ECe_climate_range_model_fit_contrast_ebayes.rds"
)

annotations_file <- file.path(
  annotation_dir,
  "Anastatica_hierochuntica_ANAS.annotations.xlsx"
)


# ---------------------------------------------------------------------------
# Biological/statistical terms to test
# ---------------------------------------------------------------------------

climate_coeffs <- c(
  "Climate_PC1",
  "Climate_PC2",
  "Climate_PC3"
)

range_contrasts <- c(
  "East_vs_Center",
  "West_vs_Center",
  "East_vs_West"
)

all_test_labels <- c(
  climate_coeffs,
  range_contrasts
)


# ---------------------------------------------------------------------------
# Gene-set analysis parameters
# ---------------------------------------------------------------------------

# Primary FDR threshold for interpreting zenith results.
fdr_cutoff <- 0.05

# Minimum number of genes required in a gene set.
min_geneset_size <- 10L


###############################################################################
# 4. Create output directories
###############################################################################

annotation_collections <- c(
  "GO_BP",
  "GO_MF",
  "GO_CC",
  "InterPro",
  "PFAM",
  "COG",
  "CAZy"
)

for (collection in annotation_collections) {

  dir.create(
    file.path(output_root, collection),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


###############################################################################
# 5. Helper functions
###############################################################################

# ---------------------------------------------------------------------------
# GO parser
# ---------------------------------------------------------------------------

parse_go_terms <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0 ||
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(character(0))
  }

  unique(
    unlist(
      str_extract_all(
        x,
        "GO:\\d{7}"
      ),
      use.names = FALSE
    )
  )
}


# ---------------------------------------------------------------------------
# PFAM parser
# ---------------------------------------------------------------------------

parse_pfam <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0 ||
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(character(0))
  }

  unique(
    unlist(
      str_extract_all(
        x,
        "PF\\d{5}"
      ),
      use.names = FALSE
    )
  )
}


# ---------------------------------------------------------------------------
# InterPro parser
# ---------------------------------------------------------------------------

parse_interpro <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0 ||
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(character(0))
  }

  hits <- str_extract_all(
    x,
    "(?<![A-Z0-9])IPR\\d{6}(?!\\d)"
  )[[1]]

  unique(hits)
}


# ---------------------------------------------------------------------------
# COG parser
# ---------------------------------------------------------------------------
#
# Extracts one-letter COG functional categories.
#
# Example:
#   G:(Glycosyl transferase);K:(...)
#
# gives:
#   G, K
#
# ---------------------------------------------------------------------------

parse_cog <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0 ||
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(character(0))
  }

  unique(
    str_extract_all(
      x,
      "(?<![A-Z])[A-Z](?=:)"
    )[[1]]
  )
}


# ---------------------------------------------------------------------------
# CAZy parser
# ---------------------------------------------------------------------------

parse_cazy <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0 ||
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(character(0))
  }

  y <- trimws(
    strsplit(
      x,
      ","
    )[[1]]
  )

  unique(
    y[
      y != ""
    ]
  )
}


# ---------------------------------------------------------------------------
# Build gene sets for one annotation column
# ---------------------------------------------------------------------------

build_one_gene_set_collection <- function(
    annotation_df,
    gene_universe,
    column,
    parser,
    min_size = 10L
) {

  if (!column %in% colnames(annotation_df)) {

    stop(
      "Annotation column '",
      column,
      "' is missing."
    )
  }

  out <- list()

  for (i in seq_len(nrow(annotation_df))) {

    gene <- annotation_df$GeneID[i]

    if (
      is.na(gene) ||
      gene == ""
    ) {
      next
    }

    terms <- parser(
      annotation_df[[column]][i]
    )

    if (length(terms) == 0) {
      next
    }

    for (term in terms) {

      out[[term]] <- c(
        out[[term]],
        gene
      )
    }
  }

  out <- lapply(
    out,
    function(x) {
      intersect(
        unique(x),
        gene_universe
      )
    }
  )

  out[
    vapply(
      out,
      length,
      integer(1)
    ) >= min_size
  ]
}


# ---------------------------------------------------------------------------
# Filter GO sets to Biological Process
# ---------------------------------------------------------------------------

filter_go_ontology <- function(
    go_sets,
    ontology
) {

  if (length(go_sets) == 0) {
    return(go_sets)
  }

  if (!ontology %in% c("BP", "MF", "CC")) {

    stop(
      "Invalid GO ontology: ",
      ontology,
      ". Expected one of BP, MF, CC."
    )
  }

  go_ids <- names(go_sets)

  ontology_table <- AnnotationDbi::select(
    GO.db,
    keys = go_ids,
    columns = "ONTOLOGY",
    keytype = "GOID"
  ) %>%
    distinct(
      GOID,
      .keep_all = TRUE
    )

  ontology_ids <- ontology_table %>%
    filter(
      ONTOLOGY == ontology
    ) %>%
    pull(
      GOID
    ) %>%
    unique()

  go_sets[
    intersect(
      go_ids,
      ontology_ids
    )
  ]
}


# ---------------------------------------------------------------------------
# Create GeneSetCollection
# ---------------------------------------------------------------------------

make_gene_set_collection <- function(
    gene_sets,
    collection_name
) {

  if (length(gene_sets) == 0) {
    stop(
      "No gene sets available for collection: ",
      collection_name
    )
  }

  GeneSetCollection(
    lapply(
      names(gene_sets),
      function(term) {

        GeneSet(
          geneIds = gene_sets[[term]],
          setName = paste(
            collection_name,
            term,
            sep = ":"
          ),
          geneIdType = NullIdentifier()
        )
      }
    )
  )
}


# ---------------------------------------------------------------------------
# Validate saved fit
# ---------------------------------------------------------------------------

check_fit <- function(
    fit,
    fit_path
) {

  if (!inherits(fit, "MArrayLM")) {

    stop(
      "Saved fit is not an MArrayLM object.\n",
      "Class: ",
      paste(
        class(fit),
        collapse = ", "
      ),
      "\nFile: ",
      fit_path
    )
  }

  if (is.null(fit$coef)) {
    stop("Saved fit contains no coefficient matrix.")
  }

  if (is.null(rownames(fit$coef))) {
    stop("Saved fit coefficient matrix has no gene IDs.")
  }

  if (is.null(fit$design)) {
    stop("Saved fit does not contain a design matrix.")
  }

  # zenith_gsa(..., inter.gene.cor = NA) requires residual information.
  if (is.null(fit$residuals)) {

    stop(
      "Saved fit does not contain residuals.\n",
      "The fit must retain residuals for zenith_gsa(..., ",
      "inter.gene.cor = NA).\n",
      "Refit the dream model with computeResiduals = TRUE."
    )
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Run one zenith test
# ---------------------------------------------------------------------------

run_zenith_test <- function(
    fit,
    gene_sets,
    label,
    collection_name,
    output_file
) {

  message("")
  message("------------------------------------------------------------")
  message(
    "Running zenith_gsa(): ",
    collection_name,
    " / ",
    label
  )
  message("------------------------------------------------------------")

  res <- zenith_gsa(
    fit = fit,
    geneSets = gene_sets,
    coefs = label,
    use.ranks = FALSE,
    n_genes_min = min_geneset_size,
    inter.gene.cor = NA,
    progressbar = TRUE
  )

  res$Test <- label
  res$Collection <- collection_name
  res$Model <- "ECe_climate_range_model"

  if (collection_name == "GO_BP") {

    res$Ontology <- "GO Biological Process"

  } else if (collection_name == "GO_MF") {

    res$Ontology <- "GO Molecular Function"

  } else if (collection_name == "GO_CC") {

    res$Ontology <- "GO Cellular Component"

  } else {

    res$Ontology <- NA_character_

  }

  expected <- c(
    "Geneset",
    "NGenes",
    "delta",
    "se",
    "PValue",
    "Direction",
    "FDR"
  )

  missing <- setdiff(
    expected,
    colnames(res)
  )

  if (length(missing) > 0) {

    stop(
      "Zenith result for ",
      collection_name,
      " / ",
      label,
      " is missing expected columns: ",
      paste(
        missing,
        collapse = ", "
      ),
      "\nAvailable columns: ",
      paste(
        colnames(res),
        collapse = ", "
      )
    )
  }

  res <- res %>%
    arrange(
      FDR,
      PValue
    )

  write.csv(
    res,
    output_file,
    row.names = FALSE
  )

  message(
    "Gene sets tested: ",
    nrow(res)
  )

  message(
    "FDR < 0.05: ",
    sum(
      res$FDR < fdr_cutoff,
      na.rm = TRUE
    )
  )

  message(
    "Minimum raw P: ",
    format(
      min(
        res$PValue,
        na.rm = TRUE
      ),
      scientific = TRUE
    )
  )

  message(
    "Minimum FDR: ",
    format(
      min(
        res$FDR,
        na.rm = TRUE
      ),
      scientific = TRUE
    )
  )

  invisible(res)
}


###############################################################################
# 6. Check input files
###############################################################################

required_files <- c(
  fit_rds,
  annotations_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "The following required input files do not exist:\n  ",
    paste(
      missing_files,
      collapse = "\n  "
    )
  )
}


###############################################################################
# 7. Load exact saved dream/eBayes fit
###############################################################################

message("")
message("============================================================")
message("Loading saved combined dream/eBayes fit")
message("============================================================")
message("File: ", fit_rds)

fit <- readRDS(
  fit_rds
)

check_fit(
  fit = fit,
  fit_path = fit_rds
)

message(
  "Genes in fit: ",
  nrow(fit$coef)
)

message(
  "Coefficients in fit: ",
  ncol(fit$coef)
)

message("Coefficient names:")
print(
  colnames(fit$coef)
)

message("============================================================")


###############################################################################
# 8. Inspect saved design and coefficient names
###############################################################################

message("")
message("Saved design matrix columns:")
print(
  colnames(fit$design)
)

message("")
message("Saved coefficient names:")
print(
  colnames(fit$coef)
)


###############################################################################
# 9. Verify climate coefficients
###############################################################################

available_coefficients <- colnames(
  fit$coef
)

missing_climate <- setdiff(
  climate_coeffs,
  available_coefficients
)

if (length(missing_climate) > 0) {

  stop(
    "The saved full model does not contain the requested climate coefficients:\n  ",
    paste(
      missing_climate,
      collapse = "\n  "
    ),
    "\n\nAvailable coefficients are:\n  ",
    paste(
      available_coefficients,
      collapse = "\n  "
    )
  )
}


###############################################################################
# 10. Verify range contrasts
###############################################################################

missing_range <- setdiff(
  range_contrasts,
  available_coefficients
)

if (length(missing_range) > 0) {

  stop(
    "The saved full model does not contain all requested range contrasts:\n  ",
    paste(
      missing_range,
      collapse = "\n  "
    ),
    "\n\nAvailable coefficients are:\n  ",
    paste(
      available_coefficients,
      collapse = "\n  "
    ),
    "\n\nCheck that the saved fit was generated with the intended ",
    "East/Centre/West contrasts."
  )
}


###############################################################################
# 11. Define exact gene universe from saved dream fit
###############################################################################

message("")
message("============================================================")
message("Defining exact gene universe from saved dream fit")
message("============================================================")

fit_genes <- rownames(
  fit$coef
)

if (
  is.null(fit_genes) ||
  length(fit_genes) == 0
) {

  stop(
    "Could not obtain gene IDs from rownames(fit$coef)."
  )
}

fit_genes <- as.character(
  fit_genes
)

if (anyDuplicated(fit_genes)) {

  stop(
    "The saved dream fit contains duplicated gene IDs."
  )
}

gene_universe <- fit_genes

message(
  "Genes in exact saved dream fit universe: ",
  length(gene_universe)
)


###############################################################################
# 12. Validate saved full-model fixed effects
###############################################################################

required_fixed_terms <- c(
  "ECe_mean",
  "Climate_PC1",
  "Climate_PC2",
  "Climate_PC3"
)

design_columns <- colnames(
  fit$design
)

missing_fixed_terms <- setdiff(
  required_fixed_terms,
  design_columns
)

if (length(missing_fixed_terms) > 0) {

  stop(
    "The saved fit design matrix is missing required fixed-effect columns:\n  ",
    paste(
      missing_fixed_terms,
      collapse = "\n  "
    ),
    "\n\nAvailable design columns:\n  ",
    paste(
      design_columns,
      collapse = "\n  "
    )
  )
}

range_design_columns <- grep(
  "^Site_of_distribution_range",
  design_columns,
  value = TRUE
)

if (length(range_design_columns) < 2) {

  stop(
    "Could not identify the expected three-level distribution-range terms ",
    "in the saved design.\n",
    "Detected range columns: ",
    paste(
      range_design_columns,
      collapse = ", "
    ),
    "\nAvailable design columns:\n  ",
    paste(
      design_columns,
      collapse = "\n  "
    )
  )
}

message("Detected range design columns:")
print(
  range_design_columns
)


###############################################################################
# 13. Load and validate annotation
###############################################################################

message("")
message("============================================================")
message("Loading Anastatica annotation")
message("============================================================")

annotations <- read_excel(
  annotations_file,
  sheet = 1,
  col_types = "text"
)

message("Annotation columns:")
print(
  colnames(annotations)
)

required_annotation_columns <- c(
  "GeneID",
  "GO Terms",
  "PFAM",
  "InterPro",
  "COG",
  "CAZyme"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annotations)
)

if (length(missing_annotation_columns) > 0) {

  stop(
    "Annotation file is missing required columns:\n  ",
    paste(
      missing_annotation_columns,
      collapse = "\n  "
    ),
    "\n\nAvailable columns:\n  ",
    paste(
      colnames(annotations),
      collapse = "\n  "
    )
  )
}

annotations$GeneID <- as.character(
  annotations$GeneID
)

message("")
message("First annotation row:")
print(
  annotations[
    1,
    c(
      "GeneID",
      "PFAM",
      "InterPro",
      "COG",
      "GO Terms",
      "CAZyme"
    ),
    drop = FALSE
  ]
)

annotation_match_fraction <- mean(
  gene_universe %in% annotations$GeneID
)

message(
  "Fraction of fit genes with direct GeneID annotation match: ",
  round(
    annotation_match_fraction,
    4
  )
)

if (annotation_match_fraction < 0.50) {

  warning(
    "Less than 50% of the model gene universe directly matches ",
    "annotation GeneID values. Check gene-ID conventions carefully."
  )
}


###############################################################################
# 14. Restrict annotations to the exact dream gene universe
###############################################################################

annotations_filtered <- annotations %>%
  filter(
    GeneID %in% gene_universe
  )

message(
  "Annotation rows matching exact fit gene universe: ",
  nrow(annotations_filtered)
)


###############################################################################
# 15. Annotation parsing diagnostics
###############################################################################

message("")
message("============================================================")
message("Annotation parsing diagnostics")
message("============================================================")

# GO

go_annotation_count <- sum(
  vapply(
    annotations_filtered$`GO Terms`,
    function(x) {
      length(
        parse_go_terms(x)
      ) > 0
    },
    logical(1)
  )
)

message(
  "Genes with at least one parsed GO ID: ",
  go_annotation_count,
  " / ",
  nrow(annotations_filtered)
)


# PFAM

pfam_annotation_count <- sum(
  vapply(
    annotations_filtered$PFAM,
    function(x) {
      length(
        parse_pfam(x)
      ) > 0
    },
    logical(1)
  )
)

message(
  "Genes with at least one parsed PFAM ID: ",
  pfam_annotation_count,
  " / ",
  nrow(annotations_filtered)
)


# InterPro

interpro_annotation_count <- sum(
  vapply(
    annotations_filtered$InterPro,
    function(x) {
      length(
        parse_interpro(x)
      ) > 0
    },
    logical(1)
  )
)

message(
  "Genes with at least one parsed InterPro ID: ",
  interpro_annotation_count,
  " / ",
  nrow(annotations_filtered)
)

example_interpro <- unique(
  unlist(
    lapply(
      annotations_filtered$InterPro,
      parse_interpro
    )
  )
)

message(
  "Unique InterPro IDs parsed: ",
  length(example_interpro)
)

if (length(example_interpro) > 0) {

  message(
    "First parsed InterPro IDs: ",
    paste(
      head(
        example_interpro,
        10
      ),
      collapse = ", "
    )
  )
}


# COG

cog_annotation_count <- sum(
  vapply(
    annotations_filtered$COG,
    function(x) {
      length(
        parse_cog(x)
      ) > 0
    },
    logical(1)
  )
)

message(
  "Genes with at least one parsed COG category: ",
  cog_annotation_count,
  " / ",
  nrow(annotations_filtered)
)


# CAZy

cazy_annotation_count <- sum(
  vapply(
    annotations_filtered$CAZyme,
    function(x) {
      length(
        parse_cazy(x)
      ) > 0
    },
    logical(1)
  )
)

message(
  "Genes with at least one parsed CAZy annotation: ",
  cazy_annotation_count,
  " / ",
  nrow(annotations_filtered)
)


###############################################################################
# 16. Build all annotation gene sets
###############################################################################

message("")
message("============================================================")
message("Building annotation gene sets")
message("============================================================")

# ---------------------------------------------------------------------------
# GO
# ---------------------------------------------------------------------------

go_sets_all <- build_one_gene_set_collection(
  annotation_df = annotations_filtered,
  gene_universe = gene_universe,
  column = "GO Terms",
  parser = parse_go_terms,
  min_size = min_geneset_size
)

message(
  "GO sets before ontology filtering: ",
  length(go_sets_all)
)


# GO Biological Process

go_sets_bp <- filter_go_ontology(
  go_sets_all,
  ontology = "BP"
)

message(
  "GO Biological Process sets: ",
  length(go_sets_bp)
)


# GO Molecular Function

go_sets_mf <- filter_go_ontology(
  go_sets_all,
  ontology = "MF"
)

message(
  "GO Molecular Function sets: ",
  length(go_sets_mf)
)


# GO Cellular Component

go_sets_cc <- filter_go_ontology(
  go_sets_all,
  ontology = "CC"
)

message(
  "GO Cellular Component sets: ",
  length(go_sets_cc)
)

# ---------------------------------------------------------------------------
# PFAM
# ---------------------------------------------------------------------------

pfam_sets <- build_one_gene_set_collection(
  annotation_df = annotations_filtered,
  gene_universe = gene_universe,
  column = "PFAM",
  parser = parse_pfam,
  min_size = min_geneset_size
)

message(
  "PFAM sets: ",
  length(pfam_sets)
)


# ---------------------------------------------------------------------------
# InterPro
# ---------------------------------------------------------------------------

interpro_sets <- build_one_gene_set_collection(
  annotation_df = annotations_filtered,
  gene_universe = gene_universe,
  column = "InterPro",
  parser = parse_interpro,
  min_size = min_geneset_size
)

message(
  "InterPro sets: ",
  length(interpro_sets)
)


# ---------------------------------------------------------------------------
# COG
# ---------------------------------------------------------------------------

cog_sets <- build_one_gene_set_collection(
  annotation_df = annotations_filtered,
  gene_universe = gene_universe,
  column = "COG",
  parser = parse_cog,
  min_size = min_geneset_size
)

message(
  "COG sets: ",
  length(cog_sets)
)


# ---------------------------------------------------------------------------
# CAZy
# ---------------------------------------------------------------------------

cazy_sets <- build_one_gene_set_collection(
  annotation_df = annotations_filtered,
  gene_universe = gene_universe,
  column = "CAZyme",
  parser = parse_cazy,
  min_size = min_geneset_size
)

message(
  "CAZy sets: ",
  length(cazy_sets)
)


###############################################################################
# 17. Store annotation collections
###############################################################################

gene_sets_by_collection <- list(
  GO_BP = go_sets_bp,
  GO_MF = go_sets_mf,
  GO_CC = go_sets_cc,
  InterPro = interpro_sets,
  PFAM = pfam_sets,
  COG = cog_sets,
  CAZy = cazy_sets
)

# Remove completely empty collections.

gene_sets_by_collection <- gene_sets_by_collection[
  vapply(
    gene_sets_by_collection,
    length,
    integer(1)
  ) > 0
]

message("")
message("Final annotation collections:")
print(
  data.frame(
    Collection = names(gene_sets_by_collection),
    N_gene_sets = vapply(
      gene_sets_by_collection,
      length,
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
)


###############################################################################
# 18. Summarize each annotation collection
###############################################################################

message("")
message("============================================================")
message("Writing annotation summaries")
message("============================================================")

collection_summary_list <- list()

for (collection_name in names(gene_sets_by_collection)) {

  sets <- gene_sets_by_collection[[collection_name]]

  gene_counts <- lengths(
    sets
  )

  summary_df <- data.frame(
    Collection = collection_name,
    N_gene_sets = length(sets),
    N_unique_annotated_genes = length(
      unique(
        unlist(
          sets,
          use.names = FALSE
        )
      )
    ),
    Minimum_set_size = min(
      gene_counts
    ),
    Median_set_size = median(
      gene_counts
    ),
    Maximum_set_size = max(
      gene_counts
    ),
    stringsAsFactors = FALSE
  )

  collection_summary_list[[collection_name]] <- summary_df

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  write.csv(
    summary_df,
    file.path(
      collection_dir,
      paste0(
        collection_name,
        "_gene_sets_summary.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      Collection = collection_name,
      GeneSet = names(sets),
      N_genes = lengths(sets),
      stringsAsFactors = FALSE
    ),
    file.path(
      collection_dir,
      paste0(
        collection_name,
        "_gene_set_sizes.csv"
      )
    ),
    row.names = FALSE
  )

  print(
    summary_df
  )
}

all_collection_summary <- bind_rows(
  collection_summary_list
)

write.csv(
  all_collection_summary,
  file.path(
    output_root,
    "all_annotation_collections_summary.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 19. Create GeneSetCollections
###############################################################################

message("")
message("============================================================")
message("Creating GeneSetCollections")
message("============================================================")

gene_set_collections <- list()

for (collection_name in names(gene_sets_by_collection)) {

  gene_set_collections[[collection_name]] <-
    make_gene_set_collection(
      gene_sets = gene_sets_by_collection[[collection_name]],
      collection_name = collection_name
    )

  message(
    collection_name,
    ": ",
    length(
      gene_set_collections[[collection_name]]
    ),
    " gene sets"
  )
}


###############################################################################
# 20. Run all six zenith tests for every annotation collection
###############################################################################
#
# Each annotation collection is analyzed independently.
#
# This means:
#
#   GO_BP:
#       6 zenith tests
#
#   InterPro:
#       6 zenith tests
#
#   PFAM:
#       6 zenith tests
#
#   COG:
#       6 zenith tests
#
#   CAZy:
#       6 zenith tests
#
# Total:
#       30 zenith analyses.
#
###############################################################################

all_zenith_results <- list()

for (collection_name in names(gene_set_collections)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  gene_sets <- gene_set_collections[[collection_name]]

  collection_results <- list()

  message("")
  message("============================================================")
  message(
    "ANNOTATION COLLECTION: ",
    collection_name
  )
  message("============================================================")

  for (test_label in all_test_labels) {

    output_file <- file.path(
      collection_dir,
      paste0(
        "zenith_",
        test_label,
        ".csv"
      )
    )

    result <- run_zenith_test(
      fit = fit,
      gene_sets = gene_sets,
      label = test_label,
      collection_name = collection_name,
      output_file = output_file
    )

    collection_results[[test_label]] <- result
  }

  all_zenith_results[[collection_name]] <-
    collection_results
}


###############################################################################
# 21. Combine results separately within each annotation collection
###############################################################################

message("")
message("============================================================")
message("Combining zenith results by annotation collection")
message("============================================================")

combined_results_by_collection <- list()

for (collection_name in names(all_zenith_results)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  combined <- bind_rows(
    all_zenith_results[[collection_name]]
  )

  combined_results_by_collection[[collection_name]] <-
    combined

  write.csv(
    combined,
    file.path(
      collection_dir,
      "zenith_all_tests_combined.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 22. Secondary global BH adjustment
###############################################################################
#
# IMPORTANT:
#
# The standard zenith FDR remains the PRIMARY reported quantity.
#
# The following adjustment pools all raw P-values across the six tests WITHIN
# EACH ANNOTATION COLLECTION.
#
# It therefore asks a more conservative question:
#
#   "Which gene sets remain significant when all six tests within this
#    annotation collection are treated as one family?"
#
# It does NOT pool GO + InterPro + PFAM + COG + CAZy together.
#
###############################################################################

message("")
message("============================================================")
message("Calculating secondary global BH corrections")
message("============================================================")

global_bh_results_by_collection <- list()

for (collection_name in names(combined_results_by_collection)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  combined <- combined_results_by_collection[[collection_name]]

  combined <- combined %>%
    mutate(
      Global_BH_FDR_all_six_tests = p.adjust(
        PValue,
        method = "BH"
      )
    )

  global_bh_results_by_collection[[collection_name]] <-
    combined

  write.csv(
    combined,
    file.path(
      collection_dir,
      "zenith_all_tests_global_BH.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 23. Diagnostic summary for every annotation collection
###############################################################################

message("")
message("============================================================")
message("Creating diagnostic summaries")
message("============================================================")

diagnostic_by_collection <- list()

for (collection_name in names(global_bh_results_by_collection)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  combined <- global_bh_results_by_collection[[collection_name]]

  diagnostic <- combined %>%
    group_by(
      Test
    ) %>%
    summarise(
      Collection = first(Collection),
      N_gene_sets_tested = n(),

      N_FDR05 = sum(
        FDR < fdr_cutoff,
        na.rm = TRUE
      ),

      N_global_BH05 = sum(
        Global_BH_FDR_all_six_tests < fdr_cutoff,
        na.rm = TRUE
      ),

      Minimum_PValue = min(
        PValue,
        na.rm = TRUE
      ),

      Minimum_FDR = min(
        FDR,
        na.rm = TRUE
      ),

      Minimum_Global_BH_FDR = min(
        Global_BH_FDR_all_six_tests,
        na.rm = TRUE
      ),

      .groups = "drop"
    )

  diagnostic_by_collection[[collection_name]] <-
    diagnostic

  write.csv(
    diagnostic,
    file.path(
      collection_dir,
      "zenith_diagnostic_summary.csv"
    ),
    row.names = FALSE
  )

  print(diagnostic)
}

all_diagnostics <- bind_rows(
  diagnostic_by_collection
)

write.csv(
  all_diagnostics,
  file.path(
    output_root,
    "all_annotation_collections_diagnostic_summary.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 24. Significant gene-set tables
###############################################################################

message("")
message("============================================================")
message("Writing significant gene-set tables")
message("============================================================")

significant_by_collection <- list()

for (collection_name in names(global_bh_results_by_collection)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  combined <- global_bh_results_by_collection[[collection_name]]

  significant <- combined %>%
    filter(
      FDR < fdr_cutoff
    ) %>%
    arrange(
      Test,
      FDR,
      PValue
    )

  significant_by_collection[[collection_name]] <-
    significant

  write.csv(
    significant,
    file.path(
      collection_dir,
      "zenith_significant_FDR05.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 25. Descriptive overlap of significant terms WITHIN each collection
###############################################################################
#
# This is descriptive only.
#
# A gene set significant for climate and range does NOT imply a formal
# climate-by-range interaction.
#
###############################################################################

message("")
message("============================================================")
message("Creating within-collection significant-set overlap tables")
message("============================================================")

for (collection_name in names(global_bh_results_by_collection)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  combined <- global_bh_results_by_collection[[collection_name]]

  significant_sets <- lapply(
    unique(
      combined$Test
    ),
    function(test_name) {

      combined %>%
        filter(
          Test == test_name,
          FDR < fdr_cutoff
        ) %>%
        pull(
          Geneset
        ) %>%
        unique()
    }
  )

  names(significant_sets) <- unique(
    combined$Test
  )

  all_significant_terms <- sort(
    unique(
      unlist(
        significant_sets,
        use.names = FALSE
      )
    )
  )

  overlap_table <- data.frame(
    GeneSet = all_significant_terms,
    stringsAsFactors = FALSE
  )

  for (test_name in names(significant_sets)) {

    overlap_table[[test_name]] <-
      all_significant_terms %in%
      significant_sets[[test_name]]
  }

  write.csv(
    overlap_table,
    file.path(
      collection_dir,
      "zenith_significant_overlap_matrix.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 26. Climate-vs-range overlap summaries
###############################################################################
#
# Again, these are DESCRIPTIVE ONLY.
#
# Shared gene sets do not constitute evidence for a formal climate-by-range
# interaction.
#
###############################################################################

message("")
message("============================================================")
message("Creating climate-vs-range overlap summaries")
message("============================================================")

for (collection_name in names(global_bh_results_by_collection)) {

  collection_dir <- file.path(
    output_root,
    collection_name
  )

  combined <- global_bh_results_by_collection[[collection_name]]

  significant_sets <- lapply(
    unique(
      combined$Test
    ),
    function(test_name) {

      combined %>%
        filter(
          Test == test_name,
          FDR < fdr_cutoff
        ) %>%
        pull(
          Geneset
        ) %>%
        unique()
    }
  )

  names(significant_sets) <- unique(
    combined$Test
  )

  climate_term_union <- sort(
    unique(
      unlist(
        significant_sets[
          intersect(
            climate_coeffs,
            names(significant_sets)
          )
        ],
        use.names = FALSE
      )
    )
  )

  range_term_union <- sort(
    unique(
      unlist(
        significant_sets[
          intersect(
            range_contrasts,
            names(significant_sets)
          )
        ],
        use.names = FALSE
      )
    )
  )

  shared_terms <- intersect(
    climate_term_union,
    range_term_union
  )

  union_terms <- union(
    climate_term_union,
    range_term_union
  )

  term_overlap_summary <- data.frame(
    Collection = collection_name,

    N_climate_gene_sets_FDR05 =
      length(climate_term_union),

    N_range_gene_sets_FDR05 =
      length(range_term_union),

    N_shared_gene_sets =
      length(shared_terms),

    Jaccard =
      ifelse(
        length(union_terms) == 0,
        NA_real_,
        length(shared_terms) /
          length(union_terms)
      ),

    stringsAsFactors = FALSE
  )

  write.csv(
    term_overlap_summary,
    file.path(
      collection_dir,
      "zenith_climate_vs_range_overlap_summary.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      GeneSet = shared_terms,
      stringsAsFactors = FALSE
    ),
    file.path(
      collection_dir,
      "zenith_shared_climate_range_terms.csv"
    ),
    row.names = FALSE
  )
}


###############################################################################
# 27. Create a cross-collection summary of gene-set counts
###############################################################################

message("")
message("============================================================")
message("Creating cross-collection summary")
message("============================================================")

cross_collection_summary <- all_diagnostics %>%
  arrange(
    Collection,
    Test
  )

write.csv(
  cross_collection_summary,
  file.path(
    output_root,
    "all_collections_zenith_summary.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 28. Save reproducibility information
###############################################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_root,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    "================================================================",
    "ZENITH MULTI-ANNOTATION ANALYSIS PARAMETERS",
    "================================================================",
    "",
    "Species:",
    "Anastatica hierochuntica",
    "",
    "Input model:",
    fit_rds,
    "",
    "Input annotation:",
    annotations_file,
    "",
    "Model represented by saved fit:",
    "~ ECe_mean + Climate_PC1 + Climate_PC2 + Climate_PC3 +",
    "  Site_of_distribution_range + (1 | Sequencing_batch) +",
    "  (1 | Population)",
    "",
    "Primary gene-set method:",
    "zenith_gsa()",
    "",
    "inter.gene.cor:",
    "NA",
    "",
    "use.ranks:",
    "FALSE",
    "",
    paste0(
      "Minimum gene-set size: ",
      min_geneset_size
    ),
    "",
    paste0(
      "Primary within-test FDR threshold: ",
      fdr_cutoff
    ),
    "",
    "Annotation collections:",
    paste(
      names(gene_sets_by_collection),
      collapse = ", "
    ),
    "",
    "Primary annotation collections:",
    "GO Biological Process",
    "GO Molecular Function",
    "GO Cellular Component",
    "",
    "Secondary annotation collections:",
    "InterPro",
    "PFAM",
    "COG",
    "CAZy",
    "",
    "Climate tests:",
    paste(
      climate_coeffs,
      collapse = ", "
    ),
    "",
    "Range contrasts:",
    paste(
      range_contrasts,
      collapse = ", "
    ),
    "",
    "FDR interpretation:",
    "Standard zenith FDR is controlled within each test and annotation",
    "collection.",
    "",
    "Secondary global BH:",
    "Raw P-values are pooled across the six tests within each annotation",
    "collection only.",
    "",
    "Annotation parsing:",
    "GO       = GO:1234567 identifiers",
    "PFAM     = PF12345 identifiers",
    "InterPro = IPR123456 identifiers",
    "COG      = one-letter COG categories",
    "CAZy     = comma-separated CAZyme identifiers",
    "",
    "Important:",
    "Overlap between annotation collections is descriptive and does not",
    "represent independent statistical evidence.",
    "",
    "================================================================"
  ),
  con = file.path(
    output_root,
    "analysis_parameters.txt"
  )
)


###############################################################################
# 29. Final message
###############################################################################

message("")
message("============================================================")
message("ZENITH MULTI-ANNOTATION ANALYSIS COMPLETED")
message("============================================================")

message("")
message("Output root:")
message(output_root)

message("")
message("Annotation collections analyzed:")

for (collection_name in names(gene_sets_by_collection)) {

  message(
    "  ",
    collection_name,
    ": ",
    length(
      gene_sets_by_collection[[collection_name]]
    ),
    " gene sets"
  )
}

message("")
message("Each collection contains:")
message("  - six individual zenith result files")
message("  - combined six-test results")
message("  - secondary global-BH results")
message("  - diagnostic summary")
message("  - significant FDR < 0.05 results")
message("  - significant-set overlap matrix")
message("  - climate-vs-range overlap summary")

message("")
message("Primary GO collections:")
message("  GO_BP")
message("  GO_MF")
message("  GO_CC")

message("")
message("Secondary collections:")
message("  InterPro")
message("  PFAM")
message("  COG")
message("  CAZy")

message("")
message("============================================================")