#!/usr/bin/env R

###############################################################################
# Author   : Iris Sammarco
# Date     : 09/09/2026
#
# GO / FUNCTIONAL ENRICHMENT ANALYSIS
# Anastatica hierochuntica
#
# PURPOSE
# -------
# This script performs functional enrichment analysis of the gene classes
# obtained from the multivariable differential-expression analysis:
#
#   1. climate_only
#   2. range_only
#   3. both
#
# The "neither" class is NOT treated as a biological target set. Instead,
# all genes that were eligible for the H4 analysis are used as the statistical
# background/universe.
#
# PRIMARY ANALYSES
# ----------------
# 1. GO Biological Process (GO-BP) over-representation analysis
# 2. GO Molecular Function (GO-MF) over-representation analysis
# 3. GO Cellular Component (GO-CC) over-representation analysis
# 4. InterPro domain enrichment
# 5. PFAM domain enrichment
#
# The analysis uses the species-specific functional annotations provided in:
#
#   Anastatica_hierochuntica_ANAS.annotations.xlsx
#
# rather than relying on annotations transferred from another species.
#
# IMPORTANT STATISTICAL PRINCIPLE
# --------------------------------
# The H4 analysis universe contains 16,872 eligible genes. By default,
# enrichment uses the subset of that universe that has at least one annotation
# in the collection being tested (annotation-aware background). This is the
# recommended setting because unannotated genes cannot be classified as members
# of a functional category. Set BACKGROUND_MODE = "all" to reproduce the
# original all-H4-universe approach.
#
# This is important because enrichment asks:
#
#   "Are particular functions over-represented among H4-associated genes
#    compared with the genes that could have been selected in this analysis?"
#
# rather than:
#
#   "Are they over-represented compared with every gene in the genome?"
#
# For each target gene class, enrichment is tested using Fisher's exact test.
# By default, the statistical background is annotation-aware: genes without
# any annotation in the relevant collection are excluded from that collection's
# background. P-values are corrected using Benjamini-Hochberg FDR separately
# within each target class and annotation collection.
#
# IMPORTANT: the minimum-overlap rule is applied only after BH correction and
# therefore does not alter the multiple-testing correction.
#
# TARGET CLASSES
# --------------
# climate_only
# range_only
# both
#
# The "neither" genes are retained in the background universe but are not
# tested as a target gene set.
#
# ANNOTATIONS AVAILABLE IN THE INPUT FILE
# ----------------------------------------
# GO Terms
# PFAM
# InterPro
# COG
# CAZyme
# EggNog
#
# This script focuses on GO-BP, GO-MF, GO-CC, InterPro and PFAM.
#
# GO terms are separated into the three standard Gene Ontology namespaces:
#
#   GO-BP = Biological Process
#   GO-MF = Molecular Function
#   GO-CC = Cellular Component
#
# The GO ontology assignment is obtained from GO.db.
#
# KEGG
# ----
# The supplied annotation structure shown by the user does not contain
# KEGG pathway or KO identifiers. Therefore this script does NOT fabricate
# KEGG mappings. KEGG analysis should be performed separately if a valid
# KEGG/KO annotation or mapping becomes available.
#
# INPUT FILES
# -----------
# 1. H4 gene classification:
#
# /Anastatica/climate_range_variance_analysis/tables/gene_classification.csv
#
# Expected columns:
#   GeneID
#   Climate_associated
#   Range_associated
#   Category
#
# 2. Anastatica functional annotation:
#
# /Anastatica/genome_annotation/Anastatica_hierochuntica_ANAS.annotations.xlsx
#
# Important annotation columns:
#   GeneID
#   GO Terms
#   PFAM
#   InterPro
#
# OUTPUT DIRECTORY
# ----------------
# A new directory is created:
#
# /Anastatica/climate_range_variance_analysis/GO_pathway_analysis/
#
# Output structure:
#
# GO_pathway_analysis/
#
#   GO_BP/
#       GO_BP_climate_only.csv
#       GO_BP_range_only.csv
#       GO_BP_both.csv
#       GO_BP_climate_only_significant.csv
#       GO_BP_range_only_significant.csv
#       GO_BP_both_significant.csv
#
#   GO_MF/
#       GO_MF_climate_only.csv
#       GO_MF_range_only.csv
#       GO_MF_both.csv
#       GO_MF_climate_only_significant.csv
#       GO_MF_range_only_significant.csv
#       GO_MF_both_significant.csv
#
#   GO_CC/
#       GO_CC_climate_only.csv
#       GO_CC_range_only.csv
#       GO_CC_both.csv
#       GO_CC_climate_only_significant.csv
#       GO_CC_range_only_significant.csv
#       GO_CC_both_significant.csv
#
#   InterPro/
#       InterPro_climate_only.csv
#       InterPro_range_only.csv
#       InterPro_both.csv
#       ...
#
#   PFAM/
#       PFAM_climate_only.csv
#       PFAM_range_only.csv
#       PFAM_both.csv
#       ...
#
#   annotation_coverage.csv
#   InterPro_term_descriptions.csv
#   PFAM_term_descriptions.csv
#   InterPro_metadata_cache.csv
#   PFAM_metadata_cache.csv
#   H4_enrichment_summary.csv
#   annotation_row_diagnostics.csv
#   H4_all_enrichment_results.csv
#
# INSTALLATION
# ------------
# Required Bioconductor packages:
#
#   BiocManager
#   GO.db
#
# Required CRAN packages:
#
#   readxl
#   dplyr
#   tidyr
#   stringr
#   ggplot2
#
# If packages are not installed, run ONCE in R:
#
# install.packages(c(
#   "BiocManager",
#   "readxl",
#   "dplyr",
#   "tidyr",
#   "stringr",
#   "ggplot2"
# ))
#
# BiocManager::install("GO.db")
#
# This script assumes R >= 4.3.
#
###############################################################################


###############################################################################
# 1. LIBRARIES
###############################################################################

suppressPackageStartupMessages({

  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(GO.db)
  library(AnnotationDbi)
  library(jsonlite)

})


###############################################################################
# 2. PARAMETERS
###############################################################################

# ---------------------------------------------------------------------------
# Main directory
# ---------------------------------------------------------------------------

root <- paste0(
  "/Anastatica/",
  "climate_range_variance_analysis"
)


# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------

classification_file <- file.path(
  root,
  "tables",
  "gene_classification.csv"
)


annotation_file <- paste0(
  "/Anastatica/",
  "genome_annotation/",
  "Anastatica_hierochuntica_ANAS.annotations.xlsx"
)


# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------

output_root <- file.path(
  root,
  "GO_pathway_analysis"
)


# ---------------------------------------------------------------------------
# Enrichment parameters
# ---------------------------------------------------------------------------

# Minimum number of annotated genes in a functional category required for
# testing.
#
# This avoids unstable enrichment tests based on extremely small categories.

MIN_CATEGORY_SIZE <- 10L


# Minimum number of genes from the target list required to overlap an
# annotation category for that category to be retained in the final results.

MIN_TARGET_OVERLAP <- 3L


# FDR threshold for significant enrichment.

FDR_CUTOFF <- 0.05

# Statistical background mode.
# "annotated_by_collection" is recommended for enrichment because genes with
# no annotation in a given collection cannot be classified into its categories.
# "all" reproduces the original analysis using all 16,872 genes.
BACKGROUND_MODE <- "annotated_by_collection"

# InterPro/PFAM API settings.
# Metadata are cached locally so that subsequent runs do not re-query the API.
API_SLEEP_SECONDS <- 0.2


# Target classes.

target_classes <- c(
  "climate_only",
  "range_only",
  "both"
)


# Annotation collections analysed.

annotation_collections <- c(
  "GO_BP",
  "GO_MF",
  "GO_CC",
  "InterPro",
  "PFAM"
)


###############################################################################
# 3. CREATE OUTPUT DIRECTORIES
###############################################################################

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

for (collection in annotation_collections) {

  dir.create(
    file.path(
      output_root,
      collection
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )

}


###############################################################################
# 4. CHECK INPUT FILES
###############################################################################

if (!file.exists(classification_file)) {

  stop(
    "Classification file not found:\n",
    classification_file
  )

}

if (!file.exists(annotation_file)) {

  stop(
    "Annotation file not found:\n",
    annotation_file
  )

}


###############################################################################
# 5. LOAD H4 GENE CLASSIFICATION
###############################################################################

message("")
message("============================================================")
message("Loading gene classification")
message("============================================================")

h4 <- read.csv(
  classification_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_h4_columns <- c(
  "GeneID",
  "Climate_associated",
  "Range_associated",
  "Category"
)


missing_h4_columns <- setdiff(
  required_h4_columns,
  colnames(h4)
)


if (length(missing_h4_columns) > 0) {

  stop(
    "H4 classification file is missing columns:\n",
    paste(
      missing_h4_columns,
      collapse = ", "
    )
  )

}


h4$GeneID <- as.character(h4$GeneID)


# Check duplicated genes.

if (anyDuplicated(h4$GeneID)) {

  stop(
    "Duplicated GeneID values detected in H4 classification."
  )

}


# Exact H4 universe.

gene_universe <- unique(
  h4$GeneID
)


message(
  "Genes in H4 universe: ",
  length(gene_universe)
)


message("")
message("H4 class counts:")

print(
  h4 %>%
    count(Category) %>%
    mutate(
      Proportion = n / sum(n)
    )
)


###############################################################################

# 6. LOAD ANNOTATION FILE

###############################################################################

message("")
message("============================================================")
message("Loading Anastatica functional annotations")
message("============================================================")

annotations_raw <- read_excel(
annotation_file,
sheet = 1,
col_types = "text"
)

required_annotation_columns <- c(
"GeneID",
"GO Terms",
"PFAM",
"InterPro"
)

missing_annotation_columns <- setdiff(
required_annotation_columns,
colnames(annotations_raw)
)

if (length(missing_annotation_columns) > 0) {

stop(
"Annotation file is missing required columns:\n",
paste(
missing_annotation_columns,
collapse = ", "
)
)

}

annotations_raw$GeneID <- trimws(
as.character(
annotations_raw$GeneID
)
)

###############################################################################

# 6A. CHECK TRANSCRIPT-LEVEL DUPLICATION

###############################################################################

message("")
message("Checking annotation structure")

n_annotation_rows <- nrow(
annotations_raw
)

n_annotation_genes <- length(
unique(
annotations_raw$GeneID
)
)

n_duplicated_geneids <- sum(
duplicated(
annotations_raw$GeneID
)
)

message(
"Annotation rows: ",
n_annotation_rows
)

message(
"Unique GeneIDs: ",
n_annotation_genes
)

message(
"Additional rows belonging to duplicated GeneIDs: ",
n_duplicated_geneids
)

if (n_duplicated_geneids > 0) {

message(
"Multiple annotation rows per GeneID detected."
)

message(
"These rows will be merged so that annotations from all transcripts",
" are retained."
)

}

###############################################################################

# 6B. KEEP GENES IN THE EXACT H4 UNIVERSE

###############################################################################

annotations_h4 <- annotations_raw %>%

filter(
GeneID %in% gene_universe
)

message(
"Annotation rows belonging to H4 genes: ",
nrow(annotations_h4)
)

message(
"Unique H4 GeneIDs represented in annotation file: ",
n_distinct(
annotations_h4$GeneID
)
)

###############################################################################

# 6C. COLLAPSE TRANSCRIPT-LEVEL ANNOTATIONS TO GENE LEVEL

###############################################################################

#

# A gene may have multiple transcript rows.

#

# For enrichment we need one row per GeneID.

#

# GO, PFAM and InterPro annotations from ALL transcript rows are combined.

#

# Existing separators are retained and a semicolon is inserted between

# transcript-level annotation strings. The downstream parsers use regular

# expressions to extract unique identifiers, so duplicated identifiers are

# removed automatically.

#

collapse_annotations <- function(x) {

x <- as.character(x)

x <- x[
!is.na(x) &
trimws(x) != ""
]

if (length(x) == 0) {

return(
  NA_character_
)

}

paste(
unique(
trimws(x)
),
collapse = ";"
)

}

annotations_gene <- annotations_h4 %>%

dplyr::group_by(
GeneID
) %>%

dplyr::summarise(

`GO Terms` =
  collapse_annotations(
    `GO Terms`
  ),

PFAM =
  collapse_annotations(
    PFAM
  ),

InterPro =
  collapse_annotations(
    InterPro
  ),

.groups = "drop"

)

message(
"Genes after transcript-level annotation aggregation: ",
nrow(annotations_gene)
)

if (
anyDuplicated(
annotations_gene$GeneID
)
) {

stop(
"ERROR: duplicated GeneID values remain after annotation aggregation."
)

}


###############################################################################

# 6D. MAKE SURE EVERY H4 GENE IS REPRESENTED

###############################################################################

#

# The exact H4 universe contains 16,872 genes.

#

# Some H4 genes may have no annotation.

# They must nevertheless remain in the statistical universe.

#

# Therefore we explicitly construct a gene-level annotation table containing

# ALL H4 genes and attach annotations where available.

#

annotations <- data.frame(
GeneID = gene_universe,
stringsAsFactors = FALSE
) %>%

dplyr::left_join(
annotations_gene,
by = "GeneID"
)

message(
"Final gene-level annotation table: ",
nrow(annotations),
" genes"
)

if (
nrow(annotations) !=
length(gene_universe)
) {

stop(
"ERROR: final annotation table does not contain exactly the H4 universe."
)

}

if (
anyDuplicated(
annotations$GeneID
)
) {

stop(
"ERROR: duplicated GeneID values remain after annotation aggregation."
)

}

message(
"Verified: one annotation row per H4 GeneID."
)

###############################################################################

# 6E. SAVE A DIAGNOSTIC TABLE OF ANNOTATION ROW COUNTS

###############################################################################

annotation_row_diagnostics <- annotations_raw %>%

dplyr::group_by(
GeneID
) %>%

dplyr::summarise(
N_annotation_rows = dplyr::n(),
.groups = "drop"
) %>%

dplyr::arrange(
dplyr::desc(N_annotation_rows)
)

write.csv(
annotation_row_diagnostics,
file.path(
output_root,
"annotation_row_diagnostics.csv"
),
row.names = FALSE
)

message(
"Annotation row diagnostics saved."
)

###############################################################################
# 7. ANNOTATION COVERAGE
###############################################################################

message("")
message("============================================================")
message("Calculating annotation coverage")
message("============================================================")


has_annotation <- function(x) {

  !is.na(x) &
    trimws(x) != ""

}


coverage <- data.frame(

  Annotation = c(
    "GO",
    "InterPro",
    "PFAM"
  ),

  N_annotated = c(

    sum(
      has_annotation(
        annotations$`GO Terms`
      )
    ),

    sum(
      has_annotation(
        annotations$InterPro
      )
    ),

    sum(
      has_annotation(
        annotations$PFAM
      )
    )

  ),

  stringsAsFactors = FALSE

) %>%

  mutate(

    N_H4_genes = length(
      gene_universe
    ),

    Proportion_annotated =
      N_annotated / N_H4_genes

  )


write.csv(
  coverage,
  file.path(
    output_root,
    "annotation_coverage.csv"
  ),
  row.names = FALSE
)


print(coverage)


###############################################################################
# 8. HELPER FUNCTIONS
###############################################################################

# ---------------------------------------------------------------------------
# Parse GO IDs
# ---------------------------------------------------------------------------

parse_go <- function(x) {

  if (
    is.na(x) ||
    trimws(x) == ""
  ) {

    return(
      character(0)
    )

  }

  unique(
    str_extract_all(
      x,
      "GO:\\d{7}"
    )[[1]]
  )

}


# ---------------------------------------------------------------------------
# Parse InterPro IDs
# ---------------------------------------------------------------------------

parse_interpro <- function(x) {

  if (
    is.na(x) ||
    trimws(x) == ""
  ) {

    return(
      character(0)
    )

  }

  unique(
    str_extract_all(
      x,
      "IPR\\d{6}"
    )[[1]]
  )

}


# ---------------------------------------------------------------------------
# Parse PFAM IDs
# ---------------------------------------------------------------------------

parse_pfam <- function(x) {

  if (
    is.na(x) ||
    trimws(x) == ""
  ) {

    return(
      character(0)
    )

  }

  unique(
    str_extract_all(
      x,
      "PF\\d{5}"
    )[[1]]
  )

}


# ---------------------------------------------------------------------------
# Build gene -> annotation mapping
# ---------------------------------------------------------------------------

build_mapping <- function(
    annotation_df,
    annotation_column,
    parser
) {

  mapping <- list()

  for (i in seq_len(nrow(annotation_df))) {

    gene <- annotation_df$GeneID[i]

    terms <- parser(
      annotation_df[[annotation_column]][i]
    )

    if (length(terms) == 0) {
      next
    }

    for (term in terms) {

      mapping[[term]] <- c(
        mapping[[term]],
        gene
      )

    }

  }


  mapping <- lapply(
    mapping,
    unique
  )


  mapping

}


# ---------------------------------------------------------------------------
# Filter GO mapping by ontology
# ---------------------------------------------------------------------------
#
# GO ontology namespaces:
#
#   BP = Biological Process
#   MF = Molecular Function
#   CC = Cellular Component
#

filter_go_ontology <- function(
    go_mapping,
    ontology_code
) {

  if (length(go_mapping) == 0) {

    return(
      go_mapping
    )

  }


  go_ids <- names(
    go_mapping
  )


  ontology <- AnnotationDbi::select(
    GO.db,
    keys = go_ids,
    columns = "ONTOLOGY",
    keytype = "GOID"
  ) %>%
    dplyr::distinct(
      GOID,
      .keep_all = TRUE
    )


  selected_ids <- ontology %>%
    filter(
      ONTOLOGY == ontology_code
    ) %>%
    pull(
      GOID
    )


  go_mapping[
    intersect(
      go_ids,
      selected_ids
    )
  ]

}

# ---------------------------------------------------------------------------
# Obtain GO term names
# ---------------------------------------------------------------------------

get_go_names <- function(go_ids) {

  if (length(go_ids) == 0) {

    return(
      data.frame(
        GOID = character(0),
        Term = character(0),
        stringsAsFactors = FALSE
      )
    )

  }


  result <- AnnotationDbi::select(
    GO.db,
    keys = go_ids,
    columns = "TERM",
    keytype = "GOID"
  ) %>%
    dplyr::distinct(
      GOID,
      .keep_all = TRUE
    )


  result %>%
    dplyr::select(
      GOID,
      TERM
    ) %>%
    dplyr::rename(
      Term = TERM
    )

}



# ---------------------------------------------------------------------------
# Retrieve InterPro / PFAM entry descriptions
# ---------------------------------------------------------------------------
#
# InterPro provides programmatic access to both InterPro and Pfam entries.
# The longer database text is stored in Detailed_Description.
#
# Metadata are cached in the output directory to make reruns reproducible and
# to avoid repeated API requests.

# Convert API fields that may be scalar, vectors, or nested lists into one
# safe character string. The InterPro API can return description fields with
# more than one element, so testing is.na(x) before collapsing can cause
# "length > 1 in coercion to logical(1)" errors.
scalar_text <- function(x) {

  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  if (is.list(x)) {
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }

  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = " ")
}


strip_html <- function(x) {

  x <- scalar_text(x)

  if (is.na(x)) {
    return(NA_character_)
  }

  x <- gsub("<[^>]+>", " ", x)
  x <- gsub("\\[[^]]+\\]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}


fetch_domain_metadata_one <- function(accession, source_db) {

  api_url <- paste0(
    "https://www.ebi.ac.uk/interpro/api/entry/",
    source_db,
    "/",
    accession,
    "/"
  )

  for (attempt in seq_len(3)) {

    ans <- tryCatch(
      jsonlite::fromJSON(
        api_url,
        simplifyVector = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(ans) && !is.null(ans$metadata)) {

      meta <- ans$metadata

      entry_name <- NA_character_

      if (!is.null(meta$name)) {
        if (is.list(meta$name)) {
          if (!is.null(meta$name$name)) {
            entry_name <- scalar_text(meta$name$name)
          } else if (!is.null(meta$name$short)) {
            entry_name <- scalar_text(meta$name$short)
          }
        } else {
          entry_name <- scalar_text(meta$name)
        }
      }

      detailed_description <- strip_html(meta$description)

      integrated <- scalar_text(meta$integrated)

      return(
        data.frame(
          Term_ID = accession,
          Description = entry_name,
          Detailed_Description = detailed_description,
          InterPro_Integrated = integrated,
          Metadata_Status = "retrieved",
          Database_URL = paste0(
            "https://www.ebi.ac.uk/interpro/entry/",
            source_db,
            "/",
            accession,
            "/"
          ),
          stringsAsFactors = FALSE
        )
      )
    }

    Sys.sleep(1)
  }

  data.frame(
    Term_ID = accession,
    Description = accession,
    Detailed_Description = NA_character_,
    InterPro_Integrated = NA_character_,
    Metadata_Status = "not_retrieved",
    Database_URL = paste0(
      "https://www.ebi.ac.uk/interpro/entry/",
      source_db,
      "/",
      accession,
      "/"
    ),
    stringsAsFactors = FALSE
  )
}


get_domain_metadata <- function(accessions, source_db, cache_file) {

  accessions <- unique(
    accessions[
      !is.na(accessions) &
      trimws(accessions) != ""
    ]
  )

  if (length(accessions) == 0) {
    return(
      data.frame(
        Term_ID = character(0),
        Description = character(0),
        Detailed_Description = character(0),
        InterPro_Integrated = character(0),
        Metadata_Status = character(0),
        Database_URL = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  if (file.exists(cache_file)) {

    cache <- tryCatch(
      read.csv(
        cache_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )

  } else {
    cache <- NULL
  }

  if (is.null(cache)) {

    cache <- data.frame(
      Term_ID = character(0),
      Description = character(0),
      Detailed_Description = character(0),
      InterPro_Integrated = character(0),
      Metadata_Status = character(0),
      Database_URL = character(0),
      stringsAsFactors = FALSE
    )

  }

  missing_ids <- setdiff(
    accessions,
    cache$Term_ID
  )

  if (length(missing_ids) > 0) {

    new_metadata <- vector(
      "list",
      length(missing_ids)
    )

    for (i in seq_along(missing_ids)) {

      message(
        "  Retrieving ",
        source_db,
        " metadata: ",
        missing_ids[i],
        " (",
        i,
        "/",
        length(missing_ids),
        ")"
      )

      new_metadata[[i]] <- fetch_domain_metadata_one(
        missing_ids[i],
        source_db
      )

      Sys.sleep(API_SLEEP_SECONDS)
    }

    cache <- dplyr::bind_rows(
      cache,
      dplyr::bind_rows(new_metadata)
    ) %>%
      dplyr::distinct(
        Term_ID,
        .keep_all = TRUE
      )

    write.csv(
      cache,
      cache_file,
      row.names = FALSE
    )
  }

  cache %>%
    dplyr::filter(
      Term_ID %in% accessions
    )
}


# ---------------------------------------------------------------------------
# Fisher enrichment for one target gene list
# ---------------------------------------------------------------------------

run_enrichment <- function(
    target_genes,
    universe_genes,
    annotation_mapping,
    annotation_name,
    target_name
) {

  target_genes <- intersect(
    unique(target_genes),
    universe_genes
  )

  universe_genes <- unique(universe_genes)

  if (length(target_genes) == 0) {
    return(data.frame())
  }

  # Calculate category sizes within the analysis background.
  category_sizes <- sapply(
    annotation_mapping,
    function(x) {
      length(intersect(x, universe_genes))
    }
  )

  # Only categories sufficiently represented in the analysis background.
  valid_categories <- names(
    category_sizes[category_sizes >= MIN_CATEGORY_SIZE]
  )

  if (length(valid_categories) == 0) {
    return(data.frame())
  }

  results <- vector("list", length(valid_categories))

  for (i in seq_along(valid_categories)) {

    category <- valid_categories[i]

    annotated_genes <- intersect(
      annotation_mapping[[category]],
      universe_genes
    )

    overlap_genes <- intersect(
      target_genes,
      annotated_genes
    )

    n_overlap <- length(overlap_genes)

    # 2 x 2 contingency table:
    #
    #                     In category   Not in category
    # Target genes             a              b
    # Background only          c              d
    #
    # Fisher's exact test is one-sided for over-representation.

    a <- n_overlap
    b <- length(target_genes) - a
    c <- length(annotated_genes) - a
    d <- length(universe_genes) - length(target_genes) - c

    contingency <- matrix(
      c(a, b, c, d),
      nrow = 2,
      byrow = TRUE
    )

    fisher_result <- fisher.test(
      contingency,
      alternative = "greater"
    )

    results[[i]] <- data.frame(
      Annotation = annotation_name,
      Target = target_name,
      Term_ID = category,
      N_Target = length(target_genes),
      N_Annotated_Background = length(annotated_genes),
      N_Overlap = n_overlap,
      Fold_Enrichment =
        (
          n_overlap / length(target_genes)
        ) /
        (
          length(annotated_genes) / length(universe_genes)
        ),
      Odds_Ratio = unname(fisher_result$estimate),
      P_value = fisher_result$p.value,
      stringsAsFactors = FALSE
    )
  }

  results <- dplyr::bind_rows(results)

  if (nrow(results) == 0) {
    return(results)
  }

  # IMPORTANT:
  # Multiple-testing correction is applied to ALL categories that were tested.
  # The minimum-overlap threshold is a reporting filter and is therefore
  # applied only AFTER BH correction. This avoids making the FDR calculation
  # conditional on observed overlap.
  results <- results %>%
    dplyr::mutate(
      FDR = p.adjust(
        P_value,
        method = "BH"
      )
    ) %>%
    dplyr::filter(
      N_Overlap >= MIN_TARGET_OVERLAP
    ) %>%
    dplyr::arrange(
      FDR,
      P_value
    )

  results
}


###############################################################################
# 9. BUILD ANNOTATION MAPPINGS
###############################################################################

message("")
message("============================================================")
message("Building annotation mappings")
message("============================================================")


go_mapping_all <- build_mapping(
  annotations,
  "GO Terms",
  parse_go
)


go_bp_mapping <- filter_go_ontology(
  go_mapping_all,
  "BP"
)


go_mf_mapping <- filter_go_ontology(
  go_mapping_all,
  "MF"
)


go_cc_mapping <- filter_go_ontology(
  go_mapping_all,
  "CC"
)


interpro_mapping <- build_mapping(
  annotations,
  "InterPro",
  parse_interpro
)


pfam_mapping <- build_mapping(
  annotations,
  "PFAM",
  parse_pfam
)


annotation_mappings <- list(

  GO_BP = go_bp_mapping,

  GO_MF = go_mf_mapping,

  GO_CC = go_cc_mapping,

  InterPro = interpro_mapping,

  PFAM = pfam_mapping

)


message(
  "GO-BP categories: ",
  length(go_bp_mapping)
)

message(
  "GO-MF categories: ",
  length(go_mf_mapping)
)

message(
  "GO-CC categories: ",
  length(go_cc_mapping)
)

message(
  "InterPro categories: ",
  length(interpro_mapping)
)

message(
  "PFAM categories: ",
  length(pfam_mapping)
)


###############################################################################
# 10. GO TERM NAMES
###############################################################################

go_names <- get_go_names(
  names(go_mapping_all)
)


###############################################################################
# 10B. INTERPRO / PFAM METADATA CACHE
###############################################################################

message("")
message("============================================================")
message("Preparing InterPro / PFAM metadata cache")
message("============================================================")

interpro_metadata_file <- file.path(
  output_root,
  "InterPro_metadata_cache.csv"
)

pfam_metadata_file <- file.path(
  output_root,
  "PFAM_metadata_cache.csv"
)

###############################################################################
# 11. RUN ENRICHMENT
###############################################################################

message("")
message("============================================================")
message("Running enrichment analyses")
message("============================================================")


all_results <- list()


for (collection in names(annotation_mappings)) {


  message("")
  message(
    "Annotation collection: ",
    collection
  )


  mapping <- annotation_mappings[[collection]]


  for (target in target_classes) {


    message(
      "  Target: ",
      target
    )


    target_genes <- h4 %>%

      dplyr::filter(
        Category == target
      ) %>%

      dplyr::pull(
        GeneID
      )


    # Use an annotation-aware background by default.
    # This avoids treating missing annotation as biological absence.
    if (BACKGROUND_MODE == "annotated_by_collection") {

      collection_universe <- unique(
        unlist(
          lapply(
            mapping,
            function(x) intersect(x, gene_universe)
          ),
          use.names = FALSE
        )
      )

    } else if (BACKGROUND_MODE == "all") {

      collection_universe <- gene_universe

    } else {

      stop(
        "Unknown BACKGROUND_MODE: ",
        BACKGROUND_MODE
      )

    }

    result <- run_enrichment(

      target_genes =
        target_genes,

      universe_genes =
        collection_universe,

      annotation_mapping =
        mapping,

      annotation_name =
        collection,

      target_name =
        target

    )


    if (nrow(result) == 0) {

      warning(
        "No enrichment categories were retained for ",
        collection,
        " / ",
        target
      )

      next

    }


    # Add GO term names.

    if (collection %in% c("GO_BP", "GO_MF", "GO_CC")) {

      result <- result %>%

        dplyr::left_join(
          go_names,
          by = c(
            "Term_ID" = "GOID"
          )
        ) %>%

        dplyr::mutate(
          Description = Term
        ) %>%

        dplyr::select(
          -Term
        )

    } else if (collection == "InterPro") {

      interpro_metadata <- get_domain_metadata(
        result$Term_ID,
        "interpro",
        interpro_metadata_file
      )

      result <- result %>%
        dplyr::left_join(
          interpro_metadata,
          by = "Term_ID"
        ) %>%
        dplyr::mutate(
          Description = dplyr::if_else(
            is.na(Description) | Description == "",
            Term_ID,
            Description
          )
        )

    } else if (collection == "PFAM") {

      pfam_metadata <- get_domain_metadata(
        result$Term_ID,
        "pfam",
        pfam_metadata_file
      )

      result <- result %>%
        dplyr::left_join(
          pfam_metadata,
          by = "Term_ID"
        ) %>%
        dplyr::mutate(
          Description = dplyr::if_else(
            is.na(Description) | Description == "",
            Term_ID,
            Description
          )
        )

    }


    # Save complete result table.

    collection_dir <- file.path(
      output_root,
      collection
    )


    output_file <- file.path(

      collection_dir,

      paste0(
        collection,
        "_",
        target,
        ".csv"
      )

    )


    write.csv(
      result,
      output_file,
      row.names = FALSE
    )


    # Save significant results.

    significant <- result %>%

      dplyr::filter(
        FDR < FDR_CUTOFF
      ) %>%

      dplyr::arrange(
        FDR,
        P_value
      )


    significant_file <- file.path(

      collection_dir,

      paste0(
        collection,
        "_",
        target,
        "_significant.csv"
      )

    )


    write.csv(
      significant,
      significant_file,
      row.names = FALSE
    )


    all_results[[
      paste(
        collection,
        target,
        sep = "__"
      )
    ]] <- result


    message(
      "    Categories tested: ",
      nrow(result)
    )

    message(
      "    Significant at FDR < ",
      FDR_CUTOFF,
      ": ",
      nrow(significant)
    )

  }

}


###############################################################################
# 12. COMBINED SUMMARY
###############################################################################

message("")
message("============================================================")
message("Creating combined enrichment summary")
message("============================================================")


combined_results <- bind_rows(
  all_results
)


# Save concise domain-name lookup tables used in the results.
if (file.exists(interpro_metadata_file)) {
  write.csv(
    read.csv(
      interpro_metadata_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    file.path(
      output_root,
      "InterPro_term_descriptions.csv"
    ),
    row.names = FALSE
  )
}

if (file.exists(pfam_metadata_file)) {
  write.csv(
    read.csv(
      pfam_metadata_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    file.path(
      output_root,
      "PFAM_term_descriptions.csv"
    ),
    row.names = FALSE
  )
}


write.csv(
  combined_results,
  file.path(
    output_root,
    "H4_all_enrichment_results.csv"
  ),
  row.names = FALSE
)


# Summary of significant terms.

summary_table <- combined_results %>%

  group_by(
    Annotation,
    Target
  ) %>%

  summarise(

    N_categories_tested = n(),

    N_significant_FDR05 =
      sum(
        FDR < FDR_CUTOFF,
        na.rm = TRUE
      ),

    Minimum_FDR =
      min(
        FDR,
        na.rm = TRUE
      ),

    Maximum_fold_enrichment_among_significant =
      ifelse(
        any(
          FDR < FDR_CUTOFF,
          na.rm = TRUE
        ),
        max(
          Fold_Enrichment[
            FDR < FDR_CUTOFF
          ],
          na.rm = TRUE
        ),
        NA_real_
      ),

    .groups = "drop"

  )


write.csv(
  summary_table,
  file.path(
    output_root,
    "H4_enrichment_summary.csv"
  ),
  row.names = FALSE
)


print(
  summary_table
)


###############################################################################
# 13. FINAL MESSAGE
###############################################################################

message("")
message("============================================================")
message("ENRICHMENT ANALYSIS COMPLETED")
message("============================================================")