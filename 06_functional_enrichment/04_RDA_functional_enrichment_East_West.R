#!/usr/bin/env Rscript

###############################################################################
# EAST / WEST CLIMATE – GO / FUNCTIONAL ENRICHMENT ANALYSIS
# Anastatica hierochuntica
#
# Date: 05/10/2026
# Author: Iris Sammarco
#
# PURPOSE
# -------
# This script performs functional enrichment analysis of the East/West climate gene classes
# obtained from the differential-expression analysis (interaction model):
#
#   1. East_specific
#   2. West_specific
#   3. shared
#
# The "neither" class is NOT treated as a biological target set. Instead,
# all genes with valid East and West omnibus tests are used as the statistical
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
#
# IMPORTANT STATISTICAL PRINCIPLE
# --------------------------------
# The East/West climate analysis universe contains 16,872 eligible genes. By default,
# enrichment uses the subset of that universe that has at least one annotation
# in the collection being tested (annotation-aware background). This is the
# recommended setting because unannotated genes cannot be classified as members
# of a functional category. Set BACKGROUND_MODE = "all_eligible" to use the complete eligible-gene universe.
#
# This is important because enrichment asks:
#
#   "Are particular functions over-represented among East/West climate-associated genes
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
# East_specific
# West_specific
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
# 1. Gene classification:
# /Anastatica/degs/ECe_climate_range_interaction_east_west_model/
# ECe_climate_range_interaction_east_west_model_East_West_climate_gene_classification.csv
#
# Expected columns:
#   FeatureID
#   West_climate_FDR
#   East_climate_FDR
#   West_climate_associated
#   East_climate_associated
#   Climate_class
#   Model_failed
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
# /Anastatica/degs/ECe_climate_range_interaction_east_west_model/East_West_climate_functional_annotation
#
# Output structure:
#
# GO_pathway_analysis/
#
#   GO_BP/
#       GO_BP_East_specific.csv
#       GO_BP_West_specific.csv
#       GO_BP_shared.csv
#       GO_BP_East_specific_significant.csv
#       GO_BP_West_specific_significant.csv
#       GO_BP_shared_significant.csv
#       GO_BP_East_specific_dotplot.png
#       GO_BP_West_specific_dotplot.png
#       GO_BP_shared_dotplot.png
#       GO_BP_shared_dotplot.pdf
#
#   GO_MF/
#       GO_MF_East_specific.csv
#       GO_MF_West_specific.csv
#       GO_MF_shared.csv
#       GO_MF_East_specific_significant.csv
#       GO_MF_West_specific_significant.csv
#       GO_MF_shared_significant.csv
#       GO_MF_East_specific_dotplot.png
#       GO_MF_West_specific_dotplot.png
#       GO_MF_shared_dotplot.png
#
#   GO_CC/
#       GO_CC_East_specific.csv
#       GO_CC_West_specific.csv
#       GO_CC_shared.csv
#       GO_CC_East_specific_significant.csv
#       GO_CC_West_specific_significant.csv
#       GO_CC_shared_significant.csv
#       GO_CC_East_specific_dotplot.png
#       GO_CC_West_specific_dotplot.png
#       GO_CC_shared_dotplot.png
#
#   InterPro/
#       InterPro_East_specific.csv
#       InterPro_West_specific.csv
#       InterPro_shared.csv
#       ...
#
#   PFAM/
#       PFAM_East_specific.csv
#       PFAM_West_specific.csv
#       PFAM_shared.csv
#       ...
#
#   annotation_coverage.csv
#   InterPro_term_descriptions.csv
#   PFAM_term_descriptions.csv
#   InterPro_metadata_cache.csv
#   PFAM_metadata_cache.csv
#   East_West_climate_enrichment_summary.csv
#   sessionInfo.txt
#   analysis_parameters.txt
#
# Dotplots are saved as both PNG (300 dpi) and vector PDF.
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
  library(ggplot2)
  library(GO.db)
  library(AnnotationDbi)
  library(jsonlite)

})

# Avoid masking problems between generic functions and dplyr specific functions:
select <- dplyr::select
filter  <- dplyr::filter
count   <- dplyr::count
rename  <- dplyr::rename
slice   <- dplyr::slice

###############################################################################
# 2. PARAMETERS
###############################################################################

# ---------------------------------------------------------------------------
# Main East/West climate directory
# ---------------------------------------------------------------------------

ew_root <- "/Anastatica/degs/ECe_climate_range_interaction_east_west_model"


# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------

classification_file <- file.path(
  ew_root,
  "ECe_climate_range_interaction_east_west_model_East_West_climate_gene_classification.csv"
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
  ew_root,
  "East_West_climate_functional_annotation"
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


# Number of top terms displayed in dotplots.

N_TOP_PLOT <- 20L

# Statistical background mode.
# "annotated_by_collection" is recommended for enrichment because genes with
# no annotation in a given collection cannot be classified into its categories.
# "all_eligible" uses the complete East/West climate analysis universe.
BACKGROUND_MODE <- "annotated_by_collection"

# InterPro/PFAM API settings.
# Metadata are cached locally so that subsequent runs do not re-query the API.
API_SLEEP_SECONDS <- 0.2


# Target East/West climate classes.

target_classes <- c(
  "East_specific",
  "West_specific",
  "shared"
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
# 5. LOAD EAST / WEST CLIMATE GENE CLASSIFICATION
###############################################################################

message("")
message("============================================================")
message("Loading East / West climate gene classification")
message("============================================================")

classification_raw <- read.csv(
  classification_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_classification_columns <- c(
  "FeatureID",
  "East_climate_FDR",
  "West_climate_FDR",
  "Climate_class"
)

missing_classification_columns <- setdiff(
  required_classification_columns,
  colnames(classification_raw)
)

if (length(missing_classification_columns) > 0) {
  stop(
    "East/West climate classification file is missing columns:\n",
    paste(
      missing_classification_columns,
      collapse = ", "
    )
  )
}

classification_raw$FeatureID <- as.character(
  classification_raw$FeatureID
)

if (anyDuplicated(classification_raw$FeatureID)) {
  stop(
    "Duplicated FeatureID values detected in East/West climate classification."
  )
}

# Only genes successfully tested in BOTH regions can be classified.
# Genes labelled "not_tested" are excluded from the enrichment universe.
classification <- classification_raw %>%
  filter(
    Climate_class != "not_tested"
  ) %>%
  mutate(
    GeneID = FeatureID,
    Category = case_when(
      Climate_class == "East_only" ~ "East_specific",
      Climate_class == "West_only" ~ "West_specific",
      Climate_class == "shared" ~ "shared",
      Climate_class == "neither" ~ "neither",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Category)
  ) %>%
  select(
    GeneID,
    Category,
    East_climate_FDR,
    West_climate_FDR
  )

gene_universe <- unique(
  classification$GeneID
)

message(
  "Genes with valid East and West omnibus climate tests: ",
  length(gene_universe)
)

message("")
message("East / West climate class counts:")

print(
  classification %>%
    count(Category) %>%
    mutate(
      Proportion = n / sum(n)
    )
)

# Keep the downstream enrichment code unchanged by using this object name.
h4 <- classification


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

# 6B. KEEP GENES IN THE EXACT GENE UNIVERSE

###############################################################################

annotations_h4 <- annotations_raw %>%

filter(
GeneID %in% gene_universe
)

message(
"Annotation rows belonging to East/West climate genes: ",
nrow(annotations_h4)
)

message(
"Unique GeneIDs represented in annotation file: ",
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

# 6D. MAKE SURE EVERY East/West GENE IS REPRESENTED

###############################################################################

#

# The exact East/West climate universe contains 16,872 genes.

#

# Some East/West climate genes may have no annotation.

# They must nevertheless remain in the statistical universe.

#

# Therefore we explicitly construct a gene-level annotation table containing

# ALL East/West climate genes and attach annotations where available.

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
"ERROR: final annotation table does not contain exactly the East/West climate universe."
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
"Verified: one annotation row per GeneID."
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

    N_genes = length(
      gene_universe
    ),

    Proportion_annotated =
      N_annotated / N_genes

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
# The concise entry name is stored in Description for tables and figures.
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

    } else if (BACKGROUND_MODE == "all_eligible") {

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
    "East_West_climate_all_enrichment_results.csv"
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
    "East_West_climate_enrichment_summary.csv"
  ),
  row.names = FALSE
)


print(
  summary_table
)


###############################################################################
# 13. DOTPLOTS OF TOP ENRICHED TERMS
###############################################################################

message("")
message("============================================================")
message("Creating enrichment dotplots")
message("============================================================")


for (collection in names(annotation_mappings)) {


  for (target in target_classes) {


    key <- paste(
      collection,
      target,
      sep = "__"
    )


    if (!key %in% names(all_results)) {
      next
    }


    result <- all_results[[key]]


    plot_data <- result %>%

      filter(
        FDR < FDR_CUTOFF
      ) %>%

      arrange(
        FDR,
        P_value
      ) %>%

      slice_head(
        n = N_TOP_PLOT
      )


    if (nrow(plot_data) == 0) {
      next
    }


    # Use term description for GO and ID for domain collections.

    plot_data <- plot_data %>%

      mutate(

        Label = ifelse(
          is.na(Description) |
            Description == "",
          Term_ID,
          Description
        ),

        Label = str_wrap(
          Label,
          width = 55
        ),

        Label = factor(
          Label,
          levels = rev(
            unique(Label)
          )
        )

      )


    p <- ggplot(
      plot_data,
      aes(
        x = Fold_Enrichment,
        y = Label,
        size = N_Overlap,
        colour = -log10(FDR)
      )
    ) +

      geom_point() +

      labs(

        title =
          paste(
            collection,
            "enrichment:",
            target
          ),

        subtitle =
          paste0(
            "Top enriched categories; FDR < ",
            FDR_CUTOFF
          ),

        x =
          "Fold enrichment",

        y =
          NULL,

        size =
          "Genes in target set",

        colour =
          expression(
            -log[10](FDR)
          )

      ) +

      theme_bw() +

      theme(

        plot.title =
          element_text(
            face = "bold"
          ),

        axis.text.y =
          element_text(
            size = 9
          )

      )


    output_file <- file.path(

      output_root,

      collection,

      paste0(
        collection,
        "_",
        target,
        "_dotplot.png"
      )

    )


    ggsave(

      filename =
        output_file,

      plot =
        p,

      width =
        10,

      height =
        8,

      dpi =
        300

    )

    # Also save a vector PDF.
    output_pdf <- file.path(

      output_root,

      collection,

      paste0(
        collection,
        "_",
        target,
        "_dotplot.pdf"
      )

    )

    ggsave(

      filename =
        output_pdf,

      plot =
        p,

      width =
        10,

      height =
        8,

      device =
        grDevices::pdf

    )

  }

}


###############################################################################
# 14. CREATE A SIMPLE TOP-TERM TABLE
###############################################################################

message("")
message("============================================================")
message("Creating top-term tables")
message("============================================================")


ppt_top_terms <- combined_results %>%

  filter(
    FDR < FDR_CUTOFF
  ) %>%

  group_by(
    Annotation,
    Target
  ) %>%

  arrange(
    FDR,
    P_value
  ) %>%

  slice_head(
    n = 20
  ) %>%

  ungroup()


write.csv(
  ppt_top_terms,
  file.path(
    output_root,
    "East_West_climate_top_enriched_terms_for_PPT.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 15. SAVE PARAMETERS AND SESSION INFORMATION
###############################################################################

writeLines(

  c(

    "============================================================",

    "EAST / WEST CLIMATE GO / FUNCTIONAL ENRICHMENT ANALYSIS",

    "============================================================",

    "",

    "Species:",

    "Anastatica hierochuntica",

    "",

    "East/West climate classification input:",

    classification_file,

    "",

    "Functional annotation input:",

    annotation_file,

    "",

    "East/West climate gene universe:",

    paste0(
      length(gene_universe),
      " genes"
    ),

    "",

    "Target gene classes:",

    paste(
      target_classes,
      collapse = ", "
    ),

    "",

    "Annotation collections:",

    paste(
      annotation_collections,
      collapse = ", "
    ),

    "",

    "GO ontology:",

    "Biological Process (BP), Molecular Function (MF), Cellular Component (CC)",
	
    "",

    "Enrichment test:",

    "One-sided Fisher's exact test (over-representation)",

    "",

    "Multiple-testing correction:",

    "Benjamini-Hochberg FDR within each target and annotation collection",

    "",

    paste0(
      "Minimum annotation category size: ",
      MIN_CATEGORY_SIZE
    ),

    "",

    paste0(
      "Minimum target overlap: ",
      MIN_TARGET_OVERLAP
    ),

    "",

    paste0(
      "Significance threshold: FDR < ",
      FDR_CUTOFF
    ),

    "",

    "Background:",

    paste0(
      "East/West climate universe mode: ",
      BACKGROUND_MODE
    ),

    ifelse(
      BACKGROUND_MODE == "annotated_by_collection",
      "For each annotation collection, only East/West climate genes with at least one annotation in that collection are used as the statistical universe.",
      "All genes in the East/West climate analysis universe are used as the statistical universe."
    ),

    "",

    "Important:",

    "Genes labelled not_tested are excluded from the enrichment universe.",

    "Genes with valid East and West omnibus climate tests, including neither genes, form the enrichment universe.",

    "",

    "KEGG:",

    "Not analysed because the supplied annotation file does not contain",

    "KEGG/KO pathway identifiers.",

    "",

    "============================================================"

  ),

  con =
    file.path(
      output_root,
      "analysis_parameters.txt"
    )

)


writeLines(

  capture.output(
    sessionInfo()
  ),

  con =
    file.path(
      output_root,
      "sessionInfo.txt"
    )

)


###############################################################################
# 16. FINAL MESSAGE
###############################################################################

message("")
message("============================================================")
message("ENRICHMENT ANALYSIS COMPLETED")
message("============================================================")

message("")
message("Output directory:")
message(output_root)

message("")
message("Main combined result:")
message(
  file.path(
    output_root,
    "East_West_climate_all_enrichment_results.csv"
  )
)

message("")
message("PPT-ready top terms:")
message(
  file.path(
    output_root,
    "East_West_climate_top_enriched_terms_for_PPT.csv"
  )
)

message("")
message("Significant enrichment is defined as:")
message(
  "FDR < ",
  FDR_CUTOFF
)

message("")
message("============================================================")