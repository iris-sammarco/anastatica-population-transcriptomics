#!/usr/bin/env Rscript

# Date: 28/08/2026
#
# Author: Iris Sammarco
#
# Purpose:
#   This script reads sample metadata from a CSV file, detects Salmon
#   quantification outputs in sample subdirectories, matches metadata
#   sample names to Salmon directories, builds a tx2gene table directly
#   from a GFF3 annotation, checks transcript ID compatibility between
#   Salmon quant.sf and the annotation, and imports Salmon quantifications
#   directly at transcript and gene level using tximeta.
#
# Input:
#   - CSV metadata file with column:
#       Sample_name
#   - Salmon quantification outputs (quant.sf files) in directories named like:
#       IS17-6_S34/quant.sf
#   - Custom transcriptome FASTA file
#   - Custom annotation GFF3 file
#   - Salmon index directory
#
# Output:
#   Metadata directory:
#     - Anastatica_coldata_for_tximeta.csv
#
#   tximeta directory:
#     - Anastatica_linked_txome.json
#     - Anastatica_tx2gene.csv
#     - tximeta_transcript_level.rds
#     - tximeta_gene_level.rds
#     - tximeta_gene_counts.csv
#     - tximeta_import_summary.txt
#
# Usage:
#   mamba activate tximeta_env
#   Rscript 04_make_metadata_and_run_tximeta.r

###############################################################################
# tximeta / BiocFileCache setup
###############################################################################

# Clear tximeta cache
unlink("/scratch/R/BiocFileCache", recursive = TRUE)

# Set and create tximeta cache directory automatically
tximeta_cache <- "/scratch/R/BiocFileCache"

if (!dir.exists(tximeta_cache)) {
  dir.create(tximeta_cache, recursive = TRUE, showWarnings = FALSE)
}

Sys.setenv(TXIMETA_HUB_CACHE = tximeta_cache)

###############################################################################
# Libraries
###############################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tximeta)
  library(SummarizedExperiment)
  library(rtracklayer)
})

tximeta::setTximetaBFC(tximeta_cache, quiet = TRUE)

message("tximeta BiocFileCache directory: ", tximeta::getTximetaBFC())

###############################################################################
# Configuration
###############################################################################

input_metadata_csv <- "/Anastatica/metadata/Anastatica_RNA_metadata_tximeta.csv"
quant_dir <- "/Anastatica/salmon_output/quant"
index_dir <- "/Anastatica/salmon_output/salmon_index"
genome_dir <- "/Anastatica/genome_annotation"
transcripts_fa <- file.path(genome_dir, "Anastatica_hierochuntica_ANAS.mrna-transcripts.fa")
gff3_file <- file.path(genome_dir, "Anastatica_hierochuntica_ANAS.gff3")

# Metadata output directory
outdir <- "/Anastatica/metadata"

# tximeta output directory
tximeta_outdir <- "/Anastatica/tximeta"

metadata_csv <- file.path(outdir, "Anastatica_metadata.csv")
coldata_csv <- file.path(outdir, "Anastatica_coldata_for_tximeta.csv")

transcript_rds <- file.path(tximeta_outdir, "tximeta_transcript_level.rds")
gene_rds <- file.path(tximeta_outdir, "tximeta_gene_level.rds")
gene_counts_csv <- file.path(tximeta_outdir, "tximeta_gene_counts.csv")
tx2gene_csv <- file.path(tximeta_outdir, "Anastatica_tx2gene.csv")
import_summary_txt <- file.path(tximeta_outdir, "tximeta_import_summary.txt")

# Kept for record only; main import below does not rely on linkedTxome
linked_json <- file.path(tximeta_outdir, "Anastatica_linked_txome.json")

###############################################################################
# Input checks
###############################################################################

stopifnot(file.exists(input_metadata_csv))
stopifnot(dir.exists(quant_dir))
stopifnot(dir.exists(index_dir))
stopifnot(file.exists(transcripts_fa))
stopifnot(file.exists(gff3_file))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tximeta_outdir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Read metadata CSV file
###############################################################################

meta <- read_csv(input_metadata_csv, show_col_types = FALSE)

if (!"Sample_name" %in% colnames(meta)) {
  stop("Metadata CSV must contain a column named 'Sample_name'")
}

meta <- meta %>%
  mutate(Sample_name = str_trim(as.character(Sample_name)))

if (anyDuplicated(meta$Sample_name)) {
  dup <- unique(meta$Sample_name[duplicated(meta$Sample_name)])
  stop("Duplicate Sample_name values found in metadata: ", paste(dup, collapse = ", "))
}

###############################################################################
# Detect Salmon quantification directories
###############################################################################

quant_dirs <- list.dirs(quant_dir, recursive = FALSE, full.names = TRUE)

if (length(quant_dirs) == 0) {
  stop("No Salmon sample directories found in: ", quant_dir)
}

quant_df <- tibble(
  quant_dir_name = basename(quant_dirs),
  Sample_name = str_remove(basename(quant_dirs), "_S\\d+$"),
  files = file.path(quant_dirs, "quant.sf")
) %>%
  filter(file.exists(files))

if (nrow(quant_df) == 0) {
  stop("No quant.sf files found under: ", quant_dir)
}

if (anyDuplicated(quant_df$quant_dir_name)) {
  dup <- unique(quant_df$quant_dir_name[duplicated(quant_df$quant_dir_name)])
  stop("Duplicate quantification directory names detected: ", paste(dup, collapse = ", "))
}

if (anyDuplicated(quant_df$Sample_name)) {
  dup <- unique(quant_df$Sample_name[duplicated(quant_df$Sample_name)])
  stop(
    "Multiple Salmon directories map to the same Sample_name after removing _S<number>: ",
    paste(dup, collapse = ", ")
  )
}

###############################################################################
# Join metadata with detected quant files
###############################################################################

coldata <- meta %>%
  left_join(quant_df, by = "Sample_name") %>%
  arrange(quant_dir_name)

missing_files <- coldata %>%
  filter(is.na(files)) %>%
  pull(Sample_name)

if (length(missing_files) > 0) {
  warning(
    "These metadata samples have no matching Salmon quant.sf file and will be excluded: ",
    paste(missing_files, collapse = ", ")
  )

  coldata <- coldata %>%
    filter(!Sample_name %in% missing_files)
}

if (nrow(coldata) == 0) {
  stop("After excluding samples with missing quant.sf files, coldata is empty.")
}

###############################################################################
# Prepare and write metadata outputs
###############################################################################

write.csv(as.data.frame(coldata), metadata_csv, row.names = FALSE, quote = TRUE)

coldata <- coldata %>%
  mutate(
    names = quant_dir_name,
    Population = if ("Population" %in% colnames(.)) factor(Population) else NULL,
    Country = if ("Country" %in% colnames(.)) factor(Country) else NULL,
    Country_code = if ("Country_code" %in% colnames(.)) factor(Country_code) else NULL,
    Site_of_distribution_range = if ("Site_of_distribution_range" %in% colnames(.)) factor(Site_of_distribution_range) else NULL,
	Site_of_distribution_range_East_West = if ("Site_of_distribution_range_East_West" %in% colnames(.))
    factor(Site_of_distribution_range_East_West) else NULL,	
    Sequencing_batch = if ("Sequencing_batch" %in% colnames(.)) factor(Sequencing_batch) else NULL,
    Genetic_cluster = if ("Genetic_cluster" %in% colnames(.)) factor(Genetic_cluster) else NULL
  ) %>%
  select(names, files, everything()) %>%
  as.data.frame(stringsAsFactors = FALSE)

stopifnot(all(c("names", "files") %in% colnames(coldata)))

write.csv(coldata, coldata_csv, row.names = FALSE, quote = TRUE)

###############################################################################
# Build tx2gene table from GFF3
###############################################################################

message("Importing GFF3 annotation: ", gff3_file)

gff <- rtracklayer::import(gff3_file)
gff_df <- as.data.frame(gff)

message("GFF3 feature types detected:")
print(table(gff_df$type))

tx2gene <- gff_df %>%
  filter(type %in% c("mRNA", "transcript")) %>%
  transmute(
    TXNAME = as.character(ID),
    GENEID = as.character(Parent)
  ) %>%
  filter(!is.na(TXNAME), !is.na(GENEID)) %>%
  distinct()

if (nrow(tx2gene) == 0) {
  stop(
    "tx2gene table is empty. ",
    "The GFF3 file may not use ID/Parent on mRNA/transcript features. ",
    "Please inspect the GFF3 attributes."
  )
}

write.csv(tx2gene, tx2gene_csv, row.names = FALSE, quote = TRUE)

message("tx2gene rows: ", nrow(tx2gene))
message("Unique transcripts in tx2gene: ", length(unique(tx2gene$TXNAME)))
message("Unique genes in tx2gene: ", length(unique(tx2gene$GENEID)))

###############################################################################
# Check compatibility between quant.sf transcript names and tx2gene
###############################################################################

first_quant <- coldata$files[1]

quant_header <- read.delim(first_quant, nrows = 5, stringsAsFactors = FALSE)

if (!"Name" %in% colnames(quant_header)) {
  stop("The first quant.sf file does not contain a 'Name' column: ", first_quant)
}

quant_tx <- quant_header$Name

message("Example transcript names from first quant.sf:")
print(head(quant_tx))

message("Example transcript names from tx2gene:")
print(head(tx2gene$TXNAME))

first_quant_all <- read.delim(first_quant, stringsAsFactors = FALSE)
missing_from_tx2gene <- setdiff(first_quant_all$Name, tx2gene$TXNAME)

message("Transcripts in first quant.sf: ", length(unique(first_quant_all$Name)))
message("Transcripts from first quant.sf missing in tx2gene: ",
        length(unique(missing_from_tx2gene)))

if (length(unique(missing_from_tx2gene)) > 0) {
  message("First few transcript IDs missing from tx2gene:")
  print(head(unique(missing_from_tx2gene), 20))
  stop(
    "Transcript IDs in quant.sf do not fully match TXNAME in tx2gene. ",
    "You probably need to adjust transcript IDs parsed from the GFF3 file."
  )
}

###############################################################################
# Optional: create linkedTxome JSON for record
###############################################################################
# If the Salmon index lacks digest metadata, this may fail. The main import
# below does not depend on linkedTxome and will continue with tx2gene.

tryCatch({
  makeLinkedTxome(
    indexDir = index_dir,
    source = "Custom",
    organism = "Anastatica hierochuntica",
    release = "ANAS_custom",
    genome = "ANAS",
    fasta = transcripts_fa,
    gtf = gff3_file,
    jsonFile = linked_json,
    write = TRUE
  )
}, error = function(e) {
  warning("makeLinkedTxome failed, but continuing with tx2gene import. Error: ",
          conditionMessage(e))
})

###############################################################################
# Import Salmon quantifications at transcript level
###############################################################################
# Saved for diagnostic/reference purposes only.

message("Importing Salmon quantifications at transcript level...")

se <- tximeta(
  coldata,
  skipMeta = TRUE
)

saveRDS(se, transcript_rds)

###############################################################################
# Import Salmon quantifications directly at gene level
###############################################################################

message("Importing Salmon quantifications directly at gene level using tx2gene...")

gse <- tximeta(
  coldata,
  tx2gene = tx2gene,
  txOut = FALSE,
  skipMeta = TRUE
)

saveRDS(gse, gene_rds)

write.csv(
  as.data.frame(assay(gse, "counts")),
  gene_counts_csv,
  row.names = TRUE
)

###############################################################################
# Write summary
###############################################################################

summary_lines <- c(
  paste("Metadata CSV:", input_metadata_csv),
  paste("Samples imported:", ncol(gse)),
  paste("Transcript rows:", nrow(se)),
  paste("Gene rows:", nrow(gse)),
  paste("Metadata CSV output:", metadata_csv),
  paste("Coldata CSV:", coldata_csv),
  paste("tx2gene CSV:", tx2gene_csv),
  paste("Linked txome JSON:", linked_json),
  paste("Transcript-level RDS:", transcript_rds),
  paste("Gene-level RDS:", gene_rds),
  paste("Gene counts CSV:", gene_counts_csv),
  paste("tximeta cache:", tximeta::getTximetaBFC())
)

writeLines(summary_lines, import_summary_txt)

message("Done. Wrote metadata to ", outdir,
        " and tximeta outputs to ", tximeta_outdir)