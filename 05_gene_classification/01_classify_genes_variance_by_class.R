#!/usr/bin/env Rscript

###############################################################################
# Author: Iris Sammarco
# Date: 10/09/2026
#
# Climate versus distribution-range contributions to transcriptomic differentiation
#
# Primary model:
#
# expression ~ ECe_mean + Climate_PC1 + Climate_PC2 + Climate_PC3 +
#              Site_of_distribution_range +
#              (1 | Sequencing_batch) +
#              (1 | Population)
#
# IMPORTANT:
# ---------------------------------------------------------------------------
# The supplied dream/eBayes fit is used to define gene-level climate/range
# association. ECe_mean (soil salinity) is modeled separately and is not
# included in the climate classification.
#
# Gene-wise variance components are NOT reconstructed from coefficient
# magnitudes. Instead, variancePartition is refit using the same expression
# matrix, samples, genes, metadata, and model structure.
#
# Climate is defined a priori as:
#   Climate_PC1 + Climate_PC2 + Climate_PC3
#
# ECe_mean is treated separately as a salinity variable.
#
# Range is:
#   Site_of_distribution_range
#
# The four response classes are defined independently from the variance
# partitioning:
#   climate_only
#   range_only
#   both
#   neither
#
# Multiple testing:
# ---------------------------------------------------------------------------
# BH correction is performed across all three climate tests + all three range
# contrasts for each gene (6 tests/gene), providing a single family for the
# classification.
#
# Variance partition:
# ---------------------------------------------------------------------------
# V_climate = Climate_PC1 + Climate_PC2 + Climate_PC3
# V_ECe     = ECe_mean
# V_range   = Site_of_distribution_range
#
# Other components are retained separately:
#   Population
#   Sequencing_batch
#   Residuals
#
# Outputs:
# ---------------------------------------------------------------------------
# - Gene classifications (climate_only, range_only, both, neither): gene_classification.csv
# - gene-level variance fractions
# - summary statistics by class
# - statistical comparisons between classes
# - bootstrap confidence intervals
# - diagnostic and summary tables
# - diagnostic tables
#
###############################################################################

options(stringsAsFactors = FALSE)

###############################################################################
# 0. Libraries
###############################################################################

suppressPackageStartupMessages({

    library(SummarizedExperiment)
    library(edgeR)
    library(limma)
    library(variancePartition)
    library(lme4)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
})

###############################################################################
# 1. Configuration
###############################################################################

BASE_DIR <- "/Anastatica"

MODEL_FILE <- file.path(
    BASE_DIR,
    "degs",
    "ECe_climate_range_model",
    "ECe_climate_range_model_fit_contrast_ebayes.rds"
)

# Expression object used in the original analysis.
TXIMETA_FILE <- file.path(
    BASE_DIR,
    "tximeta",
    "tximeta_gene_level.rds"
)

# New independent output directory.
OUT_DIR <- file.path(
    BASE_DIR,
    "climate_range_variance_analysis"
)

TABLE_DIR <- file.path(OUT_DIR, "tables")
LOG_DIR   <- file.path(OUT_DIR, "logs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Statistical settings
###############################################################################

FDR_CUTOFF <- 0.05

# Minimum expression filter.
# This is deliberately based on the scientific fixed effects rather than
# sequencing batch alone.
MIN_COUNT <- 10

# Number of bootstrap replicates.
# Increase to 5000 for final manuscript analysis if computational resources
# permit.
N_BOOT <- 2000

set.seed(20260909)

###############################################################################
# 3. Model definition
###############################################################################

CLIMATE_FIXED <- c(
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3"
)

SALINITY_FIXED <- "ECe_mean"

RANGE_FIXED <- "Site_of_distribution_range"

RANGE_CONTRASTS <- c(
    "East_vs_Center",
    "West_vs_Center",
    "East_vs_West"
)

MODEL_FORMULA <- ~
    ECe_mean +
    Climate_PC1 +
    Climate_PC2 +
    Climate_PC3 +
    (1 | Site_of_distribution_range) +
    (1 | Sequencing_batch) +
    (1 | Population)

###############################################################################
# 4. Helper functions
###############################################################################

message_log <- function(...) {
    txt <- paste0(...)
    message(txt)
    cat(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        " ",
        txt,
        "\n",
        file = file.path(LOG_DIR, "analysis.log"),
        append = TRUE
    )
}

stop_if_missing <- function(x, label) {
    if (!file.exists(x)) {
        stop(label, " not found:\n", x)
    }
}

extract_gene_id <- function(feature_id) {

    feature_id <- as.character(feature_id)

    # Handles identifiers such as:
    # transcript.GENE
    # gene.GENE
    # GENE
    #
    # The fitted model is expected to use the same IDs as the expression
    # object. We therefore first try exact matching later.

    ifelse(
        grepl("\\.", feature_id),
        sub("^[^.]+\\.", "", feature_id),
        feature_id
    )
}

###############################################################################
# 5. Check input files
###############################################################################

message_log("Starting analysis.")
message_log("R version: ", R.version.string)

stop_if_missing(
    MODEL_FILE,
    "Saved dream/eBayes model"
)

stop_if_missing(
    TXIMETA_FILE,
    "Gene-level tximeta object"
)

message_log("Model: ", MODEL_FILE)
message_log("Expression object: ", TXIMETA_FILE)
message_log("Output directory: ", OUT_DIR)

###############################################################################
# 6. Load dream/eBayes model
###############################################################################

fit <- readRDS(MODEL_FILE)

if (!inherits(fit, "MArrayLM")) {

    stop(
        "The supplied RDS is not an MArrayLM object.\n",
        "Class: ",
        paste(class(fit), collapse = ", ")
    )
}

if (is.null(fit$coef)) {
    stop("The saved model contains no coefficient matrix.")
}

fit_genes <- rownames(fit$coef)

if (is.null(fit_genes) || length(fit_genes) == 0) {
    stop("No gene identifiers found in the saved model.")
}

message_log(
    "Genes in saved model: ",
    length(fit_genes)
)

message_log(
    "Model coefficients:\n",
    paste(colnames(fit$coef), collapse = "\n")
)

###############################################################################
# 7. Verify required coefficients
###############################################################################

# The supplied fit may have only the coefficients used for contrasts.
# We nevertheless require the three climate coefficients and the range
# coefficients needed for the gene classification.

required_coefficients <- c(
    CLIMATE_FIXED,
    RANGE_CONTRASTS
)

missing_coefficients <- setdiff(
    required_coefficients,
    colnames(fit$coef)
)

if (length(missing_coefficients) > 0) {

    stop(
        "\nThe supplied model does not contain all coefficients required.\n\n",
        "Missing:\n",
        paste(missing_coefficients, collapse = "\n"),
        "\n\nAvailable coefficients:\n",
        paste(colnames(fit$coef), collapse = "\n"),
        "\n\nThis is important: variance partitioning can still be performed, ",
        "but the independent DE-based four-class classification cannot be ",
        "reconstructed from this fit."
    )
}

###############################################################################
# 8. Load exact expression object
###############################################################################

gse <- readRDS(TXIMETA_FILE)

if (!inherits(gse, "SummarizedExperiment")) {
    stop("The tximeta RDS is not a SummarizedExperiment.")
}

if (!"counts" %in% assayNames(gse)) {
    stop(
        "The tximeta object does not contain an assay called 'counts'."
    )
}

counts <- assay(gse, "counts")
metadata <- as.data.frame(colData(gse))

message_log(
    "Expression matrix: ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " samples"
)

###############################################################################
# 9. Match samples
###############################################################################

if (!"quant_dir_name" %in% colnames(metadata)) {

    stop(
        "Metadata does not contain 'quant_dir_name', which was used ",
        "to align samples in the original analysis."
    )
}

rownames(metadata) <- metadata$quant_dir_name

if (!all(colnames(counts) %in% rownames(metadata))) {

    stop(
        "Not all count-matrix samples are present in metadata."
    )
}

metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
]

if (!identical(colnames(counts), rownames(metadata))) {

    stop(
        "Counts and metadata are not in identical sample order."
    )
}

###############################################################################
# 10. Check metadata required by the model
###############################################################################

required_metadata <- c(
    "ECe_mean",
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3",
    "Site_of_distribution_range",
    "Sequencing_batch",
    "Population"
)

missing_metadata <- setdiff(
    required_metadata,
    colnames(metadata)
)

if (length(missing_metadata) > 0) {

    stop(
        "Required metadata columns are missing:\n",
        paste(missing_metadata, collapse = "\n")
    )
}

metadata$Population <- factor(metadata$Population)
metadata$Sequencing_batch <- factor(metadata$Sequencing_batch)
metadata$Site_of_distribution_range <-
    factor(metadata$Site_of_distribution_range)

for (x in c(
    "ECe_mean",
    "Climate_PC1",
    "Climate_PC2",
    "Climate_PC3"
)) {
    metadata[[x]] <- as.numeric(metadata[[x]])
}

if (anyNA(metadata[, required_metadata])) {

    stop(
        "NA values detected in variables used by the model."
    )
}

###############################################################################
# 11. Verify range reference level
###############################################################################

if ("Center" %in% levels(metadata$Site_of_distribution_range)) {

    metadata$Site_of_distribution_range <-
        relevel(
            metadata$Site_of_distribution_range,
            ref = "Center"
        )

} else {

    warning(
        "Center is not a level of Site_of_distribution_range. ",
        "Proceeding with existing reference level."
    )
}

###############################################################################
# 12. Reproduce expression filtering
###############################################################################

# The original analysis used filterByExpr(). We use the complete scientific
# model here rather than filtering only against sequencing batch.

design_filter <- model.matrix(
    ~
        ECe_mean +
        Climate_PC1 +
        Climate_PC2 +
        Climate_PC3 +
        Site_of_distribution_range,
    data = metadata
)

keep <- filterByExpr(
    counts,
    design = design_filter
)

filtered_genes <- rownames(counts)[keep]

message_log(
    "Genes after expression filtering: ",
    length(filtered_genes)
)

###############################################################################
# 13. Restrict to the exact saved-fit universe
###############################################################################

common_genes <- intersect(
    fit_genes,
    filtered_genes
)

message_log(
    "Genes common to saved fit and filtered expression object: ",
    length(common_genes)
)

if (length(common_genes) < 100) {

    stop(
        "Fewer than 100 genes could be matched between the saved model ",
        "and expression object."
    )
}

# For scientific reproducibility, do NOT silently use a different gene
# universe. Require the model and filtering to agree.

if (!setequal(fit_genes, filtered_genes)) {

    only_fit <- setdiff(fit_genes, filtered_genes)
    only_expression <- setdiff(filtered_genes, fit_genes)

    write.csv(
        data.frame(
            GeneID = only_fit
        ),
        file.path(
            TABLE_DIR,
            "diagnostic_genes_only_in_saved_model.csv"
        ),
        row.names = FALSE
    )

    write.csv(
        data.frame(
            GeneID = only_expression
        ),
        file.path(
            TABLE_DIR,
            "diagnostic_genes_only_in_expression_filter.csv"
        ),
        row.names = FALSE
    )

    warning(
        "Saved model gene universe and independently reproduced ",
        "filtering are not identical. The analysis will use their ",
        "intersection, and the discrepancy is reported in the diagnostic files."
    )
}

analysis_genes <- intersect(
    fit_genes,
    filtered_genes
)

counts_h4 <- counts[
    analysis_genes,
    ,
    drop = FALSE
]

fit_h4 <- fit[
    match(analysis_genes, rownames(fit$coef)),
]

###############################################################################
# 14. Build DE table from the saved fit
###############################################################################

message_log(
    "Constructing gene-level association statistics."
)

de_list <- list()

for (coef_name in required_coefficients) {

    beta <- fit_h4$coef[, coef_name]

    # Prefer t/p information from the saved fit when available.
    if (!is.null(fit_h4$t)) {
        statistic <- fit_h4$t[, coef_name]
    } else {
        statistic <- beta / fit_h4$stdev.unscaled[, coef_name]
    }

    if (!is.null(fit_h4$p.value)) {
        pval <- fit_h4$p.value[, coef_name]
    } else {
        pval <- 2 * pnorm(-abs(statistic))
    }

    de_list[[coef_name]] <- data.frame(
        GeneID = rownames(fit_h4$coef),
        Test = coef_name,
        logFC = beta,
        statistic = statistic,
        P.Value = pval,
        stringsAsFactors = FALSE
    )
}

de_long <- bind_rows(de_list)

###############################################################################
# 15. Global BH correction across all tests
###############################################################################

# This is a single family consisting of:
#   4 climate effects + 3 range contrasts
#
# for the complete gene universe.

de_long$FDR_global <- p.adjust(
    de_long$P.Value,
    method = "BH"
)

de_long$Significant <- (
    is.finite(de_long$FDR_global) &
    de_long$FDR_global <= FDR_CUTOFF
)

###############################################################################
# 16. Classify genes
###############################################################################

climate_genes <- de_long %>%
    filter(
        Test %in% CLIMATE_FIXED,
        Significant
    ) %>%
    distinct(GeneID) %>%
    pull(GeneID)

range_genes <- de_long %>%
    filter(
        Test %in% RANGE_CONTRASTS,
        Significant
    ) %>%
    distinct(GeneID) %>%
    pull(GeneID)

gene_class <- data.frame(
    GeneID = analysis_genes,
    stringsAsFactors = FALSE
) %>%
    mutate(
        Climate_associated = GeneID %in% climate_genes,
        Range_associated = GeneID %in% range_genes,
        Category = case_when(
            Climate_associated & Range_associated ~ "both",
            Climate_associated & !Range_associated ~ "climate_only",
            !Climate_associated & Range_associated ~ "range_only",
            TRUE ~ "neither"
        )
    )

gene_class$Category <- factor(
    gene_class$Category,
    levels = c(
        "climate_only",
        "range_only",
        "both",
        "neither"
    )
)

###############################################################################
# 17. Save gene classification
###############################################################################

write.csv(
    gene_class,
    file.path(
        TABLE_DIR,
        "gene_classification.csv"
    ),
    row.names = FALSE
)

write.csv(
    de_long,
    file.path(
        TABLE_DIR,
        "gene_level_statistics_global_BH.csv"
    ),
    row.names = FALSE
)

class_counts <- gene_class %>%
    count(Category, name = "N_genes") %>%
    mutate(
        Proportion = N_genes / sum(N_genes)
    )

# Write summary of Category, N_genes, Proportion. This corresponds to Fig 1
write.csv(
    class_counts,
    file.path(
        TABLE_DIR,
        "gene_class_counts.csv"
    ),
    row.names = FALSE
)

message_log(
    "Gene classes:\n",
    paste(
        class_counts$Category,
        class_counts$N_genes,
        sep = ": ",
        collapse = "\n"
    )
)

###############################################################################
# 18. Variance partitioning
###############################################################################

message_log(
    "Starting variancePartition analysis."
)

# Important:
# Use the same filtered genes and exact sample metadata.
#
# VariancePartition operates on normalized expression data. We therefore
# create log-CPM values using edgeR normalization.

dge <- DGEList(
    counts = counts_h4
)

dge <- calcNormFactors(
    dge,
    method = "TMM"
)

logCPM <- cpm(
    dge,
    log = TRUE,
    prior.count = 2
)

###############################################################################
# 19. Fit variance-partition model
###############################################################################

vp_formula <- MODEL_FORMULA

message_log(
    "VariancePartition formula:\n",
    paste(deparse(vp_formula), collapse = " ")
)

# The variancePartition model is fitted gene-by-gene.
#
# BPPARAM is deliberately left at the default here for portability.
# For an HPC implementation, this can be changed to SnowParam/MulticoreParam.

vp_fit <- fitExtractVarPartModel(
    logCPM,
    vp_formula,
    metadata
)

###############################################################################
# 20. Extract variance fractions
###############################################################################

vp <- as.data.frame(vp_fit)

vp$GeneID <- rownames(vp)

if (!"ECe_mean" %in% colnames(vp)) {

    stop(
        "VariancePartition output does not contain ECe_mean."
    )
}

if (!"Climate_PC1" %in% colnames(vp) ||
    !"Climate_PC2" %in% colnames(vp) ||
    !"Climate_PC3" %in% colnames(vp)) {

    stop(
        "VariancePartition output does not contain all climate PCs."
    )
}

if (!"Site_of_distribution_range" %in% colnames(vp)) {

    stop(
        "VariancePartition output does not contain ",
        "Site_of_distribution_range."
    )
}

###############################################################################
# 21. Calculate scientifically interpretable variance components
###############################################################################

vp_h4 <- vp %>%
    mutate(

        # Climate contribution: climate PCs only.
        V_climate =
            Climate_PC1 +
            Climate_PC2 +
            Climate_PC3,

        # Salinity / edaphic contribution.
        V_ECe =
            ECe_mean,

        # Geographic/range contribution.
        V_range =
            Site_of_distribution_range,

        # Other components.
        V_population =
            if ("Population" %in% colnames(.))
                Population
            else NA_real_,

        V_batch =
            if ("Sequencing_batch" %in% colnames(.))
                Sequencing_batch
            else NA_real_,

        V_residual =
            if ("Residuals" %in% colnames(.))
                Residuals
            else NA_real_,

        # Climate versus range.
        Delta_climate_minus_range =
            V_climate - V_range,

        Climate_range_ratio =
            (V_climate + 1e-8) /
            (V_range + 1e-8),

        log2_Climate_range_ratio =
            log2(Climate_range_ratio)
    )

###############################################################################
# 22. Basic variance diagnostics
###############################################################################

if (any(
    vp_h4$V_climate < -1e-8 |
    vp_h4$V_range < -1e-8,
    na.rm = TRUE
)) {

    stop(
        "Negative climate/range variance fractions detected."
    )
}

vp_h4 <- vp_h4 %>%
    filter(
        GeneID %in% analysis_genes
    ) %>%
    left_join(
        gene_class,
        by = "GeneID"
    )

if (any(is.na(vp_h4$Category))) {

    stop(
        "Some variance-partition genes could not be matched to the gene classes."
    )
}

###############################################################################
# 23. Save complete variance table
###############################################################################

# V_climate: V_Climate_PC1 + V_Climate_PC2 + V_Climate_PC3. How much model-attributed variance is associated with the climate PCs.
# V_ECe: This tells you whether salinity contributes appreciably to expression variance even though you have zero ECe DEGs.
# Delta_climate_minus_range: Positive: climate > range. Negative: range > climate
# log2_Climate_range_ratio: Useful for looking at relative dominance

write.csv(
    vp_h4,
    file.path(
        TABLE_DIR,
        "gene_level_variance_partition.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 24. Variance summaries by class
###############################################################################

variance_summary <- vp_h4 %>%
    pivot_longer(
        cols = c(
            V_climate,
			V_ECe,
            V_range,
            V_population,
            V_batch,
            V_residual
        ),
        names_to = "Component",
        values_to = "VarianceFraction"
    ) %>%
    filter(
        is.finite(VarianceFraction)
    ) %>%
    group_by(
        Category,
        Component
    ) %>%
    summarise(
        N_genes = n(),
        Mean = mean(VarianceFraction),
        Median = median(VarianceFraction),
        SD = sd(VarianceFraction),
        Q25 = quantile(VarianceFraction, 0.25),
        Q75 = quantile(VarianceFraction, 0.75),
        Min = min(VarianceFraction),
        Max = max(VarianceFraction),
        .groups = "drop"
    )

# Summarise each variance component within each gene category.
# Prioritise median and IQR over mean because variance fractions can be highly skewed
write.csv(
    variance_summary,
    file.path(
        TABLE_DIR,
        "variance_partition_summary_by_class.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 25. Primary descriptive summary
###############################################################################

h4_summary <- vp_h4 %>%
    group_by(Category) %>%
    summarise(

        N_genes = n(),

        Mean_climate_variance =
            mean(V_climate, na.rm = TRUE),

        Median_climate_variance =
            median(V_climate, na.rm = TRUE),

        Mean_range_variance =
            mean(V_range, na.rm = TRUE),

        Median_range_variance =
            median(V_range, na.rm = TRUE),

        Mean_delta =
            mean(
                Delta_climate_minus_range,
                na.rm = TRUE
            ),

        Median_delta =
            median(
                Delta_climate_minus_range,
                na.rm = TRUE
            ),

        Mean_log2_CR =
            mean(
                log2_Climate_range_ratio,
                na.rm = TRUE
            ),

        Median_log2_CR =
            median(
                log2_Climate_range_ratio,
                na.rm = TRUE
            ),

        .groups = "drop"
    )

write.csv(
    h4_summary,
    file.path(
        TABLE_DIR,
        "climate_range_summary_by_class.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 26. Statistical comparison of climate and range variance
###############################################################################

# Every gene contributes one V_climate and one V_range.
#
# The paired comparison is therefore:
#
#   V_climate vs V_range
#
# across genes.
#
# This is a genome-wide descriptive/inferential analysis, not a claim that
# genes are independent biological replicates. We therefore also report
# bootstrap CIs below.

paired_test <- suppressWarnings(
    wilcox.test(
        vp_h4$V_climate,
        vp_h4$V_range,
        paired = TRUE,
        exact = FALSE
    )
)

paired_summary <- data.frame(
    N_genes = nrow(vp_h4),
    Wilcoxon_V_climate_vs_V_range_W =
        unname(paired_test$statistic),
    Wilcoxon_P =
        paired_test$p.value,
    Median_V_climate =
        median(vp_h4$V_climate, na.rm = TRUE),
    Median_V_range =
        median(vp_h4$V_range, na.rm = TRUE),
    Median_delta =
        median(
            vp_h4$Delta_climate_minus_range,
            na.rm = TRUE
        ),
    Mean_delta =
        mean(
            vp_h4$Delta_climate_minus_range,
            na.rm = TRUE
        )
)

# Use as supporting statistical evidence, rather than central biological conclusion
write.csv(
    paired_summary,
    file.path(
        TABLE_DIR,
        "paired_climate_vs_range_variance_test.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 27. Bootstrap CI for genome-wide climate-range difference
###############################################################################

message_log(
    "Running bootstrap with ",
    N_BOOT,
    " replicates."
)

delta_values <- vp_h4$Delta_climate_minus_range
delta_values <- delta_values[is.finite(delta_values)]

bootstrap_delta <- replicate(
    N_BOOT,
    mean(
        sample(
            delta_values,
            size = length(delta_values),
            replace = TRUE
        )
    )
)

bootstrap_result <- data.frame(
    Estimate_mean_delta =
        mean(delta_values),
    Bootstrap_CI_2.5 =
        quantile(bootstrap_delta, 0.025),
    Bootstrap_CI_97.5 =
        quantile(bootstrap_delta, 0.975),
    N_bootstrap = N_BOOT
)

# CI = negative → negative, range tends to dominate.
# CI = positive → positive, climate tends to dominate.
# CI crosses 0, there isn't strong evidence for an overall difference in the genome-wide mean.
# This is across genes, not independent biological replicates, so use only as descriptive
write.csv(
    bootstrap_result,
    file.path(
        TABLE_DIR,
        "bootstrap_climate_minus_range_CI.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 28. Bootstrap class-specific climate/range differences
###############################################################################

bootstrap_class <- lapply(
    levels(vp_h4$Category),
    function(class_name) {

        x <- vp_h4 %>%
            filter(Category == class_name)

        delta <- x$Delta_climate_minus_range
        delta <- delta[is.finite(delta)]

        if (length(delta) < 2) {

            return(
                data.frame(
                    Category = class_name,
                    N_genes = length(delta),
                    Mean_delta = NA_real_,
                    CI_2.5 = NA_real_,
                    CI_97.5 = NA_real_
                )
            )
        }

        boot <- replicate(
            N_BOOT,
            mean(
                sample(
                    delta,
                    length(delta),
                    replace = TRUE
                )
            )
        )

        data.frame(
            Category = class_name,
            N_genes = length(delta),
            Mean_delta = mean(delta),
            CI_2.5 = quantile(boot, 0.025),
            CI_97.5 = quantile(boot, 0.975)
        )
    }
)

bootstrap_class <- bind_rows(
    bootstrap_class
)

# Useful for Fig 5
# For each class (climate_only, range_only, both, neither) outputs: Mean_delta, CI_2.5, CI_97.5. It asks:
# Does the variance pattern behave as expected within each DE-defined class?
# But circularity problem, as the classification was already constructed using climate/range association. So, you can present this variancePartition analysis as a complementary quantification of the magnitude of model-attributed variation underlying the DE-based classes.

write.csv(
    bootstrap_class,
    file.path(
        TABLE_DIR,
        "bootstrap_climate_range_difference_by_class.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 29. Save all variance-component summaries
###############################################################################

# Prepare variance-component values for summary tables only.
component_long <- vp_h4 %>%
    select(
        GeneID,
        Category,
        V_climate,
        V_ECe,
        V_range,
        V_population,
        V_batch,
        V_residual
    ) %>%
    pivot_longer(
        cols = c(
            V_climate,
            V_ECe,
            V_range,
            V_population,
            V_batch,
            V_residual
        ),
        names_to = "Component",
        values_to = "VarianceFraction"
    ) %>%
    mutate(
        Component = recode(
            Component,
            V_climate = "Climate",
            V_ECe = "ECe / salinity",
            V_range = "Distribution range",
            V_population = "Population",
            V_batch = "Sequencing batch",
            V_residual = "Residual"
        )
    )

component_summary_fig6 <- component_long %>%
    group_by(
        Category,
        Component
    ) %>%
    summarise(
        Mean = mean(VarianceFraction, na.rm = TRUE),
        Median = median(VarianceFraction, na.rm = TRUE),
        Q25 = quantile(VarianceFraction, 0.25, na.rm = TRUE),
        Q75 = quantile(VarianceFraction, 0.75, na.rm = TRUE),
        .groups = "drop"
    )

write.csv(
    component_summary_fig6,
    file.path(
        TABLE_DIR,
        "Figure6_mean_median_variance_components.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 35. Statistical comparison between gene classes
###############################################################################

# Kruskal-Wallis tests are used because variance fractions are bounded and
# typically non-normal.
#
# These tests are descriptive/secondary: the classes were defined from
# differential association, so this should NOT be interpreted as an
# independent test of the original classification criterion.

kw_climate <- kruskal.test(
    V_climate ~ Category,
    data = vp_h4
)

kw_range <- kruskal.test(
    V_range ~ Category,
    data = vp_h4
)

kw_delta <- kruskal.test(
    Delta_climate_minus_range ~ Category,
    data = vp_h4
)

class_tests <- data.frame(
    Variable = c(
        "V_climate",
        "V_range",
        "Delta_climate_minus_range"
    ),
    Kruskal_Wallis_chisq = c(
        unname(kw_climate$statistic),
        unname(kw_range$statistic),
        unname(kw_delta$statistic)
    ),
    P_value = c(
        kw_climate$p.value,
        kw_range$p.value,
        kw_delta$p.value
    )
)

class_tests$FDR <- p.adjust(
    class_tests$P_value,
    method = "BH"
)

# Three tests: climate variance among classes, range variance among classes, climate-minus-range difference among classes. Use FDR.
# Useful for describing class differences.

write.csv(
    class_tests,
    file.path(
        TABLE_DIR,
        "Kruskal_Wallis_variance_comparisons.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 36. Correlation between climate and range variance
###############################################################################

spearman_test <- suppressWarnings(
    cor.test(
        vp_h4$V_climate,
        vp_h4$V_range,
        method = "spearman",
        exact = FALSE
    )
)

correlation_summary <- data.frame(
    N_genes = nrow(vp_h4),
    Spearman_rho =
        unname(spearman_test$estimate),
    P_value =
        spearman_test$p.value
)

# This tells you whether genes with high climate variance also tend to have high range variance.
# Strong positive rho: Genes that are highly climate-associated also tend to be highly range-associated. That supports some degree of shared transcriptomic differentiation.
# Near-zero rho: Climate and range contributions vary relatively independently among genes.
# Negative rho: Genes with strong climate contributions tend to have weaker range contributions, and vice versa. That could support a "partitioning" model.

write.csv(
    correlation_summary,
    file.path(
        TABLE_DIR,
        "climate_range_variance_correlation.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 37. Identify variance-dominant genes
###############################################################################

vp_h4 <- vp_h4 %>%
    mutate(
        Variance_dominance = case_when(
            V_climate > V_range ~ "Climate",
            V_range > V_climate ~ "Range",
            TRUE ~ "Equal"
        )
    )

dominance_summary <- vp_h4 %>%
    count(
        Variance_dominance,
        name = "N_genes"
    ) %>%
    mutate(
        Proportion = N_genes / sum(N_genes)
    )

# Climate, Range, Equal based on: V_climate > V_range
# Don't interpret "Climate-dominant" as causal climate adaptation. It means:
# the model-attributed climate variance fraction exceeds the model-attributed range variance fraction.

write.csv(
    dominance_summary,
    file.path(
        TABLE_DIR,
        "variance_dominance_summary.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 38. Save analysis metadata/settings
###############################################################################

settings <- data.frame(
    Parameter = c(
        "Model_file",
        "Expression_file",
        "Output_directory",
        "FDR_cutoff",
        "N_bootstrap",
        "Climate_definition",
		"Salinity_definition",
        "Range_definition",
        "Variance_partition_input",
        "Expression_transformation"
    ),
    Value = c(
        MODEL_FILE,
        TXIMETA_FILE,
        OUT_DIR,
        FDR_CUTOFF,
        N_BOOT,
        paste(CLIMATE_FIXED, collapse = " + "),
		SALINITY_FIXED,
        RANGE_FIXED,
        "TMM-normalized log2-CPM",
        "log2 CPM with prior.count = 2"
    )
)

write.csv(
    settings,
    file.path(
        TABLE_DIR,
        "analysis_settings.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 39. Final diagnostic report
###############################################################################

diagnostic <- data.frame(
    Quantity = c(
        "Samples",
        "Genes in saved fit",
        "Genes after independent expression filtering",
        "Genes used",
        "Climate-associated genes",
        "Range-associated genes",
        "Climate-only genes",
        "Range-only genes",
        "Both genes",
        "Neither genes",
        "Mean climate variance",
        "Mean range variance",
        "Median climate variance",
        "Median range variance",
        "Spearman climate-range variance",
        "Paired Wilcoxon P",
        "Bootstrap mean delta",
        "Bootstrap CI lower",
        "Bootstrap CI upper"
    ),
    Value = c(
        ncol(counts_h4),
        length(fit_genes),
        length(filtered_genes),
        length(analysis_genes),
        length(climate_genes),
        length(range_genes),
        sum(gene_class$Category == "climate_only"),
        sum(gene_class$Category == "range_only"),
        sum(gene_class$Category == "both"),
        sum(gene_class$Category == "neither"),
        mean(vp_h4$V_climate, na.rm = TRUE),
        mean(vp_h4$V_range, na.rm = TRUE),
        median(vp_h4$V_climate, na.rm = TRUE),
        median(vp_h4$V_range, na.rm = TRUE),
        unname(spearman_test$estimate),
        paired_test$p.value,
        bootstrap_result$Estimate_mean_delta,
        bootstrap_result$Bootstrap_CI_2.5,
        bootstrap_result$Bootstrap_CI_97.5
    )
)

write.csv(
    diagnostic,
    file.path(
        TABLE_DIR,
        "final_diagnostic_summary.csv"
    ),
    row.names = FALSE
)

###############################################################################
# 40. Session information
###############################################################################

writeLines(
    capture.output(
        sessionInfo()
    ),
    con = file.path(
        LOG_DIR,
        "sessionInfo.txt"
    )
)

###############################################################################
# 41. Finish
###############################################################################

message_log(
    "analysis completed successfully"
)

message_log(
    "Tables: ",
    TABLE_DIR
)

print(
    diagnostic
)