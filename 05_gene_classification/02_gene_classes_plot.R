#!/usr/bin/env Rscript

###############################################################################
# Author: Iris Sammarco
# Date: 06/10/2026
#
# Run with:
# qsub Fig_Anastatica_gene_classes_run.sh
#
# Final Figure 5
#
# Panel A: Number of genes in each transcriptomic response class
#
# Panel B: Model-attributed variance fractions, shown as TWO PAIRS of boxes:
#
#   Climate panel:
#       Climate-only       -> shown
#       Climate-and-range  -> shown
#       Range-only         -> NOT displayed
#       Non-associated     -> NOT displayed
#
#   Distribution-range panel:
#       Range-only         -> shown
#       Climate-and-range  -> shown
#       Climate-only       -> NOT displayed
#       Non-associated     -> NOT displayed
#
# This is a plotting-only script. It uses:
#   gene_classification.csv
#   gene_level_variance_partition.csv
#
# Outputs (in figures_final/):
#   Figure5_final.png / .pdf
#   Figure5_panelA_gene_counts.png / .pdf
#   Figure5_panelB_variance.png / .pdf
#   Figure5_panelB_plot_data.csv
###############################################################################

options(stringsAsFactors = FALSE)

###############################################################################
# 0. Libraries
###############################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(patchwork)
})

###############################################################################
# 1. Configuration
###############################################################################

BASE_DIR <- "/Anastatica"
H4_DIR   <- file.path(BASE_DIR, "climate_range_variance_analysis")
TABLE_DIR <- file.path(H4_DIR, "tables")
FIG_DIR <- file.path(H4_DIR, "figures_final")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. Input files
###############################################################################

CLASSIFICATION_FILE <- file.path(TABLE_DIR, "gene_classification.csv")
VARIANCE_FILE       <- file.path(TABLE_DIR, "gene_level_variance_partition.csv")

###############################################################################
# 3. Check input files
###############################################################################

if (!file.exists(CLASSIFICATION_FILE)) {
    stop("\nClassification file not found:\n", CLASSIFICATION_FILE,
         "\n\nPlease run the original gene classes analysis first.")
}

if (!file.exists(VARIANCE_FILE)) {
    stop("\nVariance-partition file not found:\n", VARIANCE_FILE,
         "\n\nPlease run the original gene classes analysis first.")
}

message("")
message("Input files found.")
message("")

###############################################################################
# 4. Read input data
###############################################################################

gene_class <- read_csv(CLASSIFICATION_FILE, show_col_types = FALSE)
vp         <- read_csv(VARIANCE_FILE, show_col_types = FALSE)

###############################################################################
# 5. Check required columns
###############################################################################

missing_classification <- setdiff(c("GeneID", "Category"), colnames(gene_class))

if (length(missing_classification) > 0) {
    stop("\nThe classification file is missing required columns:\n",
         paste(missing_classification, collapse = "\n"))
}

missing_variance <- setdiff(c("GeneID", "V_climate", "V_range"), colnames(vp))

if (length(missing_variance) > 0) {
    stop("\nThe variance file is missing required columns:\n",
         paste(missing_variance, collapse = "\n"))
}

###############################################################################
# 6. Class labels and ordering
###############################################################################

class_levels <- c("climate_only", "range_only", "both", "neither")

class_labels <- c(
    "climate_only" = "Climate-only",
    "range_only"   = "Range-only",
    "both"         = "Climate-and-range",
    "neither"      = "Non-associated"
)

###############################################################################
# 7. Consistent colours (same class = same colour in every panel)
###############################################################################

class_colours <- c(
    "Climate-only"      = "#41ab5d", # green
    "Range-only"        = "#E67E22", # orange
    "Climate-and-range" = "#7A5195", # purple
    "Non-associated"    = "#BDBDBD"  # grey
)

###############################################################################
# 8. Prepare classification data
###############################################################################

gene_class_plot <- gene_class %>%
    mutate(
        Category = factor(Category, levels = class_levels),
        Category_label = recode(as.character(Category), !!!class_labels),
        Category_label = factor(Category_label, levels = unname(class_labels))
    )

###############################################################################
# 9. Panel A – raw number of genes
###############################################################################

class_counts <- gene_class_plot %>%
    count(Category_label, name = "N_genes") %>%
    complete(
        Category_label = factor(unname(class_labels), levels = unname(class_labels)),
        fill = list(N_genes = 0)
    )

message("------------------------------------------------------------")
message("Panel A – gene classification counts")
message("------------------------------------------------------------")
print(class_counts)
message("")

###############################################################################
# 10. Prepare variance data
###############################################################################

plot_data <- vp %>%
    select(GeneID, V_climate, V_range) %>%
    left_join(
        gene_class_plot %>% select(GeneID, Category, Category_label),
        by = "GeneID"
    )

if (any(is.na(plot_data$Category))) {
    stop("\n", sum(is.na(plot_data$Category)),
         " genes in the variance table could not be matched ",
         "to the classification table.")
}

###############################################################################
# 11. Panel B data – climate (Climate-only + Climate-and-range)
###############################################################################

climate_plot <- plot_data %>%
    filter(Category %in% c("climate_only", "both")) %>%
    transmute(
        GeneID, Category, Category_label,
        Component = "Climate",
        VarianceFraction = V_climate
    ) %>%
    filter(is.finite(VarianceFraction))

###############################################################################
# 12. Panel B data – range (Range-only + Climate-and-range)
###############################################################################

range_plot <- plot_data %>%
    filter(Category %in% c("range_only", "both")) %>%
    transmute(
        GeneID, Category, Category_label,
        Component = "Distribution range",
        VarianceFraction = V_range
    ) %>%
    filter(is.finite(VarianceFraction))

###############################################################################
# 13. Combine Panel B data
###############################################################################

panel_b_data <- bind_rows(climate_plot, range_plot)

panel_b_data$Component <- factor(
    panel_b_data$Component,
    levels = c("Climate", "Distribution range")
)

panel_b_data$Category_label <- factor(
    panel_b_data$Category_label,
    levels = unname(class_labels)
)

write_csv(panel_b_data, file.path(FIG_DIR, "Figure5_panelB_plot_data.csv"))

###############################################################################
# 14. Diagnostics for Panel B
###############################################################################

message("------------------------------------------------------------")
message("Panel B – genes displayed")
message("------------------------------------------------------------")
message("Climate panel:")
message("  Climate-only: ",      sum(climate_plot$Category == "climate_only"))
message("  Climate-and-range: ", sum(climate_plot$Category == "both"))
message("Distribution-range panel:")
message("  Range-only: ",        sum(range_plot$Category == "range_only"))
message("  Climate-and-range: ", sum(range_plot$Category == "both"))
message("")

###############################################################################
# 15. Panel A – gene counts
###############################################################################

p_a <- ggplot(
    class_counts,
    aes(x = Category_label, y = N_genes, fill = Category_label)
) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
    geom_text(
        aes(label = format(N_genes, big.mark = ",")),
        vjust = -0.35, size = 4
    ) +
    scale_fill_manual(values = class_colours, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Number of genes") +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(8, 8, 8, 8)
    )

###############################################################################
# 16. Panel B – two pairs of two boxes
###############################################################################

# facet_wrap(nrow = 1, scales = "free_x") -> each facet only contains the two
# classes that are actually displayed, so there are no empty x-axis slots.
# The y-axis stays shared (fixed) so the two facets are directly comparable.

p_b <- ggplot(
    panel_b_data,
    aes(x = Category_label, y = VarianceFraction, fill = Category_label)
) +
    geom_boxplot(
        width = 0.6,
        outlier.shape = 16,
        outlier.size = 1.2,
        outlier.alpha = 0.25,
        colour = "grey25"
    ) +
    facet_wrap(~ Component, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = class_colours, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, y = "Model-attributed variance fraction") +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(8, 8, 8, 8)
    )

###############################################################################
# 17. Combine Panel A and Panel B
###############################################################################

final_figure <- p_a +
    p_b +
    plot_layout(widths = c(0.8, 1.4)) +
    plot_annotation(
        tag_levels = "a",
        theme = theme(plot.tag = element_text(face = "bold", size = 14))
    )

###############################################################################
# 18. Save combined figure
###############################################################################

ggsave(
    file.path(FIG_DIR, "Figure5_final.png"),
    final_figure, width = 12, height = 6, units = "in", dpi = 300, bg = "white"
)

ggsave(
    file.path(FIG_DIR, "Figure5_final.pdf"),
    final_figure, width = 12, height = 6, units = "in",
    device = cairo_pdf, bg = "white"
)

###############################################################################
# 19. Save individual panels
###############################################################################

ggsave(
    file.path(FIG_DIR, "Figure5_panelA_gene_counts.png"),
    p_a, width = 5.5, height = 6, units = "in", dpi = 300, bg = "white"
)

ggsave(
    file.path(FIG_DIR, "Figure5_panelA_gene_counts.pdf"),
    p_a, width = 5.5, height = 6, units = "in", device = cairo_pdf, bg = "white"
)

ggsave(
    file.path(FIG_DIR, "Figure5_panelB_variance.png"),
    p_b, width = 7.5, height = 6, units = "in", dpi = 300, bg = "white"
)

ggsave(
    file.path(FIG_DIR, "Figure5_panelB_variance.pdf"),
    p_b, width = 7.5, height = 6, units = "in", device = cairo_pdf, bg = "white"
)

###############################################################################
# 20. Display figure
###############################################################################

print(final_figure)

###############################################################################
# 21. Final message
###############################################################################

message("")
message("============================================================")
message("Figure 5 completed successfully")
message("============================================================")
message("Output directory:\n", FIG_DIR)
message("")
message("Figure structure:")
message("  (a) Raw number of genes in each response class")
message("  (b) Model-attributed variance fractions (two pairs of boxes)")
message("      Climate:            Climate-only + Climate-and-range")
message("      Distribution range: Range-only + Climate-and-range")
message("")
message("Non-associated genes are not displayed in Panel B.")
message("============================================================")
