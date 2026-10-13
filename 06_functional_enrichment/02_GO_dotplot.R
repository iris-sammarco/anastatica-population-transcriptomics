#!/usr/bin/env Rscript

# ============================================================
#
# Author: Iris Sammarco
# Date: 01/10/2026
#
# GO enrichment dot plot
#
# Includes:
#   - Range-only: top 6 GO terms per ontology by FDR
#   - Climate-only: significant GO terms
#   - Climate-and-range (both): significant GO terms
#
# Three panels:
#   A. Biological Process (BP)
#   B. Cellular Component (CC)
#   C. Molecular Function (MF)
#
# Visual encoding:
#   - X-axis        = Fold enrichment
#   - Point size    = Number of overlapping genes
#   - Point colour  = -log10(FDR)
#   - Point shape   = Gene set
#
# Gene-set categories:
#   - Range-only
#   - Climate-only
#   - Climate-and-range
#
# Outputs:
#   - PNG: 600 dpi
#   - PDF: vector graphics
#
# Run with:
# mamba activate R
# Rscript 02_GO_dotplot.R
#
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(tidyverse)
library(scales)


# ------------------------------------------------------------
# 2. Input / output paths
# ------------------------------------------------------------

# Range-only:
# Top 6 GO terms per ontology by FDR
range_file <- "Anastatica/H4_climate_range_variance_analysis/GO_pathway_analysis/top6_GO_terms_byFDR_range_only.csv"

# Climate-only:
# Significant GO terms
climate_file <- "/Anastatica/H4_climate_range_variance_analysis/GO_pathway_analysis/GO_MF/GO_MF_climate_only_significant.csv"

# Climate-and-range:
# Significant GO terms
both_file <- "/Anastatica/H4_climate_range_variance_analysis/GO_pathway_analysis/GO_CC/GO_CC_both_significant.csv"


# Output directory
output_dir <- "/Anastatica/H4_climate_range_variance_analysis/GO_pathway_analysis"


# Output filenames
png_file <- file.path(
  output_dir,
  "GO_terms_dotplot.png"
)

pdf_file <- file.path(
  output_dir,
  "GO_terms_dotplot.pdf"
)


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------

# Range-only
go_range <- read.csv(
  range_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Climate-only
go_climate <- read.csv(
  climate_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Climate-and-range
go_both <- read.csv(
  both_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Add gene-set labels
# ------------------------------------------------------------

go_range <- go_range %>%
  mutate(
    Gene_set = "Range-only"
  )


go_climate <- go_climate %>%
  mutate(
    Gene_set = "Climate-only"
  )


go_both <- go_both %>%
  mutate(
    Gene_set = "Climate-and-range"
  )


# ------------------------------------------------------------
# 5. Combine all GO results
# ------------------------------------------------------------

go <- bind_rows(
  go_range,
  go_climate,
  go_both
)


# ------------------------------------------------------------
# 6. Clean / prepare data
# ------------------------------------------------------------

go <- go %>%
  mutate(

    # Ontology
    Annotation = factor(
      Annotation,
      levels = c(
        "GO_BP",
        "GO_CC",
        "GO_MF"
      )
    ),

    Ontology = recode(
      Annotation,
      GO_BP = "Biological Process",
      GO_CC = "Cellular Component",
      GO_MF = "Molecular Function"
    ),

    # Make sure numeric columns are numeric
    Fold_Enrichment = as.numeric(Fold_Enrichment),
    N_Overlap = as.numeric(N_Overlap),
    FDR = as.numeric(FDR),
    P_value = as.numeric(P_value),

    # Statistical significance
    neg_log10_FDR = -log10(FDR),

    # Gene-set ordering
    Gene_set = factor(
      Gene_set,
      levels = c(
        "Range-only",
        "Climate-only",
        "Climate-and-range"
      )
    )
  )


# ------------------------------------------------------------
# 7. Wrap long GO-term descriptions
# ------------------------------------------------------------
#
# This prevents very long GO descriptions from making the
# figure excessively wide.
#
# The original Description column is retained.
# Description_label is used only for plotting.
# ------------------------------------------------------------

go <- go %>%
  mutate(
    Description_label = str_wrap(
      Description,
      width = 50
    )
  )


# ------------------------------------------------------------
# 8. Order GO terms within each ontology
#
# Terms are ordered from most significant (lowest FDR)
# to least significant (highest FDR).
#
# The ordering is performed separately within each ontology.
# ------------------------------------------------------------

go <- go %>%
  group_by(Annotation) %>%
  arrange(FDR, .by_group = TRUE) %>%
  mutate(
    Description_ordered = factor(
      Description_label,
      levels = rev(unique(Description_label))
    )
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 9. Set ontology order
# ------------------------------------------------------------

go$Ontology <- factor(
  go$Ontology,
  levels = c(
    "Biological Process",
    "Cellular Component",
    "Molecular Function"
  )
)


# ------------------------------------------------------------
# 10. Set shape values
# ------------------------------------------------------------
#
# Shapes:
#
#   16 = filled circle      = Range-only
#   17 = filled triangle    = Climate-only
#   15 = filled square      = Climate-and-range
#
# Using shape rather than another colour scale keeps colour
# available for statistical significance (FDR).
# ------------------------------------------------------------

gene_set_shapes <- c(
  "Range-only" = 16,
  "Climate-only" = 17,
  "Climate-and-range" = 15
)


# ------------------------------------------------------------
# 11. Plot
# ------------------------------------------------------------

p <- ggplot(
  go,
  aes(
    x = Fold_Enrichment,
    y = Description_ordered
  )
) +

  # ----------------------------------------------------------
  # Main points
  # ----------------------------------------------------------
  #
  # Size    = number of overlapping genes
  # Colour  = -log10(FDR)
  # Shape   = gene set
  #
  geom_point(
    aes(
      size = N_Overlap,
      colour = neg_log10_FDR,
      shape = Gene_set
    ),
    alpha = 0.9
  ) +

  # ----------------------------------------------------------
  # One panel for each GO ontology
  # ----------------------------------------------------------

  facet_grid(
    rows = vars(Ontology),
    scales = "free_y",
    space = "free_y"
  ) +

  # ----------------------------------------------------------
  # Colour = statistical significance
  # ----------------------------------------------------------

  scale_colour_viridis_c(
    option = "plasma",
    name = expression(-log[10]~"(FDR)"),
    begin = 0.15,
    end = 0.95
  ) +

  # ----------------------------------------------------------
  # Size = number of overlapping genes
  # ----------------------------------------------------------

  scale_size_continuous(
    name = "Number of\noverlapping genes",
    range = c(3.5, 9),
    breaks = pretty(
      go$N_Overlap,
      n = 4
    )
  ) +

  # ----------------------------------------------------------
  # Shape = gene set
  # ----------------------------------------------------------

  scale_shape_manual(
    name = "Gene set",
    values = gene_set_shapes
  ) +
  
  guides(
  shape = guide_legend(
    title.theme = element_text(
      size = 13,
      face = "bold"
    ),
    label.theme = element_text(
      size = 12
    ),
    keyheight = unit(0.8, "cm"),
    keywidth = unit(0.8, "cm")
    )
  ) +

  # ----------------------------------------------------------
  # X-axis
  # ----------------------------------------------------------

  scale_x_continuous(
    name = "Fold enrichment",
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +

  # ----------------------------------------------------------
  # Y-axis
  # ----------------------------------------------------------

  scale_y_discrete(
    name = NULL
  ) +

  # ----------------------------------------------------------
  # Publication-style theme
  # ----------------------------------------------------------

  theme_classic(
    base_size = 11
  ) +

  theme(

    # --------------------------------------------------------
    # Overall title
    # --------------------------------------------------------

    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 14
      )
    ),

    # --------------------------------------------------------
    # Facet titles
    # --------------------------------------------------------

    strip.background = element_blank(),

    strip.text = element_text(
      size = 13,
      face = "bold"
    ),

    strip.placement = "outside",

    # --------------------------------------------------------
    # GO term labels
    # --------------------------------------------------------

    axis.text.y = element_text(
      size = 10.5,
      colour = "black",
      lineheight = 0.95
    ),

    # --------------------------------------------------------
    # X-axis tick labels
    # --------------------------------------------------------

    axis.text.x = element_text(
      size = 11,
      colour = "black"
    ),

    # --------------------------------------------------------
    # X-axis title
    # --------------------------------------------------------

    axis.title.x = element_text(
      size = 12,
      margin = margin(
        t = 10
      )
    ),

    # --------------------------------------------------------
    # Grid
    # --------------------------------------------------------

    panel.grid.major.x = element_line(
      linewidth = 0.25,
      colour = "grey85"
    ),

    panel.grid.major.y = element_blank(),

    panel.grid.minor = element_blank(),

    # --------------------------------------------------------
    # Space between BP / CC / MF panels
    # --------------------------------------------------------

    panel.spacing = unit(
      1.2,
      "lines"
    ),

    # --------------------------------------------------------
    # Legend titles
    # --------------------------------------------------------

    legend.title = element_text(
      size = 12,
      face = "bold"
    ),

    # --------------------------------------------------------
    # Legend labels
    # --------------------------------------------------------

    legend.text = element_text(
      size = 11
    ),

    # --------------------------------------------------------
    # Legend position
    # --------------------------------------------------------

    legend.position = "right",

    # --------------------------------------------------------
    # Figure margins
    # --------------------------------------------------------

    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 15
    )
  )


# ------------------------------------------------------------
# 12. Display plot
# ------------------------------------------------------------

print(p)


# ------------------------------------------------------------
# 13. Save high-resolution PNG
# ------------------------------------------------------------

ggsave(
  filename = png_file,
  plot = p,
  width = 13,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 14. Save vector PDF
# ------------------------------------------------------------

ggsave(
  filename = pdf_file,
  plot = p,
  width = 13,
  height = 8.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# ------------------------------------------------------------
# 15. Confirmation message
# ------------------------------------------------------------

cat("\n")
cat("GO dot plot successfully generated.\n")
cat("\n")