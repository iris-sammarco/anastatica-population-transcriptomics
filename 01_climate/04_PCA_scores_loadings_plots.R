#!/usr/bin/env Rscript

# ============================================================
#
# Author: Iris Sammarco
# Date: 16/09/2026
#
###############################################################################
# Anastatica hierochuntica
# Supplementary Figure 2B, C
#
# Fig. 2B: Climate PCA scores of sampled populations
#           PC1 vs PC2, PC1 vs PC3, PC2 vs PC3
#
# Fig. 2C: PCA loadings for PC1-PC3
#           Top N variables by absolute loading
#
# Outputs:
#   /Anastatica/climate/pca_output/
#
# PNG + PDF
#
# Run with:
# mamba activate R
# Rscript 04_PCA_scores_loadings_plots.R
#
###############################################################################

############################
# 1. Packages
############################

required_packages <- c(
  "tidyverse",
  "ggrepel",
  "patchwork",
  "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)


############################
# 2. Input/output files
############################

metadata_file <- "/Anastatica/metadata/Anastatica_RNA_metadata_tximeta.csv"

loadings_file <- "/Anastatica/climate/pca_output/bioclim_pca_loadings.csv"

output_dir <- "/Anastatica/climate/pca_output/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


############################
# 3. General figure settings
############################

# Number of top variables to display for each PC
TOP_N <- 6

# Figure dimensions
FIG2B_WIDTH  <- 12
FIG2B_HEIGHT <- 4.8

FIG2C_WIDTH  <- 16
FIG2C_HEIGHT <- 6.0

# DPI for PNG
PNG_DPI <- 600


############################
# 4. Colours
############################

# Geographic distribution ranges
range_colors <- c(
  "East"   = "#0072B2",
  "Centre" = "#009E73",
  "West"   = "#D55E00"
)


############################
# 5. Read metadata
############################

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

cat("\nMetadata loaded:\n")
cat("Rows:", nrow(metadata), "\n")
cat("Columns:", ncol(metadata), "\n")


############################
# 6. Prepare population-level PCA data
############################

# Climate PC values are population/site-level variables, and the metadata
# contains multiple RNA samples per population.
#
# Therefore, retain one row per population.
#
# We check that climate values are consistent within populations.

check_climate <- metadata %>%
  group_by(Population) %>%
  summarise(
    n_PC1 = n_distinct(Climate_PC1),
    n_PC2 = n_distinct(Climate_PC2),
    n_PC3 = n_distinct(Climate_PC3),
    n_range = n_distinct(Site_of_distribution_range),
    .groups = "drop"
  )

if (any(check_climate$n_PC1 > 1 |
        check_climate$n_PC2 > 1 |
        check_climate$n_PC3 > 1)) {

  warning(
    "At least one population has inconsistent climate PC values across samples."
  )
}

if (any(check_climate$n_range > 1)) {

  warning(
    "At least one population has more than one distribution-range category."
  )
}


# Keep one observation per population
population_data <- metadata %>%
  select(
    Population,
    Country,
    Country_code,
    Site_of_distribution_range,
    Locality,
    Climate_PC1,
    Climate_PC2,
    Climate_PC3,
    Latitude,
    Longitude,
    ECe_mean
  ) %>%
  distinct(Population, .keep_all = TRUE) %>%
  mutate(
    # Recode "Center" to "Centre" for display in the figure
    Site_of_distribution_range = recode(
      Site_of_distribution_range,
      "Center" = "Centre"
    ),

    Site_of_distribution_range =
      factor(
        Site_of_distribution_range,
        levels = c("East", "Centre", "West")
      )
  )


cat("\nPopulation-level dataset:\n")
cat("Number of populations:", nrow(population_data), "\n")

cat("\nPopulations by range:\n")
print(table(population_data$Site_of_distribution_range))


############################
# 7. Climate PCA summary
############################

# PCA summary from the climate PCA.
#
# proportion_variance is the fraction of total variance explained
# by each principal component.
#
# cumulative_variance is the cumulative proportion explained.

pca_summary <- tibble(

  PC = paste0(
    "PC",
    1:12
  ),

  eigenvalue = c(
    4.482520,
    3.234299,
    1.904368,
    1.224843,
    0.582664,
    0.278977,
    0.158380,
    0.096841,
    0.023843,
    0.008819,
    0.004101,
    0.000343
  ),

  proportion_variance = c(
    0.373543,
    0.269525,
    0.158697,
    0.102070,
    0.048555,
    0.023248,
    0.013198,
    0.008070,
    0.001987,
    0.000735,
    0.000342,
    0.0000286
  ),

  cumulative_variance = c(
    0.373543,
    0.643068,
    0.801766,
    0.903836,
    0.952391,
    0.975639,
    0.988838,
    0.996908,
    0.998895,
    0.999630,
    0.999971,
    1.000000
  )
)


############################
# Percentage variance explained for PC1-PC3
############################

explained_variance <- setNames(
  round(
    pca_summary$proportion_variance[
      pca_summary$PC %in% c("PC1", "PC2", "PC3")
    ] * 100,
    1
  ),
  pca_summary$PC[
    pca_summary$PC %in% c("PC1", "PC2", "PC3")
  ]
)

cat("\nExplained variance used in figures:\n")
print(explained_variance)


cat("\nPCA summary:\n")
print(pca_summary)


cat("\nVariance explained by PC1-PC3:\n")
print(
  pca_summary %>%
    filter(
      PC %in% c(
        "PC1",
        "PC2",
        "PC3"
      )
    )
)

############################
# 8. Function for PCA score panels
############################

make_pca_panel <- function(data, x_pc, y_pc, x_label, y_label) {

  ggplot(
    data,
    aes(
      x = .data[[x_pc]],
      y = .data[[y_pc]],
      colour = Site_of_distribution_range
    )
  ) +

    # Population points
    geom_point(
      size = 3.2,
      stroke = 0.4
    ) +

    # Population labels
    geom_text_repel(
      aes(label = Population),
      size = 3.0,
      show.legend = FALSE,

      # Give labels more space to move
      box.padding = 0.7,
      point.padding = 0.4,

      # Increase repulsion between labels
      force = 1.5,
      force_pull = 0.2,

      # Label every population
      max.overlaps = Inf,

      # Draw a segment between point and displaced label
      min.segment.length = 0,

      # Reproducible placement
      seed = 123
    ) +

    # Zero reference lines
    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "grey65"
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "grey65"
    ) +

    scale_colour_manual(
      values = range_colors,
      name = "Distribution range"
    ) +

    # Extra space for labels
    scale_x_continuous(
      expand = expansion(mult = 0.18)
    ) +

    scale_y_continuous(
      expand = expansion(mult = 0.18)
    ) +

    labs(
      x = x_label,
      y = y_label
    ) +

    theme_classic(base_size = 11) +

    theme(

      legend.position = "bottom",

      axis.title = element_text(
        size = 11,
        colour = "black"
      ),

      axis.text = element_text(
        size = 9.5,
        colour = "black"
      ),

      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),

      plot.margin = margin(
        t = 10,
        r = 18,
        b = 10,
        l = 10
      )
    )
}

############################
# 9. Axis labels
############################

PC1_label <- paste0(
  "Climate PC1 (37.4%)"
)

PC2_label <- paste0(
  "Climate PC2 (27.0%)"
)

PC3_label <- paste0(
  "Climate PC3 (15.9%)"
)


############################
# 10. Create Fig. 2B
############################

p_PC1_PC2 <- make_pca_panel(
  population_data,
  "Climate_PC1",
  "Climate_PC2",
  PC1_label,
  PC2_label
)

p_PC1_PC3 <- make_pca_panel(
  population_data,
  "Climate_PC1",
  "Climate_PC3",
  PC1_label,
  PC3_label
)

p_PC2_PC3 <- make_pca_panel(
  population_data,
  "Climate_PC2",
  "Climate_PC3",
  PC2_label,
  PC3_label
)


# Combine panels
fig2b <- (
  p_PC1_PC2 |
  p_PC1_PC3 |
  p_PC2_PC3
) +

  plot_layout(
    guides = "collect"
  ) &

  theme(
    legend.position = "bottom"
  )

############################
# 11. Save Fig. 2B
############################

fig2b_pdf <- file.path(
  output_dir,
  "Fig2B_climate_PCA_scores.pdf"
)

fig2b_png <- file.path(
  output_dir,
  "Fig2B_climate_PCA_scores.png"
)


ggsave(
  filename = fig2b_pdf,
  plot = fig2b,
  width = FIG2B_WIDTH,
  height = FIG2B_HEIGHT,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = fig2b_png,
  plot = fig2b,
  width = FIG2B_WIDTH,
  height = FIG2B_HEIGHT,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)


cat("\nFig. 2B saved:\n")
cat(fig2b_pdf, "\n")
cat(fig2b_png, "\n")


###############################################################################
# FIGURE 2C
# PCA LOADINGS
###############################################################################


############################
# 12. Read PCA loadings
############################

loadings <- read_csv(
  loadings_file,
  show_col_types = FALSE
)

cat("\nLoadings loaded:\n")
print(loadings)


############################
# 13. Bioclim variable names
############################

# Standard WorldClim / BIOCLIM variable descriptions.

bioclim_names <- c(

  bio1  = "Annual Mean Temperature",
  bio2  = "Mean Diurnal Range",
  bio3  = "Isothermality",
  bio4  = "Temperature Seasonality",
  bio5  = "Maximum Temperature of Warmest Month",
  bio6  = "Minimum Temperature of Coldest Month",
  bio7  = "Temperature Annual Range",
  bio8  = "Mean Temperature of Wettest Quarter",
  bio9  = "Mean Temperature of Driest Quarter",
  bio10 = "Mean Temperature of Warmest Quarter",
  bio11 = "Mean Temperature of Coldest Quarter",
  bio12 = "Annual Precipitation",
  bio13 = "Precipitation of Wettest Month",
  bio14 = "Precipitation of Driest Month",
  bio15 = "Precipitation Seasonality",
  bio16 = "Precipitation of Wettest Quarter",
  bio17 = "Precipitation of Driest Quarter",
  bio18 = "Precipitation of Warmest Quarter",
  bio19 = "Precipitation of Coldest Quarter"
)


############################
# 14. Convert loadings to long format
############################

loadings_long <- loadings %>%

  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "loading"
  ) %>%

  filter(
    PC %in% c("PC1", "PC2", "PC3")
  ) %>%

  mutate(
    abs_loading = abs(loading),

    variable_label = if_else(
      variable %in% names(bioclim_names),
      paste0(
        variable,
        " — ",
        bioclim_names[variable]
      ),
      variable
    ),

    direction = if_else(
      loading >= 0,
      "Positive",
      "Negative"
    )
  )


############################
# 15. Select top N variables
############################

top_loadings <- loadings_long %>%

  group_by(PC) %>%

  slice_max(
    order_by = abs_loading,
    n = TOP_N,
    with_ties = FALSE
  ) %>%

  ungroup()


############################
# 16. Make loading labels shorter
############################

# Concise labels for the publication figure.
# The BIO number is retained so the variables remain unambiguous.

short_bioclim_names <- c(
  bio1  = "Annual mean temperature",
  bio2  = "Mean diurnal range",
  bio3  = "Isothermality",
  bio4  = "Temperature seasonality",
  bio5  = "Maximum temperature, warmest month",
  bio6  = "Minimum temperature, coldest month",
  bio7  = "Temperature annual range",
  bio8  = "Wettest-quarter temperature",
  bio9  = "Driest-quarter temperature",
  bio10 = "Warmest-quarter temperature",
  bio11 = "Coldest-quarter temperature",
  bio12 = "Annual precipitation",
  bio13 = "Wettest-month precipitation",
  bio14 = "Driest-month precipitation",
  bio15 = "Precipitation seasonality",
  bio16 = "Wettest-quarter precipitation",
  bio17 = "Driest-quarter precipitation",
  bio18 = "Warmest-quarter precipitation",
  bio19 = "Coldest-quarter precipitation"
)

top_loadings <- top_loadings %>%

  mutate(
    variable_label = if_else(
      variable %in% names(short_bioclim_names),
      paste0(
        variable,
        " — ",
        short_bioclim_names[variable]
      ),
      variable
    )
  )


############################
# 17. Order variables within each PC
############################

top_loadings <- top_loadings %>%

  group_by(PC) %>%

  arrange(
    loading,
    .by_group = TRUE
  ) %>%

  mutate(
    variable_label = factor(
      variable_label,
      levels = variable_label
    )
  ) %>%

  ungroup()

############################
# 18. Function for loading plot
############################

make_loading_panel <- function(data, pc) {

  plot_data <- data %>%
    filter(PC == pc) %>%
    arrange(loading) %>%
    mutate(
      variable_label = factor(
        variable_label,
        levels = variable_label
      ),

      # Position of loading-value labels
      label_x = loading +
        if_else(
          loading >= 0,
          0.025,
          -0.025
        )
    )

  ggplot(
    plot_data,
    aes(
      x = loading,
      y = variable_label,
      fill = direction
    )
  ) +

    geom_col(
      width = 0.65,
      colour = "black",
      linewidth = 0.25
    ) +

    # Zero reference line
    geom_vline(
      xintercept = 0,
      linewidth = 0.5,
      colour = "black"
    ) +

    # Loading values
    geom_text(
      aes(
        x = label_x,
        label = sprintf("%.2f", loading),
        hjust = if_else(
          loading >= 0,
          0,
          1
        )
      ),
      size = 3.0,
      colour = "black"
    ) +

    scale_fill_manual(
      values = c(
        "Positive" = "#D95F02",
        "Negative" = "#1B9E77"
      )
    ) +

    # Common loading scale across all three PCs
    scale_x_continuous(
        breaks = seq(
          -0.4,
          0.4,
          by = 0.2
        ),
        labels = function(x) {
          sprintf("%.1f", x)
        },
        expand = expansion(mult = 0.08)
      ) +

      coord_cartesian(
        xlim = c(-0.60, 0.60),
        clip = "off"
      ) +
	  
    labs(
      x = "PCA loading",
      y = NULL,

      title = paste0(
        pc,
        " (",
        sprintf(
          "%.1f",
          explained_variance[pc]
        ),
        "% variance explained)"
      )
    ) +

    theme_classic(base_size = 11) +

    theme(

      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5,
        margin = margin(
          t = 0,
          r = 0,
          b = 10,
          l = 0
        )
      ),

      axis.text.y = element_text(
        size = 8.5,
        colour = "black"
      ),

      axis.text.x = element_text(
        size = 9,
        colour = "black"
      ),

      axis.title.x = element_text(
        size = 10,
        colour = "black"
      ),

      legend.position = "none",

      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),

      plot.margin = margin(
        t = 10,
        r = 12,
        b = 10,
        l = 10
      )
    )
}

############################
# 19. Create the three loading panels
############################

p_load_PC1 <- make_loading_panel(
  top_loadings,
  "PC1"
)

p_load_PC2 <- make_loading_panel(
  top_loadings,
  "PC2"
)

p_load_PC3 <- make_loading_panel(
  top_loadings,
  "PC3"
)


############################
# 20. Combine the three loading panels
############################

fig2c <- (
  p_load_PC1 |
  p_load_PC2 |
  p_load_PC3
) +
  plot_layout(
    widths = c(1, 1, 1)
  )


############################
# 21. Save Fig. 2C
############################

fig2c_pdf <- file.path(
  output_dir,
  "Fig2C_climate_PCA_loadings.pdf"
)

fig2c_png <- file.path(
  output_dir,
  "Fig2C_climate_PCA_loadings.png"
)


ggsave(
  filename = fig2c_pdf,
  plot = fig2c,
  width = FIG2C_WIDTH,
  height = FIG2C_HEIGHT,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = fig2c_png,
  plot = fig2c,
  width = FIG2C_WIDTH,
  height = FIG2C_HEIGHT,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)


cat("\nFig. 2C saved:\n")
cat(fig2c_pdf, "\n")
cat(fig2c_png, "\n")


###############################################################################
# 22. Export selected loadings as a CSV
###############################################################################

selected_loadings_file <- file.path(
  output_dir,
  "Fig2C_top_loadings_PC1_PC2_PC3.csv"
)

top_loadings_export <- top_loadings %>%
  select(
    PC,
    variable,
    variable_label,
    loading,
    abs_loading,
    direction
  ) %>%
  arrange(
    PC,
    desc(abs_loading)
  )

write_csv(
  top_loadings_export,
  selected_loadings_file
)

cat("\nSelected loadings saved:\n")
cat(selected_loadings_file, "\n")