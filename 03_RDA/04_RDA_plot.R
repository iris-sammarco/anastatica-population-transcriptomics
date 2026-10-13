#!/usr/bin/env Rscript

# ============================================================
#
# Author: Iris Sammarco
# Date: 06/10/2026
#
# ============================================================

############################################################
# Separate RDAs: Environment and Spatial Structure
#
# Full models:
#   Environment: Climate PC1 + PC2 + PC3 + Salinity
#   Space:       MEM2 + MEM3
#
# The RDA models include ALL predictors.
# Only predictors significant in marginal tests are displayed
# as arrows in the figures.
#
# Run with:
# qsub 04_RDA_plot_run.sh
#
############################################################


############################
# 1. Packages
############################

library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(patchwork)
library(scales)

############################
# 2. Paths
############################

analysis_dir <- "/Anastatica/IBD_IBE_transcriptome_analysis"

rds_file <- file.path(
  analysis_dir,
  "IBD_IBE_complete_analysis_objects.rds"
)

output_dir <- analysis_dir

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


############################
# 3. Load analysis objects
############################

cat("\nLoading analysis objects...\n")

obj <- readRDS(rds_file)

cat("Objects available in RDS:\n")
print(names(obj))


############################
# 4. Extract objects
############################

rda_response <- obj$rda_response

environment_data <- obj$environment_data

overall_selected_mem_matrix <- obj$overall_selected_mem_matrix


############################
# 5. Define predictors
############################

environment_predictors <- c(
  "Climate_PC1_z",
  "Climate_PC2_z",
  "Climate_PC3_z",
  "Salinity_z"
)

spatial_predictors <- c(
  "MEM2",
  "MEM3"
)


############################
# 6. Display labels for predictors
############################

# These names are used ONLY in the figures.
# The statistical models continue to use the original
# variable names above.

predictor_labels <- c(
  "Climate_PC1_z" = "Climate PC1",
  "Climate_PC2_z" = "Climate PC2",
  "Climate_PC3_z" = "Climate PC3",
  "Salinity_z"    = "Salinity",
  "MEM2"          = "MEM2",
  "MEM3"          = "MEM3"
)


############################
# 7. Check that predictors exist
############################

missing_environment <- setdiff(
  environment_predictors,
  colnames(environment_data)
)

if (length(missing_environment) > 0) {

  stop(
    "Missing environmental predictors: ",
    paste(missing_environment, collapse = ", ")
  )

}


missing_spatial <- setdiff(
  spatial_predictors,
  colnames(overall_selected_mem_matrix)
)

if (length(missing_spatial) > 0) {

  stop(
    "Missing spatial predictors: ",
    paste(missing_spatial, collapse = ", ")
  )

}


############################
# 8. Check row alignment
############################

if (
  !is.null(rownames(rda_response)) &&
  !is.null(rownames(environment_data))
) {

  if (
    !identical(
      rownames(rda_response),
      rownames(environment_data)
    )
  ) {

    stop(
      "Row names of rda_response and environment_data do not match."
    )

  }

}


if (
  !is.null(rownames(rda_response)) &&
  !is.null(rownames(overall_selected_mem_matrix))
) {

  if (
    !identical(
      rownames(rda_response),
      rownames(overall_selected_mem_matrix)
    )
  ) {

    stop(
      "Row names of rda_response and overall_selected_mem_matrix do not match."
    )

  }

}


############################
# 9. Reconstruct full RDAs
############################

cat("\nFitting environmental RDA...\n")

environment_df <- as.data.frame(
  environment_data[
    ,
    environment_predictors,
    drop = FALSE
  ]
)

environment_formula <- as.formula(
  paste(
    "rda_response ~",
    paste(
      environment_predictors,
      collapse = " + "
    )
  )
)

environment_rda <- vegan::rda(
  environment_formula,
  data = environment_df
)


cat("Fitting spatial RDA...\n")

spatial_df <- as.data.frame(
  overall_selected_mem_matrix[
    ,
    spatial_predictors,
    drop = FALSE
  ]
)

spatial_formula <- as.formula(
  paste(
    "rda_response ~",
    paste(
      spatial_predictors,
      collapse = " + "
    )
  )
)

spatial_rda <- vegan::rda(
  spatial_formula,
  data = spatial_df
)


############################
# 10. Global permutation tests
############################

set.seed(123)

cat("\n========================================\n")
cat("ENVIRONMENTAL RDA\n")
cat("========================================\n")

environment_global_test <- anova(
  environment_rda,
  permutations = 9999
)

print(environment_global_test)

cat("\nEnvironmental adjusted R2:\n")

print(
  RsquareAdj(environment_rda)$adj.r.squared
)


cat("\n========================================\n")
cat("SPATIAL RDA\n")
cat("========================================\n")

set.seed(123)

spatial_global_test <- anova(
  spatial_rda,
  permutations = 9999
)

print(spatial_global_test)

cat("\nSpatial adjusted R2:\n")

print(
  RsquareAdj(spatial_rda)$adj.r.squared
)


############################
# 11. Marginal tests
############################

cat("\n========================================\n")
cat("ENVIRONMENTAL MARGINAL TESTS\n")
cat("========================================\n")

set.seed(123)

environment_marginal <- anova(
  environment_rda,
  by = "margin",
  permutations = 9999
)

print(environment_marginal)


cat("\n========================================\n")
cat("SPATIAL MARGINAL TESTS\n")
cat("========================================\n")

set.seed(123)

spatial_marginal <- anova(
  spatial_rda,
  by = "margin",
  permutations = 9999
)

print(spatial_marginal)


############################
# 12. Extract marginal P-values
############################

environment_marginal_df <- as.data.frame(
  environment_marginal
) %>%
  rownames_to_column("Predictor")


spatial_marginal_df <- as.data.frame(
  spatial_marginal
) %>%
  rownames_to_column("Predictor")


p_column_environment <- grep(
  "^Pr\\(>F\\)$",
  colnames(environment_marginal_df),
  value = TRUE
)

p_column_spatial <- grep(
  "^Pr\\(>F\\)$",
  colnames(spatial_marginal_df),
  value = TRUE
)


if (length(p_column_environment) != 1) {

  stop(
    "Could not identify the P-value column for environmental marginal tests."
  )

}


if (length(p_column_spatial) != 1) {

  stop(
    "Could not identify the P-value column for spatial marginal tests."
  )

}


############################
# 13. Identify significant predictors
############################

environment_significant <- environment_marginal_df %>%
  filter(
    Predictor %in% environment_predictors,
    .data[[p_column_environment]] < 0.05
  ) %>%
  pull(Predictor)


spatial_significant <- spatial_marginal_df %>%
  filter(
    Predictor %in% spatial_predictors,
    .data[[p_column_spatial]] < 0.05
  ) %>%
  pull(Predictor)


cat("\n========================================\n")
cat("SIGNIFICANT ENVIRONMENTAL PREDICTORS\n")
cat("========================================\n")

print(environment_significant)


cat("\n========================================\n")
cat("SIGNIFICANT SPATIAL PREDICTORS\n")
cat("========================================\n")

print(spatial_significant)


############################
# 14. Save statistical results
############################

write.csv(
  environment_marginal_df,
  file.path(
    output_dir,
    "RDA_environment_marginal_tests.csv"
  ),
  row.names = FALSE
)


write.csv(
  spatial_marginal_df,
  file.path(
    output_dir,
    "RDA_spatial_marginal_tests.csv"
  ),
  row.names = FALSE
)


############################
# 15. Extract adjusted R2
############################

environment_adjR2 <-
  RsquareAdj(environment_rda)$adj.r.squared

spatial_adjR2 <-
  RsquareAdj(spatial_rda)$adj.r.squared


############################
# 16. Extract RDA axis percentages
############################

environment_eigenvalues <- eigenvals(
  environment_rda,
  model = "constrained"
)

spatial_eigenvalues <- eigenvals(
  spatial_rda,
  model = "constrained"
)


# Total inertia

environment_total_variance <-
  environment_rda$tot.chi

spatial_total_variance <-
  spatial_rda$tot.chi


# Percentage of total variance explained by each
# constrained RDA axis

environment_axis_percent <-
  100 *
  environment_eigenvalues /
  environment_total_variance


spatial_axis_percent <-
  100 *
  spatial_eigenvalues /
  spatial_total_variance


############################
# 17. Extract site scores
############################

environment_sites <- scores(
  environment_rda,
  display = "sites",
  scaling = 2,
  choices = 1:2
) %>%
  as.data.frame() %>%
  rownames_to_column("Sample")


spatial_sites <- scores(
  spatial_rda,
  display = "sites",
  scaling = 2,
  choices = 1:2
) %>%
  as.data.frame() %>%
  rownames_to_column("Sample")


############################
# 18. Population metadata
############################

site_codes <- c(
  "MO",
  "IS",
  "IR",
  "SP",
  "SA",
  "OM",
  "AE",
  "EG",
  "AL"
)


country_names <- c(
  "Morocco",
  "Israel",
  "Iran",
  "Spain (Canaries)",
  "Saudi Arabia",
  "Oman",
  "Arab Emirates",
  "Egypt",
  "Algeria"
)


country_lookup <- setNames(
  country_names,
  site_codes
)


############################
# 19. Assign country from sample names
############################

assign_country <- function(sample_names) {

  result <- rep(
    NA_character_,
    length(sample_names)
  )

  for (i in seq_along(site_codes)) {

    result[
      grepl(
        paste0(
          "^",
          site_codes[i]
        ),
        sample_names
      )
    ] <- country_lookup[
      site_codes[i]
    ]

  }

  result
}


environment_sites$Country <-
  assign_country(
    environment_sites$Sample
  )


spatial_sites$Country <-
  assign_country(
    spatial_sites$Sample
  )

############################
# Population labels
############################

# Use the original population/sample IDs as the
# labels displayed next to the RDA points.

environment_sites$Population <-
  environment_sites$Sample


spatial_sites$Population <-
  spatial_sites$Sample

############################
# 20. Check unmatched samples
############################

if (
  any(
    is.na(
      environment_sites$Country
    )
  )
) {

  cat(
    "\nWARNING: Environmental RDA samples without country:\n"
  )

  print(
    environment_sites$Sample[
      is.na(
        environment_sites$Country
      )
    ]
  )

}


if (
  any(
    is.na(
      spatial_sites$Country
    )
  )
) {

  cat(
    "\nWARNING: Spatial RDA samples without country:\n"
  )

  print(
    spatial_sites$Sample[
      is.na(
        spatial_sites$Country
      )
    ]
  )

}


############################
# 21. Extract full-model biplot scores
############################

# IMPORTANT:
#
# These scores come from the FULL RDA models.
#
# We only filter which arrows are DISPLAYED.
#
# The RDA models themselves still contain all predictors.


environment_arrows <- scores(
  environment_rda,
  display = "bp",
  scaling = 2,
  choices = 1:2
) %>%
  as.data.frame() %>%
  rownames_to_column("Predictor")


spatial_arrows <- scores(
  spatial_rda,
  display = "bp",
  scaling = 2,
  choices = 1:2
) %>%
  as.data.frame() %>%
  rownames_to_column("Predictor")


############################
# 22. Keep ONLY significant arrows
############################

environment_arrows_plot <- environment_arrows %>%
  filter(
    Predictor %in% environment_significant
  ) %>%
  mutate(
    Label = predictor_labels[Predictor]
  )


spatial_arrows_plot <- spatial_arrows %>%
  filter(
    Predictor %in% spatial_significant
  ) %>%
  mutate(
    Label = predictor_labels[Predictor]
  )

############################
# 23. Arrow scaling function
############################

# This scales arrows for visual readability.
# It does NOT change the underlying RDA scores.

scale_arrows <- function(
  arrows,
  sites,
  multiplier = 0.75
) {

  if (
    nrow(arrows) == 0
  ) {

    return(arrows)

  }


  max_site <- max(
    abs(
      c(
        sites$RDA1,
        sites$RDA2
      )
    ),
    na.rm = TRUE
  )


  max_arrow <- max(
    abs(
      c(
        arrows$RDA1,
        arrows$RDA2
      )
    ),
    na.rm = TRUE
  )


  if (
    max_arrow == 0
  ) {

    return(arrows)

  }


  scale_factor <- (
    max_site * multiplier
  ) / max_arrow


  arrows$RDA1 <-
    arrows$RDA1 * scale_factor

  arrows$RDA2 <-
    arrows$RDA2 * scale_factor


  arrows

}


environment_arrows_plot <-
  scale_arrows(
    environment_arrows_plot,
    environment_sites
  )


spatial_arrows_plot <-
  scale_arrows(
    spatial_arrows_plot,
    spatial_sites
  )


############################
# 23b. Predictor label positions
############################

environment_arrows_plot <-
  environment_arrows_plot %>%
  mutate(
    Label_RDA1 = RDA1,
    Label_RDA2 = RDA2
  )


spatial_arrows_plot <-
  spatial_arrows_plot %>%
  mutate(
    Label_RDA1 = RDA1,
    Label_RDA2 = RDA2
  )
  
############################
# 24. Country factors
############################

country_levels <- country_names


environment_sites$Country <-
  factor(
    environment_sites$Country,
    levels = country_levels
  )


spatial_sites$Country <-
  factor(
    spatial_sites$Country,
    levels = country_levels
  )

############################
# 25. Plotting function
############################

make_rda_plot <- function(
  sites,
  arrows,
  x_percent,
  y_percent,
  title,
  panel_letter,
  arrow_colour
) {

  p <- ggplot(
    sites,
    aes(
      x = RDA1,
      y = RDA2
    )
  ) +

    ##########################################################
    # Sample points
    ##########################################################

    geom_point(
      aes(
        fill = Country
      ),
      shape = 21,
      size = 3.0,
      alpha = 0.85,
      colour = "black",
      stroke = 0.35
    ) +


    ##########################################################
    # Zero horizontal axis
    ##########################################################

    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      colour = "grey75"
    ) +


    ##########################################################
    # Zero vertical axis
    ##########################################################

    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = "grey75"
    )


  ############################################################
  # Add arrows only if significant predictors exist
  ############################################################

  if (
    nrow(arrows) > 0
  ) {

    p <- p +

      ########################################################
      # Predictor arrows
      ########################################################

      geom_segment(
        data = arrows,

        aes(
          x = 0,
          y = 0,
          xend = RDA1,
          yend = RDA2
        ),

        inherit.aes = FALSE,

        arrow = arrow(
          length = unit(
            0.18,
            "cm"
          ),
          type = "closed"
        ),

        linewidth = 0.65,

        colour = arrow_colour
      ) +


      ########################################################
      # Predictor labels
      ########################################################

      geom_text_repel(

        data = arrows,

        aes(
          x = Label_RDA1,
          y = Label_RDA2,
          label = Label
        ),

        inherit.aes = FALSE,

        size = 3.5,

        fontface = "bold",

        colour = arrow_colour,

        seed = 123,

        force = 1.2,

        max.overlaps = Inf,

        box.padding = 0.4,

        point.padding = 0.2,

        min.segment.length = 0,

        segment.size = 0.25,

        segment.alpha = 0.4

      )

  }


  ############################################################
  # Population labels
  ############################################################

  p <- p +

    geom_text_repel(

      aes(
        x = RDA1,
        y = RDA2,
        label = Population
      ),

      size = 3.1,

      colour = "grey15",

      fontface = "plain",

      seed = 123,

      force = 1.5,

      max.overlaps = Inf,

      box.padding = 0.25,

      point.padding = 0.15,

      min.segment.length = 0,

      segment.size = 0.25,

      segment.alpha = 0.45,

      direction = "both"

    )


  ############################################################
  # Final plot formatting
  ############################################################

  p +

    labs(

      title = title,

      x = paste0(
        "RDA1 (",
        round(
          x_percent,
          1
        ),
        "%)"
      ),

      y = paste0(
        "RDA2 (",
        round(
          y_percent,
          1
        ),
        "%)"
      ),

      fill = "Country"

    ) +


    ##########################################################
    # Panel letter
    ##########################################################

    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = panel_letter,
      hjust = 1.8,
      vjust = 1.2,
      size = 5,
      fontface = "bold"
    ) +


    ##########################################################
    # Theme
    ##########################################################

    theme_classic(
      base_size = 12
    ) +

    theme(

      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5,
        margin = margin(
          b = 10
        )
      ),

      axis.title = element_text(
        size = 12
      ),

      axis.text = element_text(
        size = 10,
        colour = "grey20"
      ),

      legend.position = "right",

      legend.title = element_text(
        face = "bold",
        size = 11
      ),

      legend.text = element_text(
        size = 10
      ),

      panel.border = element_rect(
        fill = NA,
        colour = "grey30",
        linewidth = 0.5
      ),

      plot.margin = margin(
        15,
        15,
        10,
        35
      )

    )

}

############################
# 26. Create environmental plot
############################

environment_plot <- make_rda_plot(

  sites =
    environment_sites,

  arrows =
    environment_arrows_plot,

  x_percent =
    environment_axis_percent[1],

  y_percent =
    environment_axis_percent[2],

  title =
    "Environmental RDA",

  panel_letter =
    "(b)",

  arrow_colour =
    "#41ab5d"

)


############################
# 27. Create spatial plot
############################

spatial_plot <- make_rda_plot(

  sites =
    spatial_sites,

  arrows =
    spatial_arrows_plot,

  x_percent =
    spatial_axis_percent[1],

  y_percent =
    spatial_axis_percent[2],

  title =
    "Spatial RDA",

  panel_letter =
    "(a)",

  arrow_colour =
    "#E67E22"

)

############################
# 28. Display combined figure
############################

combined_figure <- (

  spatial_plot +
  environment_plot

) +

  plot_layout(
    guides = "collect"
  ) &

  theme(
    legend.position = "right"
  )


print(
  combined_figure
)


############################
# 29. Save combined figure
############################

pdf_file <- file.path(
  output_dir,
  "separate_RDA_environment_spatial_significant_marginal.pdf"
)


png_file <- file.path(
  output_dir,
  "separate_RDA_environment_spatial_significant_marginal.png"
)


ggsave(
  pdf_file,
  combined_figure,
  width = 15,
  height = 7.5,
  units = "in",
  device = "pdf"
)


ggsave(
  png_file,
  combined_figure,
  width = 15,
  height = 7.5,
  units = "in",
  dpi = 600
)


############################
# 30. Save individual figures
############################

ggsave(
  file.path(
    output_dir,
    "RDA_environment_significant_marginal.pdf"
  ),
  environment_plot,
  width = 8,
  height = 7,
  units = "in",
  device = "pdf"
)


ggsave(
  file.path(
    output_dir,
    "RDA_environment_significant_marginal.png"
  ),
  environment_plot,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)


ggsave(
  file.path(
    output_dir,
    "RDA_spatial_significant_marginal.pdf"
  ),
  spatial_plot,
  width = 8,
  height = 7,
  units = "in",
  device = "pdf"
)


ggsave(
  file.path(
    output_dir,
    "RDA_spatial_significant_marginal.png"
  ),
  spatial_plot,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)


############################
# 31. Save RDA objects
############################

saveRDS(
  environment_rda,
  file.path(
    output_dir,
    "environment_RDA_full_model.rds"
  )
)


saveRDS(
  spatial_rda,
  file.path(
    output_dir,
    "spatial_RDA_full_model.rds"
  )
)


############################
# 32. Save site scores
############################

write.csv(
  environment_sites,
  file.path(
    output_dir,
    "environment_RDA_site_scores.csv"
  ),
  row.names = FALSE
)


write.csv(
  spatial_sites,
  file.path(
    output_dir,
    "spatial_RDA_site_scores.csv"
  ),
  row.names = FALSE
)


############################
# 33. Save plotted arrows
############################

write.csv(
  environment_arrows_plot,
  file.path(
    output_dir,
    "environment_RDA_significant_arrows.csv"
  ),
  row.names = FALSE
)


write.csv(
  spatial_arrows_plot,
  file.path(
    output_dir,
    "spatial_RDA_significant_arrows.csv"
  ),
  row.names = FALSE
)


############################
# 34. Print final summary
############################

cat("\n\n")

cat(
  "====================================================\n"
)

cat(
  "RDA ANALYSIS COMPLETE\n"
)

cat(
  "====================================================\n\n"
)


cat(
  "Environmental RDA\n"
)

cat(
  "-----------------\n"
)

cat(
  "Adjusted R2: ",
  round(
    environment_adjR2,
    4
  ),
  "\n",
  sep = ""
)


cat(
  "Significant marginal predictors: ",
  ifelse(
    length(environment_significant) > 0,

    paste(
      environment_significant,
      collapse = ", "
    ),

    "NONE"
  ),

  "\n\n",
  sep = ""
)


cat(
  "Spatial RDA\n"
)

cat(
  "-----------\n"
)

cat(
  "Adjusted R2: ",
  round(
    spatial_adjR2,
    4
  ),
  "\n",
  sep = ""
)


cat(
  "Significant marginal predictors: ",
  ifelse(
    length(spatial_significant) > 0,

    paste(
      spatial_significant,
      collapse = ", "
    ),

    "NONE"
  ),

  "\n\n",
  sep = ""
)


cat(
  "Output directory:\n"
)

cat(
  output_dir,
  "\n\n"
)


cat(
  "Combined figure:\n"
)

cat(
  pdf_file,
  "\n"
)

cat(
  png_file,
  "\n\n"
)


cat(
  "====================================================\n"
)

cat(
  "Finished.\n"
)

cat(
  "====================================================\n"
)