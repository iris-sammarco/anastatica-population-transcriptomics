#!/usr/bin/env Rscript

# ============================================================
# Author: Iris Sammarco
# Date: 29/10/2026
# ============================================================

############################################################
# Separate RDAs: Climate and Spatial Structure
#
# Input: RDA_complete_analysis_objects.rds, written by the
#        main RDA transcriptome analysis script.
#
# Full models:
#   Climate: Climate PC1 + PC2 + PC3
#   Space:   MEM2 + MEM3
#
# The RDA models include ALL predictors above.
# Figure arrows shown:
#   Climate RDA: Climate PC1 and Climate PC2 only
#   Spatial RDA: MEM2 and MEM3
#
# Figure features:
#   - Predictors displayed as "Climate PC1", "Climate PC2"
#   - No subtitles, no population centroid crosses
#   - Population labels with ggrepel
#   - Grey horizontal and vertical zero-reference lines
#   - Larger text for titles, axes and legend
#
# Run with:
# qsub 04_RDA_plot_run.sh
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


############################
# 2. Paths
############################

analysis_dir <- "/Anastatica/RDA_transcriptome_analysis"

rds_file <- file.path(analysis_dir, "RDA_complete_analysis_objects.rds")

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

rda_response                <- obj$rda_response
environment_data            <- obj$environment_data
overall_selected_mem_matrix <- obj$overall_selected_mem_matrix

if (ncol(overall_selected_mem_matrix) == 0) {
  stop("No overall spatial MEMs were selected in the main analysis.")
}


############################
# 5. Define predictors
############################

# Climate model
climate_predictors <- c(
  "Climate_PC1_z",
  "Climate_PC2_z",
  "Climate_PC3_z"
)

# Climate arrows displayed in the figure
climate_arrows_displayed <- c(
  "Climate_PC1_z",
  "Climate_PC2_z"
)

# Spatial model
spatial_predictors <- c(
  "MEM2",
  "MEM3"
)

# Spatial arrows displayed in the figure
spatial_arrows_displayed <- spatial_predictors


############################
# 6. Display labels for predictors
############################

# Used ONLY in the figures. The models use the original names.

predictor_labels <- c(
  "Climate_PC1_z" = "Climate PC1",
  "Climate_PC2_z" = "Climate PC2",
  "Climate_PC3_z" = "Climate PC3",
  "MEM2"          = "MEM2",
  "MEM3"          = "MEM3"
)


############################
# 7. Check that predictors exist
############################

missing_climate <- setdiff(climate_predictors, colnames(environment_data))

if (length(missing_climate) > 0) {
  stop(
    "Missing climate predictors: ",
    paste(missing_climate, collapse = ", ")
  )
}

missing_spatial <- setdiff(spatial_predictors, colnames(overall_selected_mem_matrix))

if (length(missing_spatial) > 0) {
  stop(
    "Missing spatial predictors: ",
    paste(missing_spatial, collapse = ", "),
    "\nAvailable MEMs: ",
    paste(colnames(overall_selected_mem_matrix), collapse = ", ")
  )
}


############################
# 8. Check row alignment
############################

if (!identical(rownames(rda_response), rownames(environment_data))) {
  stop("Row names of rda_response and environment_data do not match.")
}

if (!identical(rownames(rda_response), rownames(overall_selected_mem_matrix))) {
  stop("Row names of rda_response and overall_selected_mem_matrix do not match.")
}


############################
# 9. Fit full RDAs
############################

cat("\nFitting climate RDA...\n")

climate_df <- as.data.frame(
  environment_data[, climate_predictors, drop = FALSE]
)

climate_formula <- as.formula(
  paste("rda_response ~", paste(climate_predictors, collapse = " + "))
)

climate_rda <- vegan::rda(climate_formula, data = climate_df)


cat("Fitting spatial RDA...\n")

spatial_df <- as.data.frame(
  overall_selected_mem_matrix[, spatial_predictors, drop = FALSE]
)

spatial_formula <- as.formula(
  paste("rda_response ~", paste(spatial_predictors, collapse = " + "))
)

spatial_rda <- vegan::rda(spatial_formula, data = spatial_df)


############################
# 10. Global permutation tests
############################

cat("\n========================================\n")
cat("CLIMATE RDA\n")
cat("========================================\n")

set.seed(123)
climate_global_test <- anova(climate_rda, permutations = 9999)
print(climate_global_test)

cat("\nClimate adjusted R2:\n")
print(RsquareAdj(climate_rda)$adj.r.squared)


cat("\n========================================\n")
cat("SPATIAL RDA\n")
cat("========================================\n")

set.seed(123)
spatial_global_test <- anova(spatial_rda, permutations = 9999)
print(spatial_global_test)

cat("\nSpatial adjusted R2:\n")
print(RsquareAdj(spatial_rda)$adj.r.squared)


############################
# 11. Marginal tests
############################

cat("\n========================================\n")
cat("CLIMATE MARGINAL TESTS\n")
cat("========================================\n")

set.seed(123)
climate_marginal <- anova(climate_rda, by = "margin", permutations = 9999)
print(climate_marginal)


cat("\n========================================\n")
cat("SPATIAL MARGINAL TESTS\n")
cat("========================================\n")

set.seed(123)
spatial_marginal <- anova(spatial_rda, by = "margin", permutations = 9999)
print(spatial_marginal)


############################
# 12. Save statistical results
############################

climate_marginal_df <- as.data.frame(climate_marginal) %>%
  rownames_to_column("Predictor")

spatial_marginal_df <- as.data.frame(spatial_marginal) %>%
  rownames_to_column("Predictor")

write.csv(
  climate_marginal_df,
  file.path(output_dir, "RDA_climate_marginal_tests.csv"),
  row.names = FALSE
)

write.csv(
  spatial_marginal_df,
  file.path(output_dir, "RDA_spatial_marginal_tests.csv"),
  row.names = FALSE
)


############################
# 13. Adjusted R2
############################

climate_adjR2 <- RsquareAdj(climate_rda)$adj.r.squared
spatial_adjR2 <- RsquareAdj(spatial_rda)$adj.r.squared


############################
# 14. RDA axis percentages
############################

climate_eigenvalues <- eigenvals(climate_rda, model = "constrained")
spatial_eigenvalues <- eigenvals(spatial_rda, model = "constrained")

# Percentage of TOTAL variance explained by each constrained axis
climate_axis_percent <- 100 * climate_eigenvalues / climate_rda$tot.chi
spatial_axis_percent <- 100 * spatial_eigenvalues / spatial_rda$tot.chi


############################
# 15. Site scores
############################

climate_sites <- scores(
  climate_rda,
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
# 16. Population metadata
############################

site_codes <- c("MO", "IS", "IR", "SP", "SA", "OM", "AE", "EG", "AL")

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

country_lookup <- setNames(country_names, site_codes)


############################
# 17. Assign country from sample names
############################

assign_country <- function(sample_names) {

  result <- rep(NA_character_, length(sample_names))

  for (i in seq_along(site_codes)) {
    result[grepl(paste0("^", site_codes[i]), sample_names)] <-
      country_lookup[site_codes[i]]
  }

  result
}

climate_sites$Country <- assign_country(climate_sites$Sample)
spatial_sites$Country <- assign_country(spatial_sites$Sample)

# Population labels = original population IDs
climate_sites$Population <- climate_sites$Sample
spatial_sites$Population <- spatial_sites$Sample


############################
# 18. Check unmatched samples
############################

if (any(is.na(climate_sites$Country))) {
  cat("\nWARNING: Climate RDA samples without country:\n")
  print(climate_sites$Sample[is.na(climate_sites$Country)])
}

if (any(is.na(spatial_sites$Country))) {
  cat("\nWARNING: Spatial RDA samples without country:\n")
  print(spatial_sites$Sample[is.na(spatial_sites$Country)])
}


############################
# 19. Full-model biplot scores
############################

# Scores come from the FULL models; we only choose which
# arrows are DISPLAYED.

climate_arrows <- scores(
  climate_rda,
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
# 20. Select arrows to display
############################

climate_arrows_plot <- climate_arrows %>%
  filter(Predictor %in% climate_arrows_displayed) %>%
  mutate(Label = predictor_labels[Predictor])

spatial_arrows_plot <- spatial_arrows %>%
  filter(Predictor %in% spatial_arrows_displayed) %>%
  mutate(Label = predictor_labels[Predictor])


############################
# 21. Arrow scaling (visual only)
############################

scale_arrows <- function(arrows, sites, multiplier = 0.75) {

  if (nrow(arrows) == 0) return(arrows)

  max_site  <- max(abs(c(sites$RDA1, sites$RDA2)), na.rm = TRUE)
  max_arrow <- max(abs(c(arrows$RDA1, arrows$RDA2)), na.rm = TRUE)

  if (max_arrow == 0) return(arrows)

  scale_factor <- (max_site * multiplier) / max_arrow

  arrows$RDA1 <- arrows$RDA1 * scale_factor
  arrows$RDA2 <- arrows$RDA2 * scale_factor

  arrows
}

climate_arrows_plot <- scale_arrows(climate_arrows_plot, climate_sites)
spatial_arrows_plot <- scale_arrows(spatial_arrows_plot, spatial_sites)


############################
# 22. Predictor label positions
############################

climate_arrows_plot <- climate_arrows_plot %>%
  mutate(Label_RDA1 = RDA1, Label_RDA2 = RDA2)

spatial_arrows_plot <- spatial_arrows_plot %>%
  mutate(Label_RDA1 = RDA1, Label_RDA2 = RDA2)


############################
# 23. Country factors
############################

climate_sites$Country <- factor(climate_sites$Country, levels = country_names)
spatial_sites$Country <- factor(spatial_sites$Country, levels = country_names)


############################
# 24. Plotting function
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

  p <- ggplot(sites, aes(x = RDA1, y = RDA2)) +

    # Sample points
    geom_point(
      aes(fill = Country),
      shape  = 21,
      size   = 3.4,
      alpha  = 0.85,
      colour = "black",
      stroke = 0.35
    ) +

    # Zero reference lines
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey75")


  # Predictor arrows and labels
  if (nrow(arrows) > 0) {

    p <- p +

      geom_segment(
        data = arrows,
        aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
        inherit.aes = FALSE,
        arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
        linewidth = 0.7,
        colour = arrow_colour
      ) +

      geom_text_repel(
        data = arrows,
        aes(x = Label_RDA1, y = Label_RDA2, label = Label),
        inherit.aes = FALSE,
        size = 4.3,
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


  # Population labels
  p <- p +

    geom_text_repel(
      aes(x = RDA1, y = RDA2, label = Population),
      size = 3.4,
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


  # Final formatting
  p +

	labs(
      title = title,
      tag   = panel_letter,
      x = paste0("RDA1 (", round(x_percent, 1), "%)"),
      y = paste0("RDA2 (", round(y_percent, 1), "%)"),
      fill = "Country"
    ) +

    theme_classic(base_size = 14) +
	
	theme(
      plot.title = element_text(
        face = "bold",
        size = 19,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.tag = element_text(face = "bold", size = 22),
      plot.tag.position = c(0.02, 0.98),
      axis.title  = element_text(size = 16),
      axis.text   = element_text(size = 13, colour = "grey20"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 15),
      legend.text  = element_text(size = 13),
      panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.5),
      plot.margin  = margin(15, 15, 10, 15)
    )
}


############################
# 25. Climate plot
############################

climate_plot <- make_rda_plot(
  sites        = climate_sites,
  arrows       = climate_arrows_plot,
  x_percent    = climate_axis_percent[1],
  y_percent    = climate_axis_percent[2],
  title        = "Climate RDA",
  panel_letter = "(b)",
  arrow_colour = "#41ab5d"
)


############################
# 26. Spatial plot
############################

spatial_plot <- make_rda_plot(
  sites        = spatial_sites,
  arrows       = spatial_arrows_plot,
  x_percent    = spatial_axis_percent[1],
  y_percent    = spatial_axis_percent[2],
  title        = "Spatial RDA",
  panel_letter = "(a)",
  arrow_colour = "#E67E22"
)


############################
# 27. Combined figure
############################

combined_figure <- (spatial_plot + climate_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_figure)


############################
# 28. Save combined figure
############################

pdf_file <- file.path(output_dir, "separate_RDA_climate_spatial.pdf")
png_file <- file.path(output_dir, "separate_RDA_climate_spatial.png")

ggsave(pdf_file, combined_figure,
       width = 16, height = 8, units = "in", device = "pdf")

ggsave(png_file, combined_figure,
       width = 16, height = 8, units = "in", dpi = 600)


############################
# 29. Save individual figures
############################

ggsave(file.path(output_dir, "RDA_climate.pdf"), climate_plot,
       width = 8.5, height = 7.5, units = "in", device = "pdf")

ggsave(file.path(output_dir, "RDA_climate.png"), climate_plot,
       width = 8.5, height = 7.5, units = "in", dpi = 600)

ggsave(file.path(output_dir, "RDA_spatial.pdf"), spatial_plot,
       width = 8.5, height = 7.5, units = "in", device = "pdf")

ggsave(file.path(output_dir, "RDA_spatial.png"), spatial_plot,
       width = 8.5, height = 7.5, units = "in", dpi = 600)


############################
# 30. Save RDA objects
############################

saveRDS(climate_rda, file.path(output_dir, "climate_RDA_full_model.rds"))
saveRDS(spatial_rda, file.path(output_dir, "spatial_RDA_full_model.rds"))


############################
# 31. Save site scores
############################

write.csv(
  climate_sites,
  file.path(output_dir, "climate_RDA_site_scores.csv"),
  row.names = FALSE
)

write.csv(
  spatial_sites,
  file.path(output_dir, "spatial_RDA_site_scores.csv"),
  row.names = FALSE
)


############################
# 32. Save plotted arrows
############################

write.csv(
  climate_arrows_plot,
  file.path(output_dir, "climate_RDA_plotted_arrows.csv"),
  row.names = FALSE
)

write.csv(
  spatial_arrows_plot,
  file.path(output_dir, "spatial_RDA_plotted_arrows.csv"),
  row.names = FALSE
)


############################
# 33. Final summary
############################

cat("\n\n====================================================\n")
cat("RDA ANALYSIS COMPLETE\n")
cat("====================================================\n\n")

cat("Climate RDA (PC1 + PC2 + PC3)\n")
cat("-----------------------------\n")
cat("Adjusted R2: ", round(climate_adjR2, 4), "\n", sep = "")
cat("Arrows displayed: ",
    paste(predictor_labels[climate_arrows_displayed], collapse = ", "),
    "\n\n", sep = "")

cat("Spatial RDA (MEM2 + MEM3)\n")
cat("-------------------------\n")
cat("Adjusted R2: ", round(spatial_adjR2, 4), "\n", sep = "")
cat("Arrows displayed: ",
    paste(predictor_labels[spatial_arrows_displayed], collapse = ", "),
    "\n\n", sep = "")

cat("Output directory:\n", output_dir, "\n\n", sep = "")
cat("Combined figure:\n", pdf_file, "\n", png_file, "\n\n", sep = "")

cat("====================================================\n")
cat("Finished.\n")
cat("====================================================\n")
