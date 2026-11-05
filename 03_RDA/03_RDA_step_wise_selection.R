#!/usr/bin/env Rscript

# ============================================================
# AUTHOR: Iris Sammarco
# DATE: 03/11/2026
#
# FORWARD SELECTION OF CLIMATE PCs IN THE PARTIAL RDA
# (full dataset only)
#
# Model (full dataset):
#   response ~ Climate_PC1_z + Climate_PC2_z + Climate_PC3_z
#              + Condition(Salinity_z + selected overall dbMEMs)
#
# Steps:
#   A. Manual forward selection (permutation-based, p-value criterion)
#      - at every step, all candidate p-values / F / adj. R2 are saved
#   B. vegan::ordiR2step() (adjusted-R2 criterion + p-value),
#      run as an independent cross-check
#   C. Final model of the selected PCs (overall + marginal tests)
#
# Reads the objects saved by 01_IBD_IBE_transcriptome_RDA.R
# Writes everything to: RDA_transcriptome_analysis/step_selection/
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(readr)
  library(dplyr)
  library(tibble)
})

# ============================================================
# 1. USER SETTINGS
# ============================================================

main_output_dir <-
  "/Anastatica/RDA_transcriptome_analysis"

step_output_dir <- file.path(main_output_dir, "step_selection")
dir.create(step_output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(12345)

n_permutations <- 9999     # permutations for every test
alpha_in       <- 0.05     # entry threshold for forward selection

climate_candidates <- c("Climate_PC1_z", "Climate_PC2_z", "Climate_PC3_z")
salinity_predictor <- "Salinity_z"

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

run_forward_selection <- function(response, env, candidates, cond_names,
                                  label, outdir, nperm, alpha) {

  cat("\n============================================================\n")
  cat("FORWARD SELECTION:", label, "\n")
  cat("Populations:", nrow(response), "| Genes:", ncol(response), "\n")
  cat("Conditioning on:",
      if (length(cond_names)) paste(cond_names, collapse = " + ") else "(nothing)",
      "\n")
  cat("============================================================\n")

  cond_term <- if (length(cond_names) > 0) {
    paste0("Condition(", paste(cond_names, collapse = " + "), ")")
  } else NULL

  build_formula <- function(terms) {
    rhs <- c(terms, cond_term)
    if (length(rhs) == 0) rhs <- "1"
    as.formula(paste("response ~", paste(rhs, collapse = " + ")))
  }

  # ---------- Reference: null (conditions only) and full models ----------
  null_model <- vegan::rda(build_formula(character(0)), data = env)
  full_model <- vegan::rda(build_formula(candidates),   data = env)

  full_overall  <- anova(full_model, permutations = nperm)
  full_marginal <- anova(full_model, by = "margin", permutations = nperm)
  full_r2       <- vegan::RsquareAdj(full_model)

  write.csv(as.data.frame(full_overall),
            file.path(outdir, paste0(label, "_reference_full_model_overall_test.csv")))
  write.csv(as.data.frame(full_marginal),
            file.path(outdir, paste0(label, "_reference_full_model_marginal_tests.csv")))

  # ---------- A. Manual forward selection ----------
  selected  <- character(0)
  remaining <- candidates
  candidate_log <- list()
  path_log <- list()
  step <- 0

  while (length(remaining) > 0) {
    step <- step + 1

    cand_tbl <- bind_rows(lapply(remaining, function(v) {
      terms <- c(selected, v)
      m <- vegan::rda(build_formula(terms), data = env)
      a <- as.data.frame(anova(m, by = "margin", permutations = nperm))
      tibble(
        step = step,
        already_selected = ifelse(length(selected) > 0,
                                  paste(selected, collapse = " + "), "none"),
        candidate = v,
        Df = a[v, "Df"],
        Variance = a[v, "Variance"],
        F = a[v, "F"],
        P = a[v, "Pr(>F)"],
        adjR2_if_added = vegan::RsquareAdj(m)$adj.r.squared
      )
    }))

    cand_tbl <- cand_tbl |>
      arrange(P, desc(F))

    candidate_log[[step]] <- cand_tbl
    best <- cand_tbl[1, ]

    cat("\nStep", step, "- candidates (already in model:",
        best$already_selected, ")\n")
    print(as.data.frame(cand_tbl[, c("candidate", "F", "P",
                                     "adjR2_if_added")]),
          row.names = FALSE)

    if (is.finite(best$P) && best$P < alpha) {
      selected  <- c(selected, best$candidate)
      remaining <- setdiff(remaining, best$candidate)
      path_log[[step]] <- tibble(step = step, added = best$candidate,
                                 F = best$F, P = best$P,
                                 adjR2_after = best$adjR2_if_added,
                                 action = "added")
      cat("-> Added:", best$candidate, "(P =", format.pval(best$P, digits = 3), ")\n")
    } else {
      path_log[[step]] <- tibble(step = step, added = NA_character_,
                                 F = best$F, P = best$P,
                                 adjR2_after = NA_real_,
                                 action = paste0("stopped: best candidate (",
                                                 best$candidate, ") P >= ", alpha))
      cat("-> Stopped: best candidate", best$candidate,
          "has P =", format.pval(best$P, digits = 3), ">=", alpha, "\n")
      break
    }
  }

  candidate_tbl_all <- bind_rows(candidate_log)
  path_tbl <- bind_rows(path_log)

  write_csv(candidate_tbl_all,
            file.path(outdir, paste0(label, "_forward_all_candidate_tests.csv")))
  write_csv(path_tbl,
            file.path(outdir, paste0(label, "_forward_selection_path.csv")))

  # ---------- C. Final model ----------
  final_p <- NA_real_
  final_r2 <- tibble(R2 = NA_real_, adjusted_R2 = NA_real_)

  if (length(selected) > 0) {
    final_model <- vegan::rda(build_formula(selected), data = env)
    final_overall <- anova(final_model, permutations = nperm)
    final_marginal <- anova(final_model, by = "margin", permutations = nperm)
    rr <- vegan::RsquareAdj(final_model)
    final_r2 <- tibble(R2 = rr$r.squared, adjusted_R2 = rr$adj.r.squared)
    final_p <- as.numeric(as.data.frame(final_overall)[["Pr(>F)"]][1])

    write.csv(as.data.frame(final_overall),
              file.path(outdir, paste0(label, "_final_model_overall_test.csv")))
    write.csv(as.data.frame(final_marginal),
              file.path(outdir, paste0(label, "_final_model_marginal_tests.csv")))
    capture.output(summary(final_model),
                   file = file.path(outdir, paste0(label, "_final_model_summary.txt")))
  } else {
    writeLines("No climate PC retained by forward selection.",
               file.path(outdir, paste0(label, "_final_model_NONE_selected.txt")))
  }

  # ---------- B. ordiR2step cross-check ----------
  ordi_text <- tryCatch({
    ordi <- vegan::ordiR2step(
      null_model,
      scope = formula(full_model),
      R2scope = TRUE,
      Pin = alpha,
      permutations = nperm,
      trace = TRUE
    )
    c(capture.output(print(ordi$anova)),
      "",
      "Selected formula:",
      deparse(formula(ordi)))
  }, error = function(e) {
    paste("ordiR2step failed:", conditionMessage(e))
  })

  writeLines(ordi_text,
             file.path(outdir, paste0(label, "_ordiR2step_result.txt")))
  cat("\nordiR2step cross-check:\n")
  cat(ordi_text, sep = "\n")

  # ---------- Summary row ----------
  tibble(
    model = label,
    n_populations = nrow(response),
    conditioning = if (length(cond_names)) paste(cond_names, collapse = " + ") else "none",
    full_model_overall_P = as.numeric(as.data.frame(full_overall)[["Pr(>F)"]][1]),
    full_model_adjusted_R2 = full_r2$adj.r.squared,
    PCs_selected = if (length(selected)) paste(selected, collapse = " + ") else "none",
    final_model_overall_P = final_p,
    final_model_adjusted_R2 = final_r2$adjusted_R2,
    smallest_candidate_P_step1 = min(candidate_log[[1]]$P)
  )
}

# ============================================================
# 3. LOAD SAVED OBJECTS FROM THE MAIN ANALYSIS
# ============================================================

main_objects_path <- file.path(main_output_dir, "RDA_complete_analysis_objects.rds")
if (!file.exists(main_objects_path)) {
  stop("Could not find: ", main_objects_path,
       "\nRun 01_IBD_IBE_transcriptome_RDA.R first.")
}

main_objects     <- readRDS(main_objects_path)
rda_response     <- main_objects$rda_response
environment_data <- main_objects$environment_data
population_ids   <- rownames(rda_response)

environment_data <- environment_data[population_ids, , drop = FALSE]
stopifnot(identical(rownames(environment_data), population_ids))

required <- c(climate_candidates, salinity_predictor)
if (length(setdiff(required, colnames(environment_data))) > 0) {
  stop("Missing predictors in environment_data: ",
       paste(setdiff(required, colnames(environment_data)), collapse = ", "))
}

# Overall selected spatial dbMEMs (same ones used in the main partial RDA)
mem_matrix <- main_objects$overall_selected_mem_matrix
if (!is.null(mem_matrix) && ncol(mem_matrix) > 0) {
  mem_matrix <- as.data.frame(mem_matrix)
  mem_matrix <- mem_matrix[population_ids, , drop = FALSE]
  mem_names <- colnames(mem_matrix)
  model_data <- cbind(environment_data, mem_matrix)
} else {
  mem_names <- character(0)
  model_data <- environment_data
}
cat("Overall spatial dbMEMs used as conditions:",
    if (length(mem_names)) paste(mem_names, collapse = ", ") else "none", "\n")

# ============================================================
# 4. RUN FORWARD SELECTION (full dataset)
# ============================================================

res_full <- run_forward_selection(
  response   = rda_response,
  env        = model_data,
  candidates = climate_candidates,
  cond_names = c(salinity_predictor, mem_names),
  label      = "Full_dataset_climate_given_salinity_and_space",
  outdir     = step_output_dir,
  nperm      = n_permutations,
  alpha      = alpha_in
)

# ============================================================
# 5. SUMMARY
# ============================================================

summary_tbl <- res_full
write_csv(summary_tbl, file.path(step_output_dir, "stepwise_selection_summary.csv"))

cat("\n============================================================\n")
cat("STEPWISE SELECTION SUMMARY\n")
cat("============================================================\n")
print(as.data.frame(summary_tbl))

cat("\nNotes:\n")
cat("- *_forward_all_candidate_tests.csv lists the P value of EVERY PC at every step.\n")
cat("- Forward-selection P values are NOT corrected for multiple testing or for the\n")
cat("  selection procedure itself; treat a PC retained after selection as exploratory,\n")
cat("  and report the full model too.\n")
cat("Outputs in: ", step_output_dir, "\n", sep = "")
