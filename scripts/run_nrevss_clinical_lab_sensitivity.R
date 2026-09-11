#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

source("R/analysis_config.R")
source("R/latent_curve_gam.R")
source("R/decision_engine.R")
source("R/plotting.R")
source("R/ve_draws_from_estimates.R")
source("R/analysis_outputs.R")

n_draws <- as.integer(Sys.getenv("N_DRAWS", "1000"))
analysis_region <- Sys.getenv("NREVSS_ANALYSIS_REGION", "national")
target_virus <- match.arg(Sys.getenv("NREVSS_TARGET", "all"), c("all", "a", "b"))
output_prefix <- Sys.getenv(
  "OUTPUT_PREFIX",
  paste0(
    "gam_nrevss_clinical_", analysis_region,
    if (target_virus == "all") "" else paste0("_flu_", target_virus)
  )
)
output_dirs <- analysis_output_dirs(output_prefix, default_class = "sensitivity")
reuse_gam_models <- env_flag("REUSE_GAM_MODELS", TRUE)
timing_sd <- timing_shift_sd()
ve_reference <- ve_reference_weeks()
candidate_weeks <- c(36:52, 1:12)

nrevss_path <- "data/processed/nrevss_clinical_labs_primary_seasons.csv"
ve_path <- "data/processed/cdc_ve_estimates.csv"

dir.create("models", recursive = TRUE, showWarnings = FALSE)

nrevss <- read_csv(nrevss_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("NREVSS clinical laboratory data") %>%
  filter(.data$source_region_type == analysis_region)
ve_estimates <- read_csv(ve_path, show_col_types = FALSE)
if (nrow(nrevss) == 0L) {
  stop("No NREVSS records found for NREVSS_ANALYSIS_REGION=", analysis_region)
}

if (target_virus != "all") {
  target_col <- paste0("positive_tests_", target_virus)
  if (!target_col %in% names(nrevss)) {
    stop("NREVSS_TARGET=", target_virus, " requires column ", target_col, ".")
  }
  nrevss <- nrevss %>%
    mutate(
      ilitotal = .data[[target_col]],
      unweighted_ili = 100 * .data$ilitotal / pmax(.data$total_tests, 1)
    ) %>%
    filter(!is.na(.data$ilitotal))
}

model_data <- prepare_all_age_gam_data(nrevss)
model_path <- model_path_with_seasons(
  paste0(
    "gam_nrevss_clinical_", gsub("[^A-Za-z0-9]+", "_", analysis_region),
    "_", target_virus, "_wk36", gam_smoothing_slug()
  ),
  model_data,
  provenance_paths = c(nrevss_path, "R/latent_curve_gam.R")
)
if (reuse_gam_models && file.exists(model_path)) {
  message("Loading existing GAM model: ", model_path)
  state_models <- readRDS(model_path)
} else {
  state_models <- fit_statewise_latent_gams(model_data, family = "quasibinomial")
  saveRDS(state_models, model_path)
}

diagnostics <- imap_dfr(state_models, ~ gam_diagnostic_row(.x, .y))
write_csv(
  diagnostics,
  file.path(output_dirs$tables, paste0(output_prefix, "_diagnostics.csv"))
)

burden_draws <- draw_statewise_gam_decision_burden(
  state_models = state_models,
  n_draws = n_draws,
  seed = 20260507,
  timing_shift_sd = timing_sd
)
if (analysis_region == "national") {
  burden_draws <- burden_draws %>% mutate(state = "US")
}
ve_draws <- draw_primary_ve_waning(
  ve_estimates = ve_estimates,
  n_draws = n_draws,
  age_group = "all",
  reference_weeks = ve_reference,
  seed = 20260509
)
candidate_utilities <- evaluate_candidate_utilities(
  burden_draws = burden_draws,
  ve_draws = ve_draws,
  candidate_weeks = candidate_weeks
)
results <- summarise_decision_outputs(candidate_utilities)

write_decision_outputs(results, ve_draws, output_prefix, output_dir = output_dirs$tables)
save_decision_figures(
  results,
  output_prefix,
  age_group = "all",
  output_dir = output_dirs$figures
)
write_analysis_manifest(
  output_prefix = output_prefix,
  input_paths = c(nrevss_path, ve_path),
  model_path = model_path,
  output_dir = output_dirs$tables,
  configuration = list(
    n_draws = n_draws,
    primary_season_start_years = primary_season_start_years(),
    analysis_region = analysis_region,
    target_virus = target_virus,
    burden_scale = "latent clinical-laboratory influenza positivity proportion",
    ve_reference_weeks = ve_reference,
    waning_prior = Sys.getenv("WANING_PRIOR", "ray_2019"),
    immune_lag_prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w"),
    timing_shift_sd_weeks = timing_sd
  )
)

print(results$optimal_intervals %>% filter(.data$state == "US"))
print(results$week_cutoff_probabilities %>% filter(.data$state == "US"))
message("GAM NREVSS clinical-laboratory sensitivity complete.")
