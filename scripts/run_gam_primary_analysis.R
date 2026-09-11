#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

source("R/analysis_config.R")
source("R/latent_curve_gam.R")
source("R/outcome_time.R")
source("R/decision_engine.R")
source("R/plotting.R")
source("R/ve_draws_from_estimates.R")
source("R/analysis_outputs.R")

n_draws <- as.integer(Sys.getenv("N_DRAWS", "1000"))
output_prefix <- Sys.getenv("OUTPUT_PREFIX", "gam_primary")
output_dirs <- analysis_output_dirs(output_prefix, default_class = "primary")
reuse_gam_models <- env_flag("REUSE_GAM_MODELS", TRUE)
write_burden_draws <- env_flag("WRITE_BURDEN_DRAWS", FALSE)
timing_sd <- timing_shift_sd()
ve_reference <- ve_reference_weeks()
waning_prior <- Sys.getenv("WANING_PRIOR", "ray_2019")
outcome_date_prior <- match.arg(
  Sys.getenv("OUTCOME_DATE_PRIOR", "none"),
  c("none", "flu_ve_outpatient_onset_to_enrollment")
)
if (outcome_date_prior != "none" &&
    waning_prior != "spencer_ferdinands_fast") {
  stop(
    "OUTCOME_DATE_PRIOR=", outcome_date_prior,
    " is only defined for WANING_PRIOR=spencer_ferdinands_fast. ",
    "The Ray waning clock is indexed to PCR test date."
  )
}
candidate_weeks <- c(36:52, 1:12)

state_path <- "data/processed/ilinet_state_all_age_primary_seasons.csv"
population_path <- "data/processed/state_population.csv"
ve_path <- "data/processed/cdc_ve_estimates.csv"

dir.create("models", recursive = TRUE, showWarnings = FALSE)

ilinet_state <- read_csv(state_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("ILINet all-age state data")
state_population <- read_state_population(
  population_path,
  expected_states = unique(ilinet_state$state)
)
ve_estimates <- read_csv(ve_path, show_col_types = FALSE)

model_data <- prepare_all_age_gam_data(ilinet_state)
model_path <- model_path_with_seasons(
  paste0("gam_all_age_ilinet_statewise_quasibinomial_wk36", gam_smoothing_slug()),
  model_data,
  provenance_paths = c(state_path, "R/latent_curve_gam.R")
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
  file.path(output_dirs$tables, paste0(output_prefix, "_ilinet_statewise_diagnostics.csv"))
)

burden_draws <- draw_statewise_gam_decision_burden(
  state_models = state_models,
  state_population = state_population,
  n_draws = n_draws,
  seed = 20260507,
  timing_shift_sd = timing_sd
)
if (outcome_date_prior != "none") {
  burden_draws <- apply_outcome_date_alignment(
    burden_draws,
    prior = outcome_date_prior,
    seed = 20260531
  )
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

write_decision_outputs(
  results = results,
  ve_draws = ve_draws,
  output_prefix = output_prefix,
  burden_draws = burden_draws,
  write_burden_draws = write_burden_draws,
  output_dir = output_dirs$tables
)
save_decision_figures(
  results,
  output_prefix,
  age_group = "all",
  output_dir = output_dirs$figures
)
write_analysis_manifest(
  output_prefix = output_prefix,
  input_paths = c(
    state_path,
    population_path,
    ve_path,
    "scripts/run_gam_primary_analysis.R",
    "R/latent_curve_gam.R",
    "R/decision_engine.R",
    if (outcome_date_prior != "none") "R/outcome_time.R"
  ),
  model_path = model_path,
  output_dir = output_dirs$tables,
  configuration = list(
    n_draws = n_draws,
    primary_season_start_years = primary_season_start_years(),
    burden_scale = "state latent ILI proportion; US population-weighted mean",
    missing_state_week_handling = "paired missing outpatient counts retained only on the prediction grid with zero fitting weight",
    ve_reference_weeks = ve_reference,
    waning_prior = waning_prior,
    immune_lag_prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w"),
    timing_shift_sd_weeks = timing_sd,
    outcome_date_prior = outcome_date_prior,
    outcome_date_source = if (outcome_date_prior == "none") {
      "none"
    } else {
      "Balasubramani et al. 2020, doi:10.1111/irv.12741"
    },
    outcome_date_within_interval_assumption = if (outcome_date_prior == "none") {
      "none"
    } else {
      "discrete uniform integer day within 0-2, 3-4, and 5-7 day categories"
    },
    outcome_date_category_counts = if (outcome_date_prior == "none") {
      "none"
    } else {
      "0-2 days:8520; 3-4 days:10754; 5-7 days:7861"
    },
    outcome_date_seed = if (outcome_date_prior == "none") {
      NA_integer_
    } else {
      20260531L
    },
    gam_global_k = as.integer(Sys.getenv("GAM_GLOBAL_K", "14")),
    gam_season_k = as.integer(Sys.getenv("GAM_SEASON_K", "12")),
    gam_gamma = as.numeric(Sys.getenv("GAM_GAMMA", "1"))
  )
)

print(results$optimal_intervals %>% filter(.data$state == "US"))
print(results$week_cutoff_probabilities %>% filter(.data$state == "US"))
message("GAM latent-curve primary analysis complete.")
