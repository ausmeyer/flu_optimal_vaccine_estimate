#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/analysis_config.R")
source("R/burden_posterior.R")
source("R/decision_engine.R")
source("R/plotting.R")
source("R/ve_draws_from_estimates.R")
source("R/analysis_outputs.R")

n_draws <- as.integer(Sys.getenv("N_DRAWS", "2000"))
output_prefix <- Sys.getenv("OUTPUT_PREFIX", "empirical_primary")
output_dirs <- analysis_output_dirs(output_prefix, default_class = "sensitivity")
timing_sd <- timing_shift_sd()
ve_reference <- ve_reference_weeks()
candidate_weeks <- c(36:52, 1:12)
empirical_curve_smoothing <- "none"

national_path <- "data/processed/ilinet_national_all_age_primary_seasons.csv"
ve_path <- "data/processed/cdc_ve_estimates.csv"

ilinet_national <- read_csv(national_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("ILINet all-age national data")
ve_estimates <- read_csv(ve_path, show_col_types = FALSE)

burden_draws <- draw_us_empirical_ilinet_burden(
  ilinet_national = ilinet_national,
  n_draws = n_draws,
  seed = 20260508,
  include_observation_noise = TRUE,
  timing_shift_sd = timing_sd,
  curve_smoothing = empirical_curve_smoothing,
  sample_complete_seasons_only = TRUE
)

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
  input_paths = c(
    national_path,
    ve_path,
    "scripts/run_primary_analysis.R",
    "R/burden_posterior.R",
    "R/decision_engine.R"
  ),
  output_dir = output_dirs$tables,
  configuration = list(
    n_draws = n_draws,
    primary_season_start_years = primary_season_start_years(),
    analysis_scope = "national all-age ILINet",
    national_burden_scale = "unsmoothed CDC population-weighted national ILI proportion with beta observation draws",
    empirical_curve_smoothing = empirical_curve_smoothing,
    include_observation_noise = TRUE,
    observation_sampling = "Jeffreys-prior beta draws from weekly ILI and total outpatient visit counts",
    national_weighting = "beta draw from national unweighted counts rescaled by the observed CDC weighted-to-unweighted ILI ratio",
    season_sampling = "with replacement among eligible seasons extending through MMWR week 22",
    timing_shift_interpolation = "linear",
    timing_shift_boundary_rule = "hold the nearest observed boundary value constant",
    ve_reference_weeks = ve_reference,
    waning_prior = Sys.getenv("WANING_PRIOR", "ray_2019"),
    immune_lag_prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w"),
    timing_shift_sd_weeks = timing_sd
  )
)

print(results$optimal_intervals %>% filter(.data$state == "US"))
print(results$week_cutoff_probabilities %>% filter(.data$state == "US"))
message("Empirical ILINet sensitivity analysis complete.")
