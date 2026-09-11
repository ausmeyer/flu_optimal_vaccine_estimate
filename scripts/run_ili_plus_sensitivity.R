#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/analysis_config.R")
source("R/latent_curve_gam.R")
source("R/burden_posterior.R")
source("R/ili_plus.R")
source("R/decision_engine.R")
source("R/plotting.R")
source("R/ve_draws_from_estimates.R")
source("R/analysis_outputs.R")

n_draws <- as.integer(Sys.getenv("N_DRAWS", "5000"))
ili_weighting <- match.arg(
  Sys.getenv("ILI_PLUS_ILI_WEIGHTING", "cdc_population_weighted"),
  c("cdc_population_weighted", "national_unweighted")
)
lab_scope <- match.arg(
  Sys.getenv("ILI_PLUS_LAB_SCOPE", "all_available"),
  c("all_available", "clinical_only")
)
default_prefix <- paste(
  "empirical_ili_plus",
  if (lab_scope == "clinical_only") "clinical_only" else "all_seasons",
  if (ili_weighting == "cdc_population_weighted") "weighted" else "unweighted",
  sep = "_"
)
output_prefix <- Sys.getenv("OUTPUT_PREFIX", default_prefix)
output_dirs <- analysis_output_dirs(output_prefix, default_class = "sensitivity")
timing_sd <- timing_shift_sd()
ve_reference <- ve_reference_weeks()
candidate_weeks <- c(36:52, 1:12)

ili_plus_path <- "data/processed/national_ili_plus_primary_seasons.csv"
ilinet_path <- "data/processed/ilinet_national_all_age_primary_seasons.csv"
nrevss_path <- "data/processed/nrevss_clinical_labs_primary_seasons.csv"
ve_path <- "data/processed/cdc_ve_estimates.csv"

ili_plus <- read_csv(ili_plus_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("national ILI+ data")
analysis_data <- filter_ili_plus_lab_scope(ili_plus, lab_scope)
available_seasons <- sort(unique(analysis_data$season_start_year))
if (length(available_seasons) < 3L) {
  stop("ILI+ sensitivity requires at least three influenza seasons.")
}
ve_estimates <- read_csv(ve_path, show_col_types = FALSE)

burden_draws <- draw_us_empirical_ili_plus_burden(
  ili_plus = ili_plus,
  n_draws = n_draws,
  ili_weighting = ili_weighting,
  lab_scope = lab_scope,
  seed = 20260507,
  timing_shift_sd = timing_sd
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

regret_y_label <- paste0(
  "ILI+-weighted visit equivalents potentially averted\n",
  "by optimizing timing"
)
regret_caption <- paste0(
  "Each simulation is compared with its optimal week. Standardized to ",
  "10,000 outpatient encounters per surveillance week. ILI+ is an ",
  "influenza-specific proxy, not a count of laboratory-confirmed visits."
)
write_decision_outputs(
  results = results,
  ve_draws = ve_draws,
  output_prefix = output_prefix,
  burden_draws = burden_draws,
  write_burden_draws = env_flag("WRITE_BURDEN_DRAWS", FALSE),
  output_dir = output_dirs$tables
)
save_decision_figures(
  results = results,
  output_prefix = output_prefix,
  age_group = "all",
  output_dir = output_dirs$figures,
  regret_y_label = regret_y_label,
  regret_caption = regret_caption
)
write_analysis_manifest(
  output_prefix = output_prefix,
  input_paths = c(ili_plus_path, ilinet_path, nrevss_path, ve_path),
  output_dir = output_dirs$tables,
  configuration = list(
    n_draws = n_draws,
    primary_season_start_years = primary_season_start_years(),
    available_season_start_years = available_seasons,
    burden_scale =
      "national ILI proportion multiplied by national NREVSS positivity",
    outcome_interpretation =
      "influenza-specific outpatient-burden proxy; not confirmed visits",
    ili_weighting = ili_weighting,
    nrevss_weighting = "national laboratory-test-volume weighting",
    nrevss_lab_scope = lab_scope,
    observation_sampling =
      "independent Jeffreys-prior beta draws for ILINet and NREVSS",
    smoothing_spar = 0.65,
    ve_reference_weeks = ve_reference,
    waning_prior = Sys.getenv("WANING_PRIOR", "ray_2019"),
    immune_lag_prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w"),
    timing_shift_sd_weeks = timing_sd,
    figure_regret_y_label = regret_y_label,
    figure_regret_caption = regret_caption
  )
)

print(results$optimal_intervals %>% filter(.data$state == "US"))
print(results$week_cutoff_probabilities %>% filter(.data$state == "US"))
message("National ILI+ sensitivity analysis complete.")
