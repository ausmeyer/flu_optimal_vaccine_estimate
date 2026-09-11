#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

source("R/analysis_config.R")
source("R/hhs_regions.R")
source("R/latent_curve_gam.R")
source("R/decision_engine.R")
source("R/plotting.R")
source("R/ve_draws_from_estimates.R")
source("R/analysis_outputs.R")

n_draws <- as.integer(Sys.getenv("N_DRAWS", "1000"))
output_prefix <- Sys.getenv("OUTPUT_PREFIX", "gam_age")
output_dirs <- analysis_output_dirs(output_prefix, default_class = "primary")
reuse_gam_models <- env_flag("REUSE_GAM_MODELS", TRUE)
write_burden_draws <- env_flag("WRITE_BURDEN_DRAWS", FALSE)
timing_sd <- timing_shift_sd()
ve_reference <- ve_reference_weeks()
candidate_weeks <- c(36:52, 1:12)

age_map <- tibble::tribble(
  ~age_group, ~ve_age_group,
  "0-4", "6m-8",
  "5-24", "9-17",
  "25-49", "18-49",
  "50-64", "50-64",
  "65+", "65+"
)
requested <- trimws(strsplit(Sys.getenv("AGE_GROUPS", ""), ",", fixed = TRUE)[[1]])
requested <- requested[nzchar(requested)]
if (length(requested) > 0L) {
  unknown <- setdiff(requested, age_map$age_group)
  if (length(unknown) > 0L) {
    stop("Unknown AGE_GROUPS value(s): ", paste(unknown, collapse = ", "))
  }
  age_map <- age_map %>% filter(.data$age_group %in% requested)
}

state_path <- "data/processed/ilinet_state_all_age_primary_seasons.csv"
hhs_age_path <- "data/processed/ilinet_hhs_age_primary_seasons.csv"
population_path <- "data/processed/state_population.csv"
ve_path <- "data/processed/cdc_ve_estimates.csv"

dir.create("models", recursive = TRUE, showWarnings = FALSE)

ilinet_state <- read_csv(state_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("ILINet all-age state data")
ilinet_hhs_age <- read_csv(hhs_age_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("ILINet HHS age data")
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
  message("Loading existing all-age GAM model: ", model_path)
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

run_one_age <- function(age_group, ve_age_group) {
  age_index <- match(age_group, c("0-4", "5-24", "25-49", "50-64", "65+"))
  slug <- gsub("[^A-Za-z0-9]+", "_", age_group)
  file_prefix <- paste0(output_prefix, "_", slug)
  message("Running age group ", age_group, " with VE source stratum ", ve_age_group)

  burden_draws <- draw_statewise_gam_age_decision_burden(
    state_models = state_models,
    ilinet_hhs_age = ilinet_hhs_age,
    hhs_regions = hhs_region_crosswalk(),
    state_population = state_population,
    target_age_group = age_group,
    n_draws = n_draws,
    seed = 20260507 + age_index,
    timing_shift_sd = timing_sd
  )
  ve_draws <- draw_primary_ve_waning(
    ve_estimates = ve_estimates,
    n_draws = n_draws,
    age_group = age_group,
    ve_age_group = ve_age_group,
    reference_weeks = ve_reference,
    seed = 20260520 + age_index
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
    output_prefix = file_prefix,
    burden_draws = burden_draws,
    write_burden_draws = write_burden_draws,
    output_dir = output_dirs$tables
  )
  save_decision_figures(results, file_prefix, age_group, output_dir = output_dirs$figures)
  write_analysis_manifest(
    output_prefix = file_prefix,
    input_paths = c(
      state_path,
      hhs_age_path,
      population_path,
      ve_path,
      "scripts/run_gam_age_specific_analysis.R",
      "R/hhs_regions.R",
      "R/latent_curve_gam.R",
      "R/decision_engine.R"
    ),
    model_path = model_path,
    output_dir = output_dirs$tables,
    configuration = list(
      n_draws = n_draws,
      age_group = age_group,
      ve_age_group = ve_age_group,
      primary_season_start_years = primary_season_start_years(),
      age_model = "all-age state GAM multiplied by sampled HHS ILI age composition",
      missing_state_week_handling = "paired missing outpatient counts retained only on the prediction grid with zero fitting weight",
      hhs_age_dirichlet_offset = 0.5,
      national_aggregation = "state population-weighted mean",
      ve_reference_weeks = ve_reference,
      waning_prior = Sys.getenv("WANING_PRIOR", "ray_2019"),
      immune_lag_prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w"),
      timing_shift_sd_weeks = timing_sd,
      gam_global_k = as.integer(Sys.getenv("GAM_GLOBAL_K", "14")),
      gam_season_k = as.integer(Sys.getenv("GAM_SEASON_K", "12")),
      gam_gamma = as.numeric(Sys.getenv("GAM_GAMMA", "1"))
    )
  )

  intervals <- results$optimal_intervals
  rm(burden_draws, ve_draws, candidate_utilities, results)
  gc()
  intervals
}

intervals <- pmap_dfr(age_map, run_one_age)
write_csv(
  intervals,
  file.path(output_dirs$tables, paste0(output_prefix, "_specific_optimal_week_intervals.csv"))
)
print(intervals %>% filter(.data$state == "US"))
message("GAM latent-curve age-specific analyses complete.")
