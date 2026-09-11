#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("R/analysis_config.R")
source("R/decision_engine.R")
source("R/ve_draws_from_estimates.R")
source("R/analysis_outputs.R")

n_draws <- as.integer(Sys.getenv("N_DRAWS", "5000"))
candidate_weeks <- c(36:52, 1:12)
state_path <- "data/processed/ilinet_state_all_age_primary_seasons.csv"
population_path <- "data/processed/state_population.csv"
ve_path <- "data/processed/cdc_ve_estimates.csv"
sensitivity_root <- Sys.getenv("SENSITIVITY_ROOT", "outputs/sensitivity")
output_dir <- Sys.getenv(
  "VALIDATION_OUTPUT_DIR",
  "outputs/diagnostics/held_out_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ilinet <- read_csv(state_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("ILINet all-age state data")
population <- read_state_population(
  population_path,
  expected_states = unique(ilinet$state)
)
ve_estimates <- read_csv(ve_path, show_col_types = FALSE)

complete_seasons <- ilinet %>%
  group_by(.data$season, .data$season_start_year) %>%
  summarise(complete = any(.data$mmwr_week == 22L), .groups = "drop") %>%
  filter(.data$complete) %>%
  arrange(.data$season_start_year)

read_training_results <- function(omitted_year) {
  prefix <- paste0("gam_loso_drop_", omitted_year)
  table_dir <- file.path(sensitivity_root, prefix, "tables")
  required <- file.path(
    table_dir,
    paste0(
      prefix,
      c(
        "_regret_curve.csv",
        "_optimal_week_draws.csv",
        "_optimal_week_intervals.csv"
      )
    )
  )
  missing <- required[!file.exists(required)]
  if (length(missing) > 0L) {
    stop(
      "Held-out validation requires completed leave-one-season-out results. Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  list(
    regret = read_csv(required[[1]], show_col_types = FALSE) %>%
      filter(.data$state == "US", .data$age_group == "all"),
    optimal = read_csv(required[[2]], show_col_types = FALSE) %>%
      filter(.data$state == "US", .data$age_group == "all"),
    interval = read_csv(required[[3]], show_col_types = FALSE) %>%
      filter(.data$state == "US", .data$age_group == "all")
  )
}

build_held_out_burden <- function(held_out_season) {
  season_data <- ilinet %>%
    filter(.data$season == held_out_season) %>%
    mutate(burden = .data$ilitotal / .data$total_patients) %>%
    left_join(population, by = "state", relationship = "many-to-one")

  max_week <- if (any(season_data$mmwr_week == 53L)) 53L else 52L
  national <- season_data %>%
    filter(is_in_burden_window(.data$mmwr_week)) %>%
    group_by(.data$mmwr_week) %>%
    summarise(
      burden = weighted.mean(.data$burden, .data$population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(
      state = "US",
      age_group = "all",
      week = .data$mmwr_week,
      burden = .data$burden,
      max_week = max_week
    )

  crossing(draw = seq_len(n_draws), national)
}

score_one_season <- function(held_out_season, omitted_year) {
  message("Evaluating held-out season ", held_out_season)
  training <- read_training_results(omitted_year)
  burden_draws <- build_held_out_burden(held_out_season)
  ve_draws <- draw_primary_ve_waning(
    ve_estimates = ve_estimates,
    n_draws = n_draws,
    age_group = "all",
    reference_weeks = 8,
    seed = 20260509
  )
  utilities <- evaluate_candidate_utilities(
    burden_draws = burden_draws,
    ve_draws = ve_draws,
    candidate_weeks = candidate_weeks
  )
  held_out <- summarise_decision_outputs(utilities)

  training_week <- training$regret %>%
    slice_min(.data$mean_regret, n = 1, with_ties = FALSE) %>%
    pull(.data$vaccination_week)
  held_out_week <- held_out$regret_curve %>%
    filter(.data$state == "US", .data$age_group == "all") %>%
    slice_min(.data$mean_regret, n = 1, with_ties = FALSE) %>%
    pull(.data$vaccination_week)

  held_out_draws <- held_out$optimal_weeks %>%
    filter(.data$state == "US", .data$age_group == "all")
  held_out_interval <- held_out$optimal_intervals %>%
    filter(.data$state == "US", .data$age_group == "all")
  held_out_cutoffs <- held_out$week_cutoff_probabilities %>%
    filter(.data$state == "US", .data$age_group == "all") %>%
    select("cutoff_week", "probability") %>%
    pivot_wider(
      names_from = "cutoff_week",
      values_from = "probability",
      names_prefix = "p_by_week_"
    )

  training_counts <- training$optimal %>%
    count(.data$optimal_week, name = "n") %>%
    complete(optimal_week = candidate_weeks, fill = list(n = 0L)) %>%
    mutate(
      probability = (.data$n + 0.5) /
        (sum(.data$n) + 0.5 * length(candidate_weeks))
    )
  held_out_probabilities <- held_out_draws %>%
    count(.data$optimal_week, name = "n") %>%
    complete(optimal_week = candidate_weeks, fill = list(n = 0L)) %>%
    mutate(probability = .data$n / sum(.data$n))
  probability_comparison <- training_counts %>%
    select("optimal_week", training_probability = "probability") %>%
    left_join(
      held_out_probabilities %>%
        select("optimal_week", held_out_probability = "probability"),
      by = "optimal_week"
    )

  predictive_interval <- training$interval
  held_out_indices <- held_out_draws$optimal_season_index
  training_week_regret <- utilities %>%
    group_by(.data$draw) %>%
    mutate(best_utility = max(.data$utility)) %>%
    filter(.data$vaccination_week == training_week) %>%
    ungroup() %>%
    summarise(
      mean_regret = mean(.data$best_utility - .data$utility),
      mean_relative_regret = mean(
        (.data$best_utility - .data$utility) / .data$best_utility
      ),
      .groups = "drop"
    )

  tibble(
    held_out_season = held_out_season,
    held_out_season_start_year = omitted_year,
    training_regret_minimizing_week = training_week,
    held_out_regret_minimizing_week = held_out_week,
    held_out_median_optimal_week = held_out_interval$median_week,
    held_out_optimal_in_training_50_interval = mean(
      held_out_indices >= predictive_interval$lower_50_index &
        held_out_indices <= predictive_interval$upper_50_index
    ),
    held_out_optimal_in_training_95_interval = mean(
      held_out_indices >= predictive_interval$lower_95_index &
        held_out_indices <= predictive_interval$upper_95_index
    ),
    held_out_regret_at_training_week = training_week_regret$mean_regret,
    held_out_visits_lost_per_10000_at_training_week =
      10000 * training_week_regret$mean_regret,
    held_out_relative_regret_at_training_week =
      training_week_regret$mean_relative_regret,
    optimal_week_log_score = -sum(
      probability_comparison$held_out_probability *
        log(probability_comparison$training_probability)
    ),
    optimal_week_brier_score = sum(
      (probability_comparison$training_probability -
         probability_comparison$held_out_probability)^2
    )
  ) %>%
    bind_cols(held_out_cutoffs)
}

validation <- map2_dfr(
  complete_seasons$season,
  complete_seasons$season_start_year,
  score_one_season
)

summary <- validation %>%
  summarise(
    n_held_out_seasons = n(),
    median_visits_lost_per_10000 = median(
      .data$held_out_visits_lost_per_10000_at_training_week
    ),
    maximum_visits_lost_per_10000 = max(
      .data$held_out_visits_lost_per_10000_at_training_week
    ),
    median_relative_regret = median(
      .data$held_out_relative_regret_at_training_week
    ),
    mean_50_interval_coverage = mean(
      .data$held_out_optimal_in_training_50_interval
    ),
    mean_95_interval_coverage = mean(
      .data$held_out_optimal_in_training_95_interval
    ),
    mean_optimal_week_log_score = mean(.data$optimal_week_log_score),
    mean_optimal_week_brier_score = mean(.data$optimal_week_brier_score)
  )

write_csv(
  validation,
  file.path(output_dir, "held_out_season_decision_validation.csv")
)
write_csv(
  summary,
  file.path(output_dir, "held_out_season_decision_validation_summary.csv")
)

print(validation, n = Inf)
print(summary)
message("Held-out-season decision validation complete.")
