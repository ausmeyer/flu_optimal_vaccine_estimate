#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("R/calendar.R")

output_dir <- "outputs/diagnostics/monte_carlo"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

candidate_weeks <- c(36:52, 1:12)
optimal_paths <- list.files(
  c("outputs/primary", "outputs/sensitivity"),
  pattern = "_optimal_week_draws[.]csv$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(optimal_paths) == 0L) {
  stop("No current optimal-week draw files were found.")
}

summary_specs <- tibble(optimal_path = optimal_paths) %>%
  mutate(
    analysis = sub("_optimal_week_draws[.]csv$", "", basename(.data$optimal_path)),
    model = if_else(
      grepl("^empirical", .data$analysis),
      "empirical",
      "GAM"
    ),
    near_path = sub(
      "_optimal_week_draws[.]csv$",
      "_near_optimal_curve.csv",
      .data$optimal_path
    )
  ) %>%
  select("model", "analysis", "optimal_path", "near_path")

binomial_mcse <- function(p, n) {
  sqrt(pmax(p * (1 - p), 0) / n)
}

summarise_optimal_file <- function(model, analysis, optimal_path, near_path) {
  if (!file.exists(optimal_path)) {
    warning("Missing optimal draw file: ", optimal_path)
    return(list(probability = tibble(), split = tibble()))
  }

  draws <- read_csv(optimal_path, show_col_types = FALSE) %>%
    filter(.data$optimal_week %in% candidate_weeks)
  source_age_group <- paste(sort(unique(draws$age_group)), collapse = ",")

  n_by_group <- draws %>%
    distinct(.data$draw, .data$state, .data$age_group) %>%
    count(.data$state, .data$age_group, name = "n_draws")

  exact_summary <- draws %>%
    count(.data$state, .data$age_group, .data$optimal_week, name = "n") %>%
    tidyr::complete(
      state,
      age_group,
      optimal_week = candidate_weeks,
      fill = list(n = 0L)
    ) %>%
    left_join(n_by_group, by = c("state", "age_group")) %>%
    mutate(
      model = .env$model,
      analysis = .env$analysis,
      source_age_group = .env$source_age_group,
      probability = .data$n / .data$n_draws,
      mcse = binomial_mcse(.data$probability, .data$n_draws),
      mc95_half_width = 1.96 * .data$mcse,
      summary_type = "exact_argmax",
      week = .data$optimal_week,
      tolerance = NA_real_
    ) %>%
    select(
      "model",
      "analysis",
      "source_age_group",
      "state",
      "age_group",
      "summary_type",
      "week",
      "tolerance",
      "probability",
      "n_draws",
      "mcse",
      "mc95_half_width"
    )

  split_summary <- draws %>%
    mutate(split = if_else(.data$draw <= median(unique(.data$draw)), "first_half", "second_half")) %>%
    count(.data$state, .data$age_group, .data$split, .data$optimal_week, name = "n") %>%
    group_by(.data$state, .data$age_group, .data$split) %>%
    mutate(probability = .data$n / sum(.data$n)) %>%
    ungroup() %>%
    select("state", "age_group", "split", "optimal_week", "probability") %>%
    tidyr::pivot_wider(names_from = "split", values_from = "probability", values_fill = 0) %>%
    mutate(
      model = .env$model,
      analysis = .env$analysis,
      source_age_group = .env$source_age_group,
      abs_split_difference = abs(.data$first_half - .data$second_half)
    )

  split_out <- split_summary %>%
    group_by(
      .data$model, .data$analysis, .data$source_age_group,
      .data$state, .data$age_group
    ) %>%
    summarise(
      max_exact_argmax_split_difference = max(.data$abs_split_difference, na.rm = TRUE),
      mean_exact_argmax_split_difference = mean(.data$abs_split_difference, na.rm = TRUE),
      .groups = "drop"
    )

  near_summary <- tibble()
  if (file.exists(near_path)) {
    near_summary <- read_csv(near_path, show_col_types = FALSE) %>%
      filter(.data$vaccination_week %in% candidate_weeks) %>%
      left_join(n_by_group, by = c("state", "age_group")) %>%
      mutate(
        model = .env$model,
        analysis = .env$analysis,
        source_age_group = .env$source_age_group,
        probability = .data$probability_near_optimal,
        mcse = binomial_mcse(.data$probability, .data$n_draws),
        mc95_half_width = 1.96 * .data$mcse,
        summary_type = "near_optimal",
        week = .data$vaccination_week,
        tolerance = .data$relative_regret_tolerance
      ) %>%
      select(
        "model",
        "analysis",
        "source_age_group",
        "state",
        "age_group",
        "summary_type",
        "week",
        "tolerance",
        "probability",
        "n_draws",
        "mcse",
        "mc95_half_width"
      )
  }

  list(probability = bind_rows(exact_summary, near_summary), split = split_out)
}

audits <- pmap(summary_specs, summarise_optimal_file)
probability_audit <- map_dfr(audits, "probability")
split_audit <- map_dfr(audits, "split")

probability_rollup <- probability_audit %>%
  group_by(
    .data$model, .data$analysis, .data$source_age_group,
    .data$summary_type, .data$tolerance
  ) %>%
  summarise(
    max_mcse = max(.data$mcse, na.rm = TRUE),
    max_mc95_half_width = max(.data$mc95_half_width, na.rm = TRUE),
    median_mc95_half_width = stats::median(.data$mc95_half_width, na.rm = TRUE),
    .groups = "drop"
  )

split_rollup <- split_audit %>%
  group_by(.data$model, .data$analysis, .data$source_age_group) %>%
  summarise(
    max_state_split_difference = max(.data$max_exact_argmax_split_difference, na.rm = TRUE),
    median_state_split_difference = stats::median(.data$max_exact_argmax_split_difference, na.rm = TRUE),
    us_split_difference = .data$max_exact_argmax_split_difference[.data$state == "US"][[1]],
    .groups = "drop"
  )

write_csv(probability_audit, file.path(output_dir, "monte_carlo_probability_precision.csv"))
write_csv(probability_rollup, file.path(output_dir, "monte_carlo_probability_precision_rollup.csv"))
write_csv(split_audit, file.path(output_dir, "monte_carlo_exact_argmax_split_audit.csv"))
write_csv(split_rollup, file.path(output_dir, "monte_carlo_exact_argmax_split_rollup.csv"))

cat("Wrote Monte Carlo precision audits to ", output_dir, "\n", sep = "")
print(probability_rollup %>% filter(.data$model == "GAM"), n = Inf)
print(split_rollup %>% filter(.data$model == "GAM"), n = Inf)
