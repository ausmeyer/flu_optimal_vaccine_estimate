#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("R/calendar.R")

input_dirs <- c("outputs/primary", "outputs/sensitivity")
output_dir <- "outputs/sensitivity/summary/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_with_analysis <- function(path, suffix) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(analysis = sub(paste0(suffix, "$"), "", basename(path)), .before = 1)
}

interval_files <- list.files(
  input_dirs,
  pattern = "_optimal_week_intervals\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
interval_files <- interval_files[!grepl("_specific_optimal_week_intervals\\.csv$", interval_files)]
cutoff_files <- list.files(
  input_dirs,
  pattern = "_week_cutoff_probabilities\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
regret_files <- list.files(
  input_dirs,
  pattern = "_regret_curve\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
near_files <- list.files(
  input_dirs,
  pattern = "_near_optimal_curve\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(interval_files) == 0L) {
  stop("No optimal-week interval files found.")
}

intervals <- map_dfr(interval_files, read_with_analysis, suffix = "_optimal_week_intervals\\.csv") %>%
  filter(.data$state == "US") %>%
  select(
    "analysis",
    "state",
    "age_group",
    "median_week",
    "lower_50_week",
    "upper_50_week",
    "lower_95_week",
    "upper_95_week"
  )

cutoff_probs <- if (length(cutoff_files) > 0L) {
  map_dfr(
    cutoff_files,
    read_with_analysis,
    suffix = "_week_cutoff_probabilities\\.csv"
  ) %>%
    filter(.data$state == "US") %>%
    select("analysis", "state", "age_group", "cutoff_week", "probability") %>%
    pivot_wider(
      names_from = "cutoff_week",
      values_from = "probability",
      names_prefix = "p_optimal_by_mmwr_week_"
    )
} else {
  tibble()
}

regret_minima <- if (length(regret_files) > 0L) {
  map_dfr(regret_files, read_with_analysis, suffix = "_regret_curve\\.csv") %>%
    filter(.data$state == "US") %>%
    group_by(.data$analysis, .data$state, .data$age_group) %>%
    slice_min(.data$mean_regret, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      analysis = .data$analysis,
      state = .data$state,
      age_group = .data$age_group,
      min_mean_regret_week = .data$vaccination_week,
      min_mean_regret = .data$mean_regret
    )
} else {
  tibble()
}

near_5pct_top <- if (length(near_files) > 0L) {
  map_dfr(near_files, read_with_analysis, suffix = "_near_optimal_curve\\.csv") %>%
    filter(.data$state == "US", .data$relative_regret_tolerance == 0.05) %>%
    group_by(.data$analysis, .data$state, .data$age_group) %>%
    slice_max(.data$probability_near_optimal, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      analysis = .data$analysis,
      state = .data$state,
      age_group = .data$age_group,
      top_near_optimal_5pct_week = .data$vaccination_week,
      top_near_optimal_5pct_probability = .data$probability_near_optimal
    )
} else {
  tibble()
}

summary <- intervals %>%
  left_join(cutoff_probs, by = c("analysis", "state", "age_group")) %>%
  left_join(regret_minima, by = c("analysis", "state", "age_group")) %>%
  left_join(near_5pct_top, by = c("analysis", "state", "age_group")) %>%
  arrange(.data$analysis, factor(.data$age_group, levels = c("all", "0-4", "5-24", "25-49", "50-64", "65+")))

loso_summary <- summary %>%
  filter(grepl("^gam_loso_drop_", .data$analysis)) %>%
  mutate(dropped_season_start_year = as.integer(sub("^gam_loso_drop_", "", .data$analysis))) %>%
  select("dropped_season_start_year", everything())

write_csv(summary, file.path(output_dir, "robustness_suite_us_summary.csv"))
if (nrow(loso_summary) > 0L) {
  write_csv(loso_summary, file.path(output_dir, "gam_leave_one_season_out_us_summary.csv"))
}

cat("Wrote robustness summary tables to ", output_dir, "\n", sep = "")
print(summary %>% filter(.data$state == "US"), n = Inf)
