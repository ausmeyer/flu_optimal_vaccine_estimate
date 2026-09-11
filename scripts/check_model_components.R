#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/analysis_config.R")
source("R/hhs_regions.R")
source("R/latent_curve_gam.R")
source("R/burden_posterior.R")
source("R/ili_plus.R")
source("R/outcome_time.R")
source("R/decision_engine.R")
source("R/ve_draws_from_estimates.R")
source("R/publication_spec.R")
source("R/fluview_download.R")

# Calendar fixtures distinguish analytic week 36 from CDC download week 40,
# including the extra MMWR week in 2014.
stopifnot(identical(fluview_download_seasons(c(2010L, 2023L)),
  c(2009L, 2010L, 2022L, 2023L)))
fixture_dates <- seq(as.Date("2014-08-31"), as.Date("2015-05-31"), by = "week")
fixture_mmwr <- MMWRweek::MMWRweek(fixture_dates)
fixture_weeks <- data.frame(year = fixture_mmwr$MMWRyear, week = fixture_mmwr$MMWRweek)
stopifnot(
  nrow(fixture_weeks) == 40L,
  53L %in% fixture_weeks$week,
  fluview_raw_covers_burden_window(fixture_weeks, 2014L),
  !fluview_raw_covers_burden_window(fixture_weeks[-(1:4), ], 2014L),
  !fluview_raw_covers_burden_window(fixture_weeks[fixture_weeks$week != 53L, ], 2014L)
)

# Independent polynomial values in percentage points at elapsed weekly times.
# These checks detect use of weeks in place of the two-week polynomial unit.
stopifnot(
  max(abs(55 * spencer_relative_ve_retention(c(0, 8, 16, 26, 27), "fast") -
    c(55, 50.48, 40.20, 1.70, 0))) < 1e-10,
  max(abs(55 * spencer_relative_ve_retention(c(0, 8, 16, 36, 37), "slow") -
    c(55, 53.16, 49.08, 3.88, 0))) < 1e-10
)

state_data <- read_csv(
  "data/processed/ilinet_state_all_age_primary_seasons.csv",
  show_col_types = FALSE
)
hhs_age_data <- read_csv(
  "data/processed/ilinet_hhs_age_primary_seasons.csv",
  show_col_types = FALSE
)
national_data <- read_csv(
  "data/processed/ilinet_national_all_age_primary_seasons.csv",
  show_col_types = FALSE
)
nrevss_data <- read_csv(
  "data/processed/nrevss_clinical_labs_primary_seasons.csv",
  show_col_types = FALSE
)
ve_estimates <- read_csv(
  "data/processed/cdc_ve_estimates.csv",
  show_col_types = FALSE
)

# All-age and age-specific estimates intentionally have unequal season coverage.
# Give each synthetic source a distinct ratio to detect cross-age sampling.
ve_pool_fixture <- publication_ve_pool_coverage() %>%
  mutate(
    season = paste0(.data$season_start_year, "/",
      sprintf("%02d", (.data$season_start_year + 1L) %% 100L)),
    source_id = paste(.data$ve_age_group, .data$season, sep = ":"),
    log_or_mean = log(seq_len(n()) / 100),
    log_or_sd = 0
  )
ve_pool_sizes <- table(ve_pool_fixture$ve_age_group)
stopifnot(
  nrow(ve_pool_fixture) == 62L,
  ve_pool_sizes[["all"]] == 12L,
  all(ve_pool_sizes[names(ve_pool_sizes) != "all"] == 10L)
)
for (source_age_group in unique(ve_pool_fixture$ve_age_group)) {
  expected_pool <- filter(ve_pool_fixture, .data$ve_age_group == source_age_group)
  pool_draws <- draw_reported_ve_from_estimates(
    ve_pool_fixture,
    n_draws = 1000,
    age_group = "analysis_age",
    ve_age_group = source_age_group,
    seed = 98
  )
  source_index <- match(pool_draws$ve_source, expected_pool$source_id)
  stopifnot(
    !anyNA(source_index),
    setequal(pool_draws$ve_source, expected_pool$source_id),
    all(pool_draws$age_group == "analysis_age"),
    identical(pool_draws$ve_source_season, expected_pool$season[source_index]),
    max(abs(pool_draws$reported_effect_ratio -
      exp(expected_pool$log_or_mean[source_index]))) < 1e-12
  )
}

endpoint <- latest_observed_week(state_data)
stopifnot(endpoint$mmwr_year[[1]] == 2026L, endpoint$mmwr_week[[1]] == 22L)

latest_season <- 2025L
state_latest <- state_data %>%
  filter(.data$season_start_year == latest_season)
hhs_latest <- hhs_age_data %>%
  filter(.data$season_start_year == latest_season)
national_latest <- national_data %>%
  filter(.data$season_start_year == latest_season)
stopifnot(
  n_distinct(state_latest$state[state_latest$mmwr_week == 22L]) == 51L,
  n_distinct(hhs_latest$hhs_region[hhs_latest$mmwr_week == 22L]) == 10L,
  sum(national_latest$mmwr_week == 22L) == 1L
)

ili_plus_data <- prepare_national_ili_plus(national_data, nrevss_data)
clinical_only_ili_plus <- filter_ili_plus_lab_scope(
  ili_plus_data,
  "clinical_only"
)
expected_clinical_seasons <- ili_plus_data %>%
  summarise(
    eligible = all(.data$nrevss_lab_series == "clinical_labs") &&
      any(.data$mmwr_week == 36L) &&
      any(.data$mmwr_week == 22L),
    .by = c("season", "season_start_year")
  ) %>%
  filter(.data$eligible) %>%
  pull(.data$season_start_year)
stopifnot(
  nrow(ili_plus_data) == nrow(national_data),
  setequal(
    unique(ili_plus_data$season_start_year),
    primary_season_start_years()
  ),
  setequal(
    unique(clinical_only_ili_plus$season_start_year),
    expected_clinical_seasons
  ),
  !2015L %in% clinical_only_ili_plus$season_start_year,
  all(ili_plus_data$ili_plus_weighted >= 0 &
    ili_plus_data$ili_plus_weighted <= 1),
  all(ili_plus_data$ili_plus_unweighted >= 0 &
    ili_plus_data$ili_plus_unweighted <= 1)
)
ili_plus_draws <- draw_us_empirical_ili_plus_burden(
  ili_plus_data,
  n_draws = 3,
  ili_weighting = "cdc_population_weighted",
  lab_scope = "all_available",
  seed = 99,
  timing_shift_sd = 0.1
)
stopifnot(
  n_distinct(ili_plus_draws$draw) == 3L,
  all(is.finite(ili_plus_draws$burden)),
  all(ili_plus_draws$burden >= 0)
)

outcome_distribution <- flu_ve_onset_to_enrollment_distribution()
stopifnot(
  sum(outcome_distribution$n_observed) == 27135L,
  identical(outcome_distribution$lower_day, c(0L, 3L, 5L)),
  identical(outcome_distribution$upper_day, c(2L, 4L, 7L))
)

ve_draws <- draw_primary_ve_waning(
  ve_estimates,
  n_draws = 100,
  age_group = "all",
  reference_weeks = 8,
  seed = 101
)
ve_at_reference <- 1 - ve_draws$initial_effect_ratio_unbounded * exp(
  ve_draws$beta_wane_per_28d * ve_draws$ve_reference_weeks / 4
)
stopifnot(max(abs(ve_at_reference - ve_draws$reported_ve_draw)) < 1e-10)

for (waning_prior in c("spencer_ferdinands_fast", "spencer_slow")) {
  spencer_draws <- draw_primary_ve_waning(
    ve_estimates,
    n_draws = 100,
    age_group = "all",
    reference_weeks = 8,
    waning_prior = waning_prior,
    seed = 101
  )
  curve <- if (waning_prior == "spencer_ferdinands_fast") "fast" else "slow"
  ve_at_reference <- spencer_draws$initial_ve_unbounded *
    spencer_relative_ve_retention(spencer_draws$ve_reference_weeks, curve)
  stopifnot(max(abs(ve_at_reference - spencer_draws$reported_ve_draw)) < 1e-10)
}

small_state <- state_data %>%
  filter(
    .data$state %in% c("Alabama", "Alaska"),
    .data$season_start_year %in% c(2017L, 2018L)
  )
small_hhs_age <- hhs_age_data %>%
  filter(.data$season_start_year %in% c(2017L, 2018L))
small_population <- read_state_population(
  expected_states = unique(small_state$state)
) %>%
  filter(.data$state %in% unique(small_state$state))

model_data <- prepare_all_age_gam_data(small_state)

age_groups <- c("0-4", "5-24", "25-49", "50-64", "65+")
share_scenarios <- tibble::tibble(draw = 1:3, sampled_season = "2018/19")
age_shares <- bind_rows(lapply(age_groups, function(age_group) {
  draw_hhs_age_share(
    small_hhs_age,
    share_scenarios,
    target_age_group = age_group,
    seed = 100
  )
}))
share_sums <- age_shares %>%
  group_by(.data$draw, .data$season, .data$hhs_region, .data$week) %>%
  summarise(total = sum(.data$age_share), .groups = "drop")
stopifnot(max(abs(share_sums$total - 1)) < 1e-12)

state_models <- fit_statewise_latent_gams(
  model_data,
  family = "quasibinomial",
  global_k = 8,
  season_k = 5
)
all_age_burden <- draw_statewise_gam_decision_burden(
  state_models,
  state_population = small_population,
  n_draws = 3,
  seed = 102,
  timing_shift_sd = 0.5
)
no_shift_burden <- draw_statewise_gam_decision_burden(
  state_models,
  state_population = small_population,
  n_draws = 3,
  seed = 102,
  timing_shift_sd = 0
)

# Zero timing variance must preserve the paired seasons and coefficient draws.
# Coefficient sampling has its own seed, even though rnorm(sd = 0) uses no RNG.
paired_columns <- c(
  "draw", "state", "age_group", "season", "week", "max_week",
  "posterior_ili_prob", "model_source"
)
order_burden <- function(dat) {
  arrange(dat, .data$draw, .data$state, .data$age_group, .data$season, .data$week)
}
stopifnot(
  identical(
    select(order_burden(no_shift_burden), all_of(paired_columns)),
    select(order_burden(all_age_burden), all_of(paired_columns))
  ),
  all(no_shift_burden$timing_shift_weeks == 0),
  max(abs(no_shift_burden$burden - no_shift_burden$posterior_ili_prob)) < 1e-12,
  identical(
    order_burden(shift_burden_curves(no_shift_burden))$burden,
    order_burden(no_shift_burden)$burden
  )
)

# Applying the original shifts to the unshifted state curves must reproduce
# the paired shifted analysis, including its population-weighted US curve.
reapplied_shift_burden <- no_shift_burden %>%
  filter(.data$state != "US") %>%
  select(-"timing_shift_weeks") %>%
  left_join(
    distinct(all_age_burden, .data$draw, .data$timing_shift_weeks),
    by = "draw",
    relationship = "many-to-one"
  ) %>%
  shift_burden_curves() %>%
  add_population_weighted_us(small_population)
stopifnot(isTRUE(all.equal(
  select(order_burden(reapplied_shift_burden), all_of(names(all_age_burden))),
  order_burden(all_age_burden),
  tolerance = 1e-12
)))

age_burden <- draw_statewise_gam_age_decision_burden(
  state_models,
  ilinet_hhs_age = small_hhs_age,
  hhs_regions = hhs_region_crosswalk(),
  state_population = small_population,
  target_age_group = "65+",
  n_draws = 3,
  seed = 103,
  timing_shift_sd = 0.5
)
onset_aligned_burden <- apply_outcome_date_alignment(
  all_age_burden,
  prior = "flu_ve_outpatient_onset_to_enrollment",
  seed = 105
)

expected_states <- c("Alabama", "Alaska", "US")
stopifnot(
  setequal(unique(all_age_burden$state), expected_states),
  setequal(unique(age_burden$state), expected_states),
  max(all_age_burden$burden) < 1,
  max(age_burden$burden) < 1,
  all(onset_aligned_burden$outcome_date_delay_days %in% 0:7),
  all(onset_aligned_burden$outcome_date_offset_weeks <= 0),
  all(is.finite(onset_aligned_burden$burden))
)

test_ve <- ve_draws %>% filter(.data$draw <= 3)
test_burden <- all_age_burden %>% filter(.data$draw <= 3)
utilities <- evaluate_candidate_utilities(
  test_burden,
  test_ve,
  candidate_weeks = c(40L, 44L, 48L)
)
optimal <- posterior_optimal_weeks(utilities)
cutoffs <- prob_optimal_by_week_cutoff(optimal)
stopifnot(
  setequal(unique(cutoffs$cutoff_label), paste("MMWR week", c(39L, 44L, 48L, 52L))),
  all(cutoffs$probability >= 0 & cutoffs$probability <= 1)
)

spencer_ve <- draw_primary_ve_waning(
  ve_estimates,
  n_draws = 3,
  age_group = "all",
  reference_weeks = 8,
  waning_prior = "spencer_ferdinands_fast",
  seed = 104
)
spencer_utilities <- evaluate_candidate_utilities(
  test_burden,
  spencer_ve,
  candidate_weeks = c(40L, 44L, 48L)
)
stopifnot(all(is.finite(spencer_utilities$utility)))

message("Model component checks passed.")
