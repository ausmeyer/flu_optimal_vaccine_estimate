library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("R/calendar.R")
source("R/ve_priors.R")
source("R/protection_models.R")

evaluate_candidate_utilities <- function(
    burden_draws,
    ve_draws,
    candidate_weeks = c(36:52, 1:12),
    immune_lag_weeks = 2L) {
  required_burden <- c("draw", "state", "age_group", "week", "burden")
  required_ve <- c("draw", "age_group", "initial_ve", "beta_wane_per_28d")

  missing_burden <- setdiff(required_burden, names(burden_draws))
  missing_ve <- setdiff(required_ve, names(ve_draws))
  if (length(missing_burden) > 0) {
    stop("burden_draws missing: ", paste(missing_burden, collapse = ", "))
  }
  if (length(missing_ve) > 0) {
    stop("ve_draws missing: ", paste(missing_ve, collapse = ", "))
  }

  if ("max_week" %in% names(burden_draws)) {
    burden_eval <- burden_draws
  } else if ("season" %in% names(burden_draws)) {
    season_max_week <- burden_draws %>%
      group_by(.data$season) %>%
      summarise(max_week = if_else(any(.data$week == 53L), 53L, 52L), .groups = "drop")

    burden_eval <- burden_draws %>%
      left_join(season_max_week, by = "season")
  } else {
    burden_eval <- burden_draws %>%
      mutate(max_week = 52L)
  }

  burden_eval <- burden_eval %>%
    filter(is_in_burden_window(.data$week)) %>%
    select("draw", "state", "age_group", "week", "burden", "max_week") %>%
    group_by(.data$draw, .data$state, .data$age_group, .data$week, .data$max_week) %>%
    summarise(burden = sum(.data$burden, na.rm = TRUE), .groups = "drop")

  ve_eval <- ve_draws
  if (!"immune_lag_weeks" %in% names(ve_eval)) {
    ve_eval <- ve_eval %>% mutate(immune_lag_weeks = as.numeric(immune_lag_weeks))
  } else {
    ve_eval <- ve_eval %>%
      mutate(immune_lag_weeks = coalesce(.data$immune_lag_weeks, as.numeric(immune_lag_weeks)))
  }
  if (!"protection_model" %in% names(ve_eval)) {
    ve_eval <- ve_eval %>%
      mutate(protection_model = "exponential_effect_ratio")
  }
  ve_eval <- ve_eval %>%
    select(
      "draw", "age_group", "initial_ve", "beta_wane_per_28d",
      "immune_lag_weeks", "protection_model"
    )

  candidate_weeks <- as.integer(candidate_weeks)

  eval_group <- function(dat) {
    state_value <- dat$state[[1]]
    age_value <- dat$age_group[[1]]
    max_week_value <- dat$max_week[[1]]
    draw_ids <- sort(unique(dat$draw))
    burden_weeks <- sort(unique(dat$week))

    ve_sub <- ve_eval %>%
      filter(.data$age_group == !!age_value, .data$draw %in% draw_ids) %>%
      arrange(.data$draw)
    ve_match <- match(draw_ids, ve_sub$draw)
    if (anyNA(ve_match)) {
      stop("Missing VE draw(s) for age group ", age_value)
    }
    initial_ve <- ve_sub$initial_ve[ve_match]
    beta_wane_per_28d <- ve_sub$beta_wane_per_28d[ve_match]
    immune_lag_by_draw <- ve_sub$immune_lag_weeks[ve_match]
    protection_model <- unique(ve_sub$protection_model[ve_match])
    if (length(protection_model) != 1L) {
      stop("Each analysis must use one protection model.")
    }

    burden_matrix <- matrix(0, nrow = length(draw_ids), ncol = length(burden_weeks))
    burden_matrix[cbind(match(dat$draw, draw_ids), match(dat$week, burden_weeks))] <- dat$burden

    utility_matrix <- matrix(0, nrow = length(draw_ids), ncol = length(candidate_weeks))

    for (j in seq_along(candidate_weeks)) {
      weeks_since <- weeks_since_full_response(
        vaccination_week = candidate_weeks[[j]],
        burden_week = burden_weeks,
        immune_lag_weeks = 0L,
        max_week = max_week_value
      )
      weeks_since <- sweep(
        matrix(rep(weeks_since, each = length(draw_ids)), nrow = length(draw_ids)),
        1,
        immune_lag_by_draw,
        `-`
      )
      months_since <- pmax(weeks_since, 0) / 4
      raw_protection <- switch(
        protection_model,
        exponential_effect_ratio =
          1 - (1 - initial_ve) * exp(beta_wane_per_28d * months_since),
        spencer_fast_relative_ve =
          initial_ve * spencer_relative_ve_retention(weeks_since, "fast"),
        spencer_slow_relative_ve =
          initial_ve * spencer_relative_ve_retention(weeks_since, "slow"),
        stop("Unknown protection model: ", protection_model)
      )
      protection <- pmin(pmax(raw_protection, 0), 1)
      if (any(weeks_since < 0)) {
        protection[weeks_since < 0] <- 0
      }
      utility_matrix[, j] <- rowSums(burden_matrix * protection, na.rm = TRUE)
    }

    tibble(
      draw = rep(draw_ids, each = length(candidate_weeks)),
      state = state_value,
      age_group = age_value,
      vaccination_week = rep(candidate_weeks, times = length(draw_ids)),
      utility = as.vector(t(utility_matrix))
    )
  }

  burden_eval %>%
    group_by(.data$state, .data$age_group, .data$max_week) %>%
    group_split() %>%
    map_dfr(eval_group) %>%
    group_by(.data$draw, .data$state, .data$age_group, .data$vaccination_week) %>%
    summarise(utility = sum(.data$utility, na.rm = TRUE), .groups = "drop")
}

posterior_optimal_weeks <- function(
    candidate_utilities,
    tie_tolerance = 1e-9,
    seed = 20260507) {
  set.seed(seed)

  candidate_utilities %>%
    group_by(.data$draw, .data$state, .data$age_group) %>%
    mutate(
      max_utility = max(.data$utility, na.rm = TRUE),
      utility_tolerance = tie_tolerance * pmax(1, abs(.data$max_utility)),
      is_tied_max = abs(.data$utility - .data$max_utility) <= .data$utility_tolerance,
      n_tied_max = sum(.data$is_tied_max),
      tie_break = stats::runif(n())
    ) %>%
    filter(.data$is_tied_max) %>%
    slice_max(.data$tie_break, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      draw = .data$draw,
      state = .data$state,
      age_group = .data$age_group,
      optimal_week = .data$vaccination_week,
      optimal_season_index = week_to_season_index(.data$vaccination_week),
      optimal_utility = .data$utility,
      n_tied_max = .data$n_tied_max
    )
}

summarise_optimal_week_distribution <- function(optimal_weeks) {
  n_draws <- optimal_weeks %>%
    distinct(.data$draw, .data$state, .data$age_group) %>%
    count(.data$state, .data$age_group, name = "n_draws")

  optimal_weeks %>%
    count(.data$state, .data$age_group, .data$optimal_week, name = "n") %>%
    complete(
      state,
      age_group,
      optimal_week = 1:52,
      fill = list(n = 0L)
    ) %>%
    left_join(n_draws, by = c("state", "age_group")) %>%
    mutate(probability = .data$n / .data$n_draws)
}

summarise_optimal_week_intervals <- function(optimal_weeks) {
  optimal_weeks %>%
    group_by(.data$state, .data$age_group) %>%
    summarise(
      median_index = round(stats::median(.data$optimal_season_index)),
      lower_50_index = round(stats::quantile(.data$optimal_season_index, 0.25)),
      upper_50_index = round(stats::quantile(.data$optimal_season_index, 0.75)),
      lower_95_index = round(stats::quantile(.data$optimal_season_index, 0.025)),
      upper_95_index = round(stats::quantile(.data$optimal_season_index, 0.975)),
      .groups = "drop"
    ) %>%
    mutate(
      median_week = season_index_to_week(.data$median_index),
      lower_50_week = season_index_to_week(.data$lower_50_index),
      upper_50_week = season_index_to_week(.data$upper_50_index),
      lower_95_week = season_index_to_week(.data$lower_95_index),
      upper_95_week = season_index_to_week(.data$upper_95_index)
    )
}

prob_optimal_by_week_cutoff <- function(
    optimal_weeks,
    cutoff_weeks = c(39L, 44L, 48L, 52L)) {
  cutoffs <- tibble(cutoff_week = as.integer(cutoff_weeks)) %>%
    mutate(
      cutoff_label = paste("MMWR week", .data$cutoff_week),
      cutoff_index = week_to_season_index(.data$cutoff_week)
    ) %>%
    select("cutoff_week", "cutoff_label", "cutoff_index")

  optimal_weeks %>%
    tidyr::crossing(cutoffs) %>%
    group_by(.data$state, .data$age_group, .data$cutoff_week, .data$cutoff_label) %>%
    summarise(
      probability = mean(.data$optimal_season_index <= .data$cutoff_index),
      .groups = "drop"
    )
}

summarise_regret_curve <- function(candidate_utilities) {
  best <- candidate_utilities %>%
    group_by(.data$draw, .data$state, .data$age_group) %>%
    summarise(best_utility = max(.data$utility, na.rm = TRUE), .groups = "drop")

  candidate_utilities %>%
    left_join(best, by = c("draw", "state", "age_group")) %>%
    mutate(regret = .data$best_utility - .data$utility) %>%
    group_by(.data$state, .data$age_group, .data$vaccination_week) %>%
    summarise(
      mean_regret = mean(.data$regret, na.rm = TRUE),
      median_regret = stats::median(.data$regret, na.rm = TRUE),
      lower_50 = stats::quantile(.data$regret, 0.25, na.rm = TRUE),
      upper_50 = stats::quantile(.data$regret, 0.75, na.rm = TRUE),
      lower_95 = stats::quantile(.data$regret, 0.025, na.rm = TRUE),
      upper_95 = stats::quantile(.data$regret, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_near_optimal_curve <- function(
    candidate_utilities,
    relative_regret_tolerances = c(0.01, 0.05, 0.10)) {
  best <- candidate_utilities %>%
    group_by(.data$draw, .data$state, .data$age_group) %>%
    summarise(best_utility = max(.data$utility, na.rm = TRUE), .groups = "drop")

  candidate_utilities %>%
    left_join(best, by = c("draw", "state", "age_group")) %>%
    mutate(
      regret = .data$best_utility - .data$utility,
      relative_regret = if_else(
        .data$best_utility > 0,
        .data$regret / .data$best_utility,
        0
      )
    ) %>%
    tidyr::crossing(relative_regret_tolerance = relative_regret_tolerances) %>%
    group_by(.data$state, .data$age_group, .data$vaccination_week, .data$relative_regret_tolerance) %>%
    summarise(
      probability_near_optimal = mean(.data$relative_regret <= .data$relative_regret_tolerance, na.rm = TRUE),
      mean_relative_regret = mean(.data$relative_regret, na.rm = TRUE),
      median_relative_regret = stats::median(.data$relative_regret, na.rm = TRUE),
      lower_50_relative_regret = stats::quantile(.data$relative_regret, 0.25, na.rm = TRUE),
      upper_50_relative_regret = stats::quantile(.data$relative_regret, 0.75, na.rm = TRUE),
      lower_95_relative_regret = stats::quantile(.data$relative_regret, 0.025, na.rm = TRUE),
      upper_95_relative_regret = stats::quantile(.data$relative_regret, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

aggregate_us_population <- function(burden_draws, state_weights = NULL) {
  if (is.null(state_weights)) {
    stop("state_weights are required for national aggregation.")
  }
  if (!all(c("state", "weight") %in% names(state_weights))) {
    stop("state_weights must include state and weight.")
  }
  join_columns <- intersect(c("state", "age_group"), names(state_weights))
  group_columns <- intersect(
    c("draw", "age_group", "season", "week", "max_week"),
    names(burden_draws)
  )

  burden_draws %>%
    left_join(state_weights, by = join_columns, relationship = "many-to-one") %>%
    { if (anyNA(.$weight)) stop("Missing state weight after aggregation join.") else . } %>%
    group_by(across(all_of(group_columns))) %>%
    summarise(
      burden = stats::weighted.mean(.data$burden, .data$weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(state = "US")
}
