library(dplyr)
library(tidyr)
library(tibble)

source("R/calendar.R")
source("R/latent_curve_gam.R")

complete_common_seasons <- function(data, states = NULL, complete_only = TRUE) {
  if (is.null(states)) {
    available <- data %>%
      group_by(.data$season) %>%
      summarise(complete = any(.data$mmwr_week == 22L), .groups = "drop")
  } else {
    by_state <- data %>%
      filter(.data$state %in% states) %>%
      group_by(.data$season, .data$state) %>%
      summarise(complete = any(.data$mmwr_week == 22L), .groups = "drop")
    available <- by_state %>%
      group_by(.data$season) %>%
      summarise(
        n_states = n_distinct(.data$state),
        complete = n_states == length(states) && all(.data$complete),
        .groups = "drop"
      )
  }

  seasons <- available %>%
    filter(!complete_only | .data$complete) %>%
    pull(.data$season)
  if (length(seasons) == 0L && complete_only) {
    return(complete_common_seasons(data, states = states, complete_only = FALSE))
  }
  if (length(seasons) == 0L) {
    stop("No influenza seasons are available for empirical sampling.")
  }
  sort(as.character(seasons))
}

sample_empirical_scenarios <- function(
    data,
    n_draws,
    seed,
    states = NULL,
    complete_only = TRUE,
    timing_shift_sd = 0.75) {
  set.seed(seed)
  seasons <- complete_common_seasons(data, states = states, complete_only = complete_only)
  tibble(
    draw = seq_len(n_draws),
    sampled_season = sample(seasons, n_draws, replace = TRUE),
    timing_shift_weeks = stats::rnorm(n_draws, mean = 0, sd = timing_shift_sd)
  )
}

transform_empirical_burden_draws <- function(
    burden_draws,
    timing_shift_sd = 0.75,
    curve_smoothing = c("cubic_spline", "none"),
    spline_spar = 0.65) {
  curve_smoothing <- match.arg(curve_smoothing)

  burden_draws %>%
    group_by(.data$draw, .data$state, .data$age_group, .data$season) %>%
    group_modify(function(dat, key) {
      max_week <- if (any(dat$week == 53L)) 53L else 52L
      season_pos <- week_to_season_index(
        dat$week,
        start_week = 36L,
        end_week = 22L,
        max_week = max_week
      )
      ordered <- order(season_pos)
      shift <- if ("timing_shift_weeks" %in% names(dat)) {
        dat$timing_shift_weeks[[1]]
      } else if (timing_shift_sd > 0) {
        stats::rnorm(1L, mean = 0, sd = timing_shift_sd)
      } else {
        0
      }

      if (curve_smoothing == "none" || length(unique(season_pos)) < 4L) {
        # Preserve the sampled weekly values. Linear interpolation is used
        # only to apply the continuous timing shift; rule = 2 holds the
        # nearest sampled boundary value outside the observed week range.
        shifted <- stats::approx(
          x = season_pos[ordered],
          y = dat$burden[ordered],
          xout = season_pos - shift,
          rule = 2,
          ties = mean
        )$y
      } else {
        fit <- stats::smooth.spline(
          x = season_pos[ordered],
          y = log1p(pmax(dat$burden[ordered], 0)),
          spar = spline_spar
        )
        shifted <- expm1(stats::predict(fit, x = season_pos - shift)$y)
      }

      dat$burden <- pmax(shifted, 0)
      dat$timing_shift_weeks <- shift
      dat$max_week <- max_week
      dat
    }) %>%
    ungroup()
}

draw_empirical_ilinet_burden <- function(
    ilinet_state,
    n_draws = 2000,
    target_states = NULL,
    age_group = "all",
    seed = 20260507,
    include_observation_noise = TRUE,
    timing_shift_sd = 0.75,
    curve_smoothing = c("cubic_spline", "none"),
    spline_spar = 0.65,
    apply_timing_shift = TRUE,
    sample_complete_seasons_only = tolower(Sys.getenv(
      "SAMPLE_COMPLETE_SEASONS_ONLY", "true"
    )) %in% c("1", "true", "yes")) {
  curve_smoothing <- match.arg(curve_smoothing)
  if (is.null(target_states)) {
    target_states <- sort(unique(ilinet_state$state))
  }
  observed_ilinet <- ilinet_state %>%
    filter(
      is.finite(.data$ilitotal),
      is.finite(.data$total_patients),
      .data$total_patients > 0,
      .data$ilitotal >= 0,
      .data$ilitotal <= .data$total_patients
    )
  n_dropped <- nrow(ilinet_state) - nrow(observed_ilinet)
  if (n_dropped > 0L) {
    message(
      "Treating ", n_dropped,
      " state-week row(s) with no valid outpatient denominator as missing."
    )
  }
  missing_states <- setdiff(target_states, unique(observed_ilinet$state))
  if (length(missing_states) > 0L) {
    stop("No valid empirical observations for: ", paste(missing_states, collapse = ", "))
  }
  scenarios <- sample_empirical_scenarios(
    data = observed_ilinet,
    n_draws = n_draws,
    seed = seed,
    states = target_states,
    complete_only = sample_complete_seasons_only,
    timing_shift_sd = timing_shift_sd
  )

  set.seed(seed + 1L)
  draws <- scenarios %>%
    inner_join(
      observed_ilinet %>% filter(.data$state %in% target_states),
      by = c("sampled_season" = "season"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      posterior_ili_prob = if (include_observation_noise) {
        stats::rbeta(
          n(),
          shape1 = pmax(.data$ilitotal, 0) + 0.5,
          shape2 = pmax(.data$total_patients - .data$ilitotal, 0) + 0.5
        )
      } else {
        .data$ilitotal / .data$total_patients
      },
      burden = .data$posterior_ili_prob,
      age_group = age_group,
      season = .data$sampled_season,
      week = .data$mmwr_week,
      model_source = "empirical_state_ilinet_proportion"
    ) %>%
    select(
      "draw", "state", "age_group", "season", "week", "burden",
      "posterior_ili_prob", "model_source", "timing_shift_weeks"
    )

  if (apply_timing_shift) {
    transform_empirical_burden_draws(
      draws,
      timing_shift_sd = timing_shift_sd,
      curve_smoothing = curve_smoothing,
      spline_spar = spline_spar
    )
  } else {
    draws
  }
}

draw_us_empirical_ilinet_burden <- function(
    ilinet_national,
    n_draws = 2000,
    age_group = "all",
    seed = 20260507,
    include_observation_noise = TRUE,
    timing_shift_sd = 0.75,
    curve_smoothing = c("cubic_spline", "none"),
    spline_spar = 0.65,
    sample_complete_seasons_only = tolower(Sys.getenv(
      "SAMPLE_COMPLETE_SEASONS_ONLY", "true"
    )) %in% c("1", "true", "yes")) {
  curve_smoothing <- match.arg(curve_smoothing)
  if (!"weighted_ili" %in% names(ilinet_national)) {
    stop("National ILINet data must include CDC weighted_ili.")
  }
  invalid_counts <- !is.finite(ilinet_national$ilitotal) |
    !is.finite(ilinet_national$total_patients) |
    ilinet_national$total_patients <= 0 |
    ilinet_national$ilitotal < 0 |
    ilinet_national$ilitotal > ilinet_national$total_patients
  if (any(invalid_counts)) {
    stop("National empirical ILINet data contain invalid outpatient counts.")
  }
  scenarios <- sample_empirical_scenarios(
    data = ilinet_national,
    n_draws = n_draws,
    seed = seed,
    complete_only = sample_complete_seasons_only,
    timing_shift_sd = timing_shift_sd
  )

  set.seed(seed + 1L)
  scenarios %>%
    inner_join(
      ilinet_national,
      by = c("sampled_season" = "season"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      sampled_unweighted_prob = if (include_observation_noise) {
        stats::rbeta(
          n(),
          shape1 = pmax(.data$ilitotal, 0) + 0.5,
          shape2 = pmax(.data$total_patients - .data$ilitotal, 0) + 0.5
        )
      } else {
        .data$ilitotal / .data$total_patients
      },
      weighting_ratio = (.data$weighted_ili / 100) /
        pmax(.data$unweighted_ili / 100, 1e-8),
      posterior_ili_prob = pmin(pmax(.data$sampled_unweighted_prob * .data$weighting_ratio, 0), 1),
      burden = .data$posterior_ili_prob,
      state = "US",
      age_group = age_group,
      season = .data$sampled_season,
      week = .data$mmwr_week,
      model_source = "empirical_cdc_population_weighted_national_ilinet"
    ) %>%
    select(
      "draw", "state", "age_group", "season", "week", "burden",
      "posterior_ili_prob", "model_source", "timing_shift_weeks"
    ) %>%
    transform_empirical_burden_draws(
      timing_shift_sd = timing_shift_sd,
      curve_smoothing = curve_smoothing,
      spline_spar = spline_spar
    )
}

draw_empirical_age_burden <- function(
    ilinet_state,
    ilinet_hhs_age,
    hhs_regions,
    state_population,
    n_draws = 2000,
    target_age_group,
    seed = 20260507,
    include_observation_noise = TRUE,
    timing_shift_sd = 0.75,
    curve_smoothing = c("cubic_spline", "none"),
    spline_spar = 0.65,
    dirichlet_offset = 0.5,
    sample_complete_seasons_only = tolower(Sys.getenv(
      "SAMPLE_COMPLETE_SEASONS_ONLY", "true"
    )) %in% c("1", "true", "yes")) {
  curve_smoothing <- match.arg(curve_smoothing)
  all_age_draws <- draw_empirical_ilinet_burden(
    ilinet_state = ilinet_state,
    n_draws = n_draws,
    seed = seed,
    include_observation_noise = include_observation_noise,
    timing_shift_sd = timing_shift_sd,
    curve_smoothing = curve_smoothing,
    spline_spar = spline_spar,
    apply_timing_shift = FALSE,
    sample_complete_seasons_only = sample_complete_seasons_only
  )
  scenarios <- all_age_draws %>%
    distinct(.data$draw, sampled_season = .data$season, .data$timing_shift_weeks)
  age_share <- draw_hhs_age_share(
    ilinet_hhs_age = ilinet_hhs_age,
    scenarios = scenarios,
    target_age_group = target_age_group,
    seed = seed + 2L,
    dirichlet_offset = dirichlet_offset
  )
  mapping <- state_hhs_mapping(hhs_regions, unique(all_age_draws$state))

  age_draws <- all_age_draws %>%
    left_join(mapping, by = "state", relationship = "many-to-one") %>%
    inner_join(
      age_share,
      by = c("draw", "season", "hhs_region", "week"),
      relationship = "many-to-one"
    ) %>%
    transmute(
      draw = .data$draw,
      state = .data$state,
      age_group = target_age_group,
      season = .data$season,
      week = .data$week,
      burden = .data$burden * .data$age_share,
      posterior_ili_prob = NA_real_,
      model_source = "empirical_state_ilinet_with_hhs_age_composition",
      timing_shift_weeks = .data$timing_shift_weeks
    ) %>%
    transform_empirical_burden_draws(
      timing_shift_sd = timing_shift_sd,
      curve_smoothing = curve_smoothing,
      spline_spar = spline_spar
    )

  add_population_weighted_us(age_draws, state_population)
}
