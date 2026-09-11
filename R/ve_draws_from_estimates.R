library(dplyr)
library(tidyr)
library(tibble)

source("R/protection_models.R")

draw_reported_ve_from_estimates <- function(
    ve_estimates,
    n_draws = 2000,
    age_group = "all",
    ve_age_group = age_group,
    seed = 20260507) {
  set.seed(seed)

  pool <- ve_estimates %>%
    filter(.data$ve_age_group == !!ve_age_group)

  if (nrow(pool) == 0) {
    stop("No VE estimates available for age group: ", ve_age_group)
  }

  log_effect_mean <- if ("log_or_mean" %in% names(pool)) "log_or_mean" else "log_rr_mean"
  log_effect_sd <- if ("log_or_sd" %in% names(pool)) "log_or_sd" else "log_rr_sd"
  missing <- setdiff(c(log_effect_mean, log_effect_sd), names(pool))
  if (length(missing) > 0L) {
    stop("VE estimates missing: ", paste(missing, collapse = ", "))
  }

  sampled <- pool[sample.int(nrow(pool), size = n_draws, replace = TRUE), ]
  reported_effect_ratio <- exp(stats::rnorm(
    n_draws,
    mean = sampled[[log_effect_mean]],
    sd = sampled[[log_effect_sd]]
  ))

  tibble(
    draw = seq_len(n_draws),
    age_group = age_group,
    reported_effect_ratio = reported_effect_ratio,
    reported_ve_draw = 1 - reported_effect_ratio,
    ve_source_season = sampled$season,
    ve_source = sampled$source_id
  )
}

calibrate_initial_ve <- function(ve_draws, reference_weeks) {
  reference_weeks <- as.numeric(reference_weeks)
  if (length(reference_weeks) != 1L || is.na(reference_weeks) || reference_weeks < 0) {
    stop("reference_weeks must be one nonnegative number.")
  }

  required <- c("reported_effect_ratio", "reported_ve_draw", "beta_wane_per_28d")
  missing <- setdiff(required, names(ve_draws))
  if (length(missing) > 0L) {
    stop("ve_draws missing: ", paste(missing, collapse = ", "))
  }

  if (!"protection_model" %in% names(ve_draws)) {
    ve_draws <- ve_draws %>%
      mutate(protection_model = "exponential_effect_ratio")
  }

  ve_draws %>%
    mutate(
      ve_reference_weeks = reference_weeks,
      reference_retention = case_when(
        .data$protection_model == "spencer_fast_relative_ve" ~
          spencer_relative_ve_retention(.data$ve_reference_weeks, "fast"),
        .data$protection_model == "spencer_slow_relative_ve" ~
          spencer_relative_ve_retention(.data$ve_reference_weeks, "slow"),
        TRUE ~ NA_real_
      ),
      initial_effect_ratio_unbounded = if_else(
        .data$protection_model == "exponential_effect_ratio",
        .data$reported_effect_ratio /
          exp(.data$beta_wane_per_28d * .data$ve_reference_weeks / 4),
        NA_real_
      ),
      initial_ve_unbounded = case_when(
        .data$protection_model == "exponential_effect_ratio" ~
          1 - .data$initial_effect_ratio_unbounded,
        TRUE ~ .data$reported_ve_draw / .data$reference_retention
      ),
      initial_ve = pmin(pmax(.data$initial_ve_unbounded, 0), 1),
      initial_effect_ratio = 1 - .data$initial_ve,
      initial_ve_bounded = .data$initial_ve_unbounded < 0 | .data$initial_ve_unbounded > 1
    )
}

draw_primary_waning <- function(
    n_draws = 2000,
    seed = 20260507,
    prior = Sys.getenv("WANING_PRIOR", "ray_2019")) {
  set.seed(seed)
  prior <- match.arg(
    prior,
    c(
      "ray_2019",
      "ray_2019_lower95",
      "ray_2019_upper95",
      "spencer_ferdinands_fast",
      "spencer_slow",
      "none"
    )
  )
  if (prior == "none") {
    return(tibble(
      draw = seq_len(n_draws),
      beta_wane_per_28d = 0,
      protection_model = protection_model_from_prior(prior),
      waning_source = "sensitivity_no_waning"
    ))
  }

  if (prior %in% c("spencer_ferdinands_fast", "spencer_slow")) {
    source_label <- if (prior == "spencer_ferdinands_fast") {
      "spencer_2026_ferdinands_fast_polynomial_relative_ve"
    } else {
      "spencer_2026_slow_polynomial_relative_ve"
    }
    return(tibble(
      draw = seq_len(n_draws),
      beta_wane_per_28d = NA_real_,
      protection_model = protection_model_from_prior(prior),
      waning_source = source_label
    ))
  }

  se <- (log(1.20) - log(1.13)) / (2 * qnorm(0.975))
  if (prior == "ray_2019_lower95") {
    return(tibble(
      draw = seq_len(n_draws),
      beta_wane_per_28d = log(1.13),
      protection_model = protection_model_from_prior(prior),
      waning_source = "ray_2019_lower95_or_per_28d_fixed"
    ))
  }
  if (prior == "ray_2019_upper95") {
    return(tibble(
      draw = seq_len(n_draws),
      beta_wane_per_28d = log(1.20),
      protection_model = protection_model_from_prior(prior),
      waning_source = "ray_2019_upper95_or_per_28d_fixed"
    ))
  }

  tibble(
    draw = seq_len(n_draws),
    beta_wane_per_28d = rnorm(n_draws, mean = log(1.16), sd = se),
    protection_model = protection_model_from_prior(prior),
    waning_source = "ray_2019_continuous_or_per_28d_ci_1.13_1.20"
  )
}

draw_immune_response_lag <- function(
    n_draws = 2000,
    seed = 20260507,
    prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w")) {
  set.seed(seed)
  prior <- match.arg(prior, c("discrete_about_2w", "fixed_2w"))

  if (prior == "fixed_2w") {
    return(tibble(
      draw = seq_len(n_draws),
      immune_lag_weeks = 2,
      immune_lag_source = "cdc_acip_about_2_weeks_fixed"
    ))
  }

  tibble(
    draw = seq_len(n_draws),
    immune_lag_weeks = sample(c(1, 2, 3), size = n_draws, replace = TRUE, prob = c(0.15, 0.70, 0.15)),
    immune_lag_source = "cdc_acip_about_2_weeks_discrete_prior"
  )
}

draw_primary_ve_waning <- function(
    ve_estimates,
    n_draws = 2000,
    age_group = "all",
    ve_age_group = age_group,
    reference_weeks = as.numeric(Sys.getenv("VE_REFERENCE_WEEKS", "8")),
    waning_prior = Sys.getenv("WANING_PRIOR", "ray_2019"),
    immune_lag_prior = Sys.getenv("IMMUNE_LAG_PRIOR", "discrete_about_2w"),
    seed = 20260507) {
  reported <- draw_reported_ve_from_estimates(
    ve_estimates = ve_estimates,
    n_draws = n_draws,
    age_group = age_group,
    ve_age_group = ve_age_group,
    seed = seed
  )
  waning <- draw_primary_waning(
    n_draws = n_draws,
    seed = seed + 1L,
    prior = waning_prior
  )
  immune_lag <- draw_immune_response_lag(
    n_draws = n_draws,
    seed = seed + 2L,
    prior = immune_lag_prior
  )

  reported %>%
    left_join(waning, by = "draw") %>%
    left_join(immune_lag, by = "draw") %>%
    calibrate_initial_ve(reference_weeks = reference_weeks)
}
