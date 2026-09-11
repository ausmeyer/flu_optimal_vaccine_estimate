library(dplyr)
library(tidyr)
library(tibble)
library(mgcv)

source("R/calendar.R")

prepare_all_age_gam_data <- function(ilinet_state) {
  season_max_week <- ilinet_state %>%
    group_by(.data$season) %>%
    summarise(max_week = if_else(any(.data$mmwr_week == 53L), 53L, 52L), .groups = "drop")

  ilinet_state %>%
    left_join(season_max_week, by = "season") %>%
    filter(is_in_burden_window(.data$mmwr_week)) %>%
    mutate(
      state = factor(.data$state),
      season = factor(.data$season),
      season_pos = week_to_season_index(
        .data$mmwr_week,
        start_week = 36L,
        end_week = 22L,
        max_week = .data$max_week
      ),
      ili_cases = pmax(round(.data$ilitotal), 0),
      non_ili_cases = pmax(round(.data$total_patients - .data$ilitotal), 0),
      has_observed_denominator = is.finite(.data$total_patients) &
        .data$total_patients > 0 &
        is.finite(.data$ilitotal) &
        .data$ilitotal >= 0 &
        .data$ilitotal <= .data$total_patients,
      total_patients_model = if_else(
        .data$has_observed_denominator,
        pmax(round(.data$total_patients), 1),
        0
      ),
      ili_prop = if_else(
        .data$has_observed_denominator,
        pmin(pmax(.data$ili_cases / .data$total_patients_model, 0), 1),
        NA_real_
      ),
      model_ili_prop = coalesce(.data$ili_prop, 0),
      age_group = "all"
    )
}

fit_statewise_latent_gams <- function(
    model_data,
    family = c("quasibinomial", "binomial"),
    global_k = as.integer(Sys.getenv("GAM_GLOBAL_K", "14")),
    season_k = as.integer(Sys.getenv("GAM_SEASON_K", "12")),
    gamma = as.numeric(Sys.getenv("GAM_GAMMA", "1"))) {
  family <- match.arg(family)
  split_data <- split(
    model_data %>% mutate(state = droplevels(factor(.data$state))),
    model_data$state,
    drop = TRUE
  )

  lapply(split_data, function(dat) {
    dat <- dat %>%
      mutate(
        season = droplevels(factor(.data$season)),
        state = droplevels(factor(.data$state))
      )

    response <- if (family == "quasibinomial") {
      "model_ili_prop"
    } else {
      "cbind(ili_cases, non_ili_cases)"
    }
    formula <- stats::as.formula(sprintf(
      "%s ~ season + s(season_pos, k = %d, bs = 'cr') +
       s(season_pos, season, bs = 'fs', k = %d)",
      response,
      global_k,
      season_k
    ))

    args <- list(
      formula = formula,
      family = if (family == "quasibinomial") stats::quasibinomial() else stats::binomial(),
      data = dat,
      method = "REML",
      gamma = gamma
    )
    if (family == "quasibinomial") {
      args$weights <- dat$total_patients_model
    }

    list(
      fit = do.call(mgcv::gam, args),
      data = dat,
      smoothing_config = list(
        family = family,
        global_k = global_k,
        season_k = season_k,
        gamma = gamma
      )
    )
  })
}

make_prediction_grid <- function(model_data) {
  model_data %>%
    distinct(.data$state, .data$season, .data$mmwr_week, .data$season_pos, .data$age_group) %>%
    arrange(.data$state, .data$season, .data$season_pos) %>%
    mutate(
      state = factor(.data$state, levels = levels(model_data$state)),
      season = factor(.data$season, levels = levels(model_data$season))
    )
}

draw_latent_gam_burden <- function(
    fit,
    model_data,
    n_draws = 1000,
    seed = 20260507,
    sampled_seasons = NULL,
    source = "all_age_ilinet_gam") {
  set.seed(seed)
  pred_grid <- make_prediction_grid(model_data)
  design_matrix <- stats::predict(fit, newdata = pred_grid, type = "lpmatrix")
  coefficient_draws <- MASS::mvrnorm(
    n = n_draws,
    mu = stats::coef(fit),
    Sigma = fit$Vp
  )

  if (!is.null(sampled_seasons) && length(sampled_seasons) != n_draws) {
    stop("sampled_seasons must have one value per draw.")
  }

  bind_rows(lapply(seq_len(n_draws), function(draw_id) {
    rows <- if (is.null(sampled_seasons)) {
      seq_len(nrow(pred_grid))
    } else {
      which(as.character(pred_grid$season) == sampled_seasons[[draw_id]])
    }
    if (length(rows) == 0L) {
      stop("Sampled season is absent from the model prediction grid.")
    }
    posterior_ili_prob <- stats::plogis(
      as.numeric(design_matrix[rows, , drop = FALSE] %*% coefficient_draws[draw_id, ])
    )
    pred_grid[rows, , drop = FALSE] %>%
      transmute(
        draw = draw_id,
        state = as.character(.data$state),
        age_group = as.character(.data$age_group),
        season = as.character(.data$season),
        week = .data$mmwr_week,
        burden = posterior_ili_prob,
        posterior_ili_prob = posterior_ili_prob,
        model_source = source
      )
  }))
}

available_decision_seasons <- function(state_models, complete_only = TRUE) {
  by_state <- lapply(state_models, function(obj) {
    obj$data %>%
      mutate(season = as.character(.data$season)) %>%
      group_by(.data$season) %>%
      summarise(complete = any(.data$mmwr_week == 22L), .groups = "drop") %>%
      filter(!complete_only | .data$complete) %>%
      pull(.data$season)
  })
  seasons <- Reduce(intersect, by_state)
  if (length(seasons) == 0L && complete_only) {
    return(available_decision_seasons(state_models, complete_only = FALSE))
  }
  if (length(seasons) == 0L) {
    stop("The state models have no common influenza seasons.")
  }
  sort(seasons)
}

sample_decision_scenarios <- function(
    state_models,
    n_draws,
    seed,
    sample_complete_seasons_only = TRUE,
    timing_shift_sd = 0.75) {
  set.seed(seed)
  seasons <- available_decision_seasons(
    state_models,
    complete_only = sample_complete_seasons_only
  )
  tibble(
    draw = seq_len(n_draws),
    sampled_season = sample(seasons, n_draws, replace = TRUE),
    timing_shift_weeks = stats::rnorm(n_draws, mean = 0, sd = timing_shift_sd)
  )
}

shift_burden_curves <- function(burden_draws) {
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
      shift <- dat$timing_shift_weeks[[1]]
      ordered <- order(season_pos)
      dat$burden <- pmax(stats::approx(
        x = season_pos[ordered],
        y = dat$burden[ordered],
        xout = season_pos - shift,
        rule = 2,
        ties = mean
      )$y, 0)
      dat$max_week <- max_week
      dat
    }) %>%
    ungroup()
}

add_population_weighted_us <- function(state_draws, state_population) {
  required <- c("state", "population")
  missing <- setdiff(required, names(state_population))
  if (length(missing) > 0L) {
    stop("state_population missing: ", paste(missing, collapse = ", "))
  }

  weighted <- state_draws %>%
    left_join(
      state_population %>% select("state", "population"),
      by = "state",
      relationship = "many-to-one"
    )
  missing_states <- weighted %>%
    filter(is.na(.data$population)) %>%
    distinct(.data$state) %>%
    pull(.data$state)
  if (length(missing_states) > 0L) {
    stop("Population weights missing for: ", paste(missing_states, collapse = ", "))
  }

  us_draws <- weighted %>%
    group_by(
      .data$draw,
      .data$age_group,
      .data$season,
      .data$week,
      .data$max_week,
      .data$timing_shift_weeks
    ) %>%
    summarise(
      burden = stats::weighted.mean(.data$burden, .data$population, na.rm = TRUE),
      posterior_ili_prob = if (all(is.na(.data$posterior_ili_prob))) {
        NA_real_
      } else {
        stats::weighted.mean(.data$posterior_ili_prob, .data$population, na.rm = TRUE)
      },
      model_source = first(.data$model_source),
      .groups = "drop"
    ) %>%
    mutate(state = "US", .before = "age_group")

  bind_rows(
    weighted %>% select(-"population"),
    us_draws
  )
}

draw_statewise_gam_decision_burden <- function(
    state_models,
    state_population = NULL,
    n_draws = 1000,
    seed = 20260507,
    sample_complete_seasons_only = tolower(Sys.getenv(
      "SAMPLE_COMPLETE_SEASONS_ONLY", "true"
    )) %in% c("1", "true", "yes"),
    timing_shift_sd = as.numeric(Sys.getenv("TIMING_SHIFT_SD", "0.75"))) {
  scenarios <- sample_decision_scenarios(
    state_models = state_models,
    n_draws = n_draws,
    seed = seed,
    sample_complete_seasons_only = sample_complete_seasons_only,
    timing_shift_sd = timing_shift_sd
  )

  set.seed(seed + 1L)
  state_draws <- bind_rows(lapply(names(state_models), function(state_name) {
    obj <- state_models[[state_name]]
    draw_latent_gam_burden(
      fit = obj$fit,
      model_data = obj$data,
      n_draws = n_draws,
      seed = sample.int(.Machine$integer.max, 1L),
      sampled_seasons = scenarios$sampled_season
    )
  })) %>%
    left_join(
      scenarios %>% select("draw", "timing_shift_weeks"),
      by = "draw",
      relationship = "many-to-one"
    ) %>%
    shift_burden_curves()

  if (is.null(state_population)) {
    state_draws
  } else {
    add_population_weighted_us(state_draws, state_population)
  }
}

state_hhs_mapping <- function(hhs_regions, states) {
  hhs_regions %>%
    transmute(
      state = as.character(.data$state_or_territory),
      hhs_region = as.character(.data$region)
    ) %>%
    filter(.data$state %in% states) %>%
    distinct(.data$state, .data$hhs_region)
}

draw_hhs_age_share <- function(
    ilinet_hhs_age,
    scenarios,
    target_age_group,
    seed,
    dirichlet_offset = 0.5) {
  set.seed(seed)
  age_data <- ilinet_hhs_age %>%
    mutate(
      season = as.character(.data$season),
      hhs_region = as.character(.data$hhs_region)
    ) %>%
    group_by(.data$season, .data$hhs_region, .data$mmwr_week, .data$age_group) %>%
    summarise(ili_age_count = sum(.data$ili_age_count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = "age_group",
      values_from = "ili_age_count",
      values_fill = 0
    )

  key_columns <- c("season", "hhs_region", "mmwr_week")
  age_columns <- setdiff(names(age_data), key_columns)
  missing_age <- setdiff(target_age_group, age_columns)
  if (length(missing_age) > 0L) {
    stop("HHS ILINet data do not include age group: ", target_age_group)
  }

  joined <- scenarios %>%
    select("draw", "sampled_season") %>%
    inner_join(age_data, by = c("sampled_season" = "season"), relationship = "many-to-many")
  gamma_draws <- vapply(age_columns, function(column) {
    stats::rgamma(
      nrow(joined),
      shape = pmax(joined[[column]], 0) + dirichlet_offset,
      rate = 1
    )
  }, numeric(nrow(joined)))
  joined$age_share <- gamma_draws[, target_age_group] / rowSums(gamma_draws)

  joined %>%
    transmute(
      draw = .data$draw,
      season = .data$sampled_season,
      hhs_region = .data$hhs_region,
      week = .data$mmwr_week,
      age_group = target_age_group,
      age_share = .data$age_share
    )
}

draw_statewise_gam_age_decision_burden <- function(
    state_models,
    ilinet_hhs_age,
    hhs_regions,
    state_population,
    target_age_group,
    n_draws = 1000,
    seed = 20260507,
    dirichlet_offset = 0.5,
    sample_complete_seasons_only = tolower(Sys.getenv(
      "SAMPLE_COMPLETE_SEASONS_ONLY", "true"
    )) %in% c("1", "true", "yes"),
    timing_shift_sd = as.numeric(Sys.getenv("TIMING_SHIFT_SD", "0.75"))) {
  scenarios <- sample_decision_scenarios(
    state_models = state_models,
    n_draws = n_draws,
    seed = seed,
    sample_complete_seasons_only = sample_complete_seasons_only,
    timing_shift_sd = timing_shift_sd
  )

  set.seed(seed + 1L)
  all_age_draws <- bind_rows(lapply(names(state_models), function(state_name) {
    obj <- state_models[[state_name]]
    draw_latent_gam_burden(
      fit = obj$fit,
      model_data = obj$data,
      n_draws = n_draws,
      seed = sample.int(.Machine$integer.max, 1L),
      sampled_seasons = scenarios$sampled_season
    )
  })) %>%
    left_join(
      scenarios %>% select("draw", "timing_shift_weeks"),
      by = "draw",
      relationship = "many-to-one"
    )

  mapping <- state_hhs_mapping(hhs_regions, unique(all_age_draws$state))
  if (nrow(mapping) != length(unique(all_age_draws$state))) {
    missing_states <- setdiff(unique(all_age_draws$state), mapping$state)
    stop("HHS region mapping missing for: ", paste(missing_states, collapse = ", "))
  }
  age_share <- draw_hhs_age_share(
    ilinet_hhs_age = ilinet_hhs_age,
    scenarios = scenarios,
    target_age_group = target_age_group,
    seed = seed + 2L,
    dirichlet_offset = dirichlet_offset
  )

  age_draws <- all_age_draws %>%
    left_join(mapping, by = "state", relationship = "many-to-one") %>%
    inner_join(
      age_share,
      by = c("draw", "season", "hhs_region", "week"),
      relationship = "many-to-one"
    ) %>%
    mutate(
      burden = .data$burden * .data$age_share,
      age_group = target_age_group,
      posterior_ili_prob = NA_real_,
      model_source = "all_age_ilinet_gam_with_hhs_age_composition"
    ) %>%
    select(
      "draw",
      "state",
      "age_group",
      "season",
      "week",
      "burden",
      "posterior_ili_prob",
      "model_source",
      "timing_shift_weeks"
    ) %>%
    shift_burden_curves()

  add_population_weighted_us(age_draws, state_population)
}

summarise_gam_basis_check <- function(fit) {
  checks <- tryCatch(mgcv::k.check(fit), error = function(e) NULL)
  empty <- tibble(
    global_k_index = NA_real_,
    global_k_check_p_value = NA_real_,
    season_k_index = NA_real_,
    season_k_check_p_value = NA_real_
  )
  if (is.null(checks) || nrow(checks) == 0L) {
    return(empty)
  }

  terms <- rownames(checks)
  global_row <- which(terms == "s(season_pos)")
  season_row <- which(grepl("season_pos,season", gsub(" ", "", terms), fixed = TRUE))
  value <- function(rows, column, fun = min) {
    x <- suppressWarnings(as.numeric(checks[rows, column]))
    if (length(x) == 0L || all(is.na(x))) NA_real_ else fun(x, na.rm = TRUE)
  }

  tibble(
    global_k_index = value(global_row, "k-index"),
    global_k_check_p_value = value(global_row, "p-value"),
    season_k_index = value(season_row, "k-index"),
    season_k_check_p_value = value(season_row, "p-value")
  )
}

summarise_gam_residual_acf <- function(fit, model_data, max_lag = 4L) {
  residual_data <- model_data %>%
    mutate(deviance_residual = as.numeric(stats::residuals(fit, type = "deviance"))) %>%
    filter(.data$total_patients_model > 0) %>%
    group_by(.data$season) %>%
    arrange(.data$season_pos, .by_group = TRUE) %>%
    group_split()

  by_season <- bind_rows(lapply(residual_data, function(dat) {
    lag_max <- min(max_lag, nrow(dat) - 1L)
    if (lag_max < 1L || stats::sd(dat$deviance_residual, na.rm = TRUE) == 0) {
      return(tibble(lag1_acf = NA_real_, max_abs_acf_lag1_to_4 = NA_real_))
    }
    values <- stats::acf(
      dat$deviance_residual,
      lag.max = lag_max,
      plot = FALSE,
      na.action = stats::na.pass
    )$acf[-1L]
    tibble(
      lag1_acf = values[[1]],
      max_abs_acf_lag1_to_4 = max(abs(values), na.rm = TRUE)
    )
  }))

  tibble(
    median_lag1_residual_acf = median(by_season$lag1_acf, na.rm = TRUE),
    max_abs_lag1_residual_acf = max(abs(by_season$lag1_acf), na.rm = TRUE),
    share_seasons_abs_lag1_acf_gt_0_3 = mean(abs(by_season$lag1_acf) > 0.3, na.rm = TRUE),
    median_max_abs_residual_acf_lag1_to_4 = median(
      by_season$max_abs_acf_lag1_to_4,
      na.rm = TRUE
    )
  )
}

gam_diagnostic_row <- function(model_object, state_name) {
  fit <- model_object$fit
  config <- model_object$smoothing_config
  tibble(
    state = state_name,
    n = stats::nobs(fit),
    edf = sum(fit$edf),
    scale = fit$scale,
    deviance_explained = summary(fit)$dev.expl,
    aic = AIC(fit),
    gam_family = config$family,
    gam_global_k = config$global_k,
    gam_season_k = config$season_k,
    gam_gamma = config$gamma
  ) %>%
    bind_cols(summarise_gam_basis_check(fit)) %>%
    bind_cols(summarise_gam_residual_acf(fit, model_object$data))
}

save_gam_diagnostics <- function(fit, model_data, path_prefix) {
  writeLines(capture.output(summary(fit)), paste0(path_prefix, "_summary.txt"))
  diagnostics <- tibble(
    n = stats::nobs(fit),
    edf = sum(fit$edf),
    scale = fit$scale,
    deviance_explained = summary(fit)$dev.expl,
    aic = AIC(fit)
  ) %>%
    bind_cols(summarise_gam_basis_check(fit)) %>%
    bind_cols(summarise_gam_residual_acf(fit, model_data))
  readr::write_csv(diagnostics, paste0(path_prefix, "_diagnostics.csv"))
}
