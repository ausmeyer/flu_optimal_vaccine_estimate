library(dplyr)
library(tibble)

prepare_national_ili_plus <- function(ilinet_national, nrevss_national) {
  ilinet_required <- c(
    "season", "season_start_year", "mmwr_year", "mmwr_week",
    "week_start", "ilitotal", "total_patients",
    "unweighted_ili", "weighted_ili"
  )
  nrevss_required <- c(
    "source_region_type", "season", "season_start_year", "mmwr_year",
    "mmwr_week", "week_start", "positive_tests", "total_tests",
    "nrevss_lab_series"
  )
  missing_ilinet <- setdiff(ilinet_required, names(ilinet_national))
  missing_nrevss <- setdiff(nrevss_required, names(nrevss_national))
  if (length(missing_ilinet) > 0L) {
    stop("National ILINet data missing: ", paste(missing_ilinet, collapse = ", "))
  }
  if (length(missing_nrevss) > 0L) {
    stop("National NREVSS data missing: ", paste(missing_nrevss, collapse = ", "))
  }

  key <- c("season", "season_start_year", "mmwr_year", "mmwr_week")
  ilinet <- ilinet_national %>%
    select(
      all_of(key),
      "week_start", "ilitotal", "total_patients",
      "unweighted_ili", "weighted_ili"
    ) %>%
    rename(
      week_start_ilinet = "week_start",
      ilinet_ili_visits = "ilitotal",
      ilinet_total_visits = "total_patients",
      ilinet_unweighted_ili_percent = "unweighted_ili",
      ilinet_weighted_ili_percent = "weighted_ili"
    )
  nrevss <- nrevss_national %>%
    filter(.data$source_region_type == "national") %>%
    select(
      all_of(key),
      "week_start", "positive_tests", "total_tests", "nrevss_lab_series"
    ) %>%
    rename(
      week_start_nrevss = "week_start",
      nrevss_positive_tests = "positive_tests",
      nrevss_total_tests = "total_tests"
    )

  if (anyDuplicated(ilinet[key])) {
    stop("National ILINet data contain duplicate season-week keys.")
  }
  if (anyDuplicated(nrevss[key])) {
    stop("National NREVSS data contain duplicate season-week keys.")
  }

  joined <- full_join(
    ilinet %>% mutate(has_ilinet = TRUE),
    nrevss %>% mutate(has_nrevss = TRUE),
    by = key,
    relationship = "one-to-one"
  )
  unmatched <- joined %>%
    filter(is.na(.data$has_ilinet) | is.na(.data$has_nrevss))
  if (nrow(unmatched) > 0L) {
    stop(
      "ILINet and NREVSS do not align for ", nrow(unmatched),
      " national season-week rows."
    )
  }
  mismatched_dates <- joined %>%
    filter(
      !is.na(.data$week_start_ilinet),
      !is.na(.data$week_start_nrevss),
      as.Date(.data$week_start_ilinet) != as.Date(.data$week_start_nrevss)
    )
  if (nrow(mismatched_dates) > 0L) {
    stop(
      "ILINet and NREVSS assign different start dates to ",
      nrow(mismatched_dates), " national season-week rows."
    )
  }

  joined %>%
    mutate(
      ilinet_unweighted_ili_prob =
        .data$ilinet_ili_visits / .data$ilinet_total_visits,
      ilinet_weighted_ili_prob = .data$ilinet_weighted_ili_percent / 100,
      nrevss_positive_prob =
        .data$nrevss_positive_tests / .data$nrevss_total_tests,
      ili_plus_weighted =
        .data$ilinet_weighted_ili_prob * .data$nrevss_positive_prob,
      ili_plus_unweighted =
        .data$ilinet_unweighted_ili_prob * .data$nrevss_positive_prob,
      week_start = dplyr::coalesce(
        as.Date(.data$week_start_ilinet),
        as.Date(.data$week_start_nrevss)
      ),
      state = "US",
      age_group = "all"
    ) %>%
    select(
      all_of(key), "week_start", "state", "age_group",
      "ilinet_ili_visits", "ilinet_total_visits",
      "ilinet_unweighted_ili_percent", "ilinet_weighted_ili_percent",
      "ilinet_unweighted_ili_prob", "ilinet_weighted_ili_prob",
      "nrevss_positive_tests", "nrevss_total_tests", "nrevss_positive_prob",
      "nrevss_lab_series", "ili_plus_weighted", "ili_plus_unweighted"
    ) %>%
    arrange(.data$season_start_year, .data$mmwr_year, .data$mmwr_week)
}

filter_ili_plus_lab_scope <- function(
    ili_plus,
    lab_scope = c("all_available", "clinical_only")) {
  lab_scope <- match.arg(lab_scope)
  out <- if (lab_scope == "clinical_only") {
    coverage <- ili_plus %>%
      summarise(
        all_weeks_clinical = all(
          .data$nrevss_lab_series == "clinical_labs"
        ),
        has_window_start = any(.data$mmwr_week == 36L),
        has_window_end = any(.data$mmwr_week == 22L),
        .by = c("season", "season_start_year")
      )
    eligible <- coverage %>%
      filter(
        .data$all_weeks_clinical,
        .data$has_window_start,
        .data$has_window_end
      ) %>%
      select("season", "season_start_year")
    ili_plus %>%
      semi_join(
        eligible,
        by = c("season", "season_start_year")
      )
  } else {
    ili_plus
  }
  if (nrow(out) == 0L) {
    stop("No ILI+ rows remain for lab scope ", lab_scope, ".")
  }
  out
}

draw_us_empirical_ili_plus_burden <- function(
    ili_plus,
    n_draws = 2000,
    ili_weighting = c("cdc_population_weighted", "national_unweighted"),
    lab_scope = c("all_available", "clinical_only"),
    age_group = "all",
    seed = 20260507,
    timing_shift_sd = 0.75,
    smooth_spar = 0.65,
    sample_complete_seasons_only = TRUE) {
  ili_weighting <- match.arg(ili_weighting)
  lab_scope <- match.arg(lab_scope)
  analysis_data <- filter_ili_plus_lab_scope(ili_plus, lab_scope)

  scenarios <- sample_empirical_scenarios(
    data = analysis_data,
    n_draws = n_draws,
    seed = seed,
    complete_only = sample_complete_seasons_only,
    timing_shift_sd = timing_shift_sd
  )

  set.seed(seed + 1L)
  draws <- scenarios %>%
    inner_join(
      analysis_data,
      by = c("sampled_season" = "season"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      sampled_unweighted_ili_prob = stats::rbeta(
        n(),
        shape1 = pmax(.data$ilinet_ili_visits, 0) + 0.5,
        shape2 = pmax(
          .data$ilinet_total_visits - .data$ilinet_ili_visits,
          0
        ) + 0.5
      ),
      ilinet_weighting_ratio = if_else(
        .data$ilinet_unweighted_ili_prob > 0,
        .data$ilinet_weighted_ili_prob /
          .data$ilinet_unweighted_ili_prob,
        1
      ),
      sampled_ili_prob = if (ili_weighting == "cdc_population_weighted") {
        pmin(
          pmax(
            .data$sampled_unweighted_ili_prob *
              .data$ilinet_weighting_ratio,
            0
          ),
          1
        )
      } else {
        .data$sampled_unweighted_ili_prob
      },
      sampled_nrevss_positive_prob = stats::rbeta(
        n(),
        shape1 = pmax(.data$nrevss_positive_tests, 0) + 0.5,
        shape2 = pmax(
          .data$nrevss_total_tests - .data$nrevss_positive_tests,
          0
        ) + 0.5
      ),
      posterior_ili_plus = .data$sampled_ili_prob *
        .data$sampled_nrevss_positive_prob,
      burden = .data$posterior_ili_plus,
      posterior_ili_prob = .data$sampled_ili_prob,
      state = "US",
      age_group = age_group,
      season = .data$sampled_season,
      week = .data$mmwr_week,
      model_source = paste(
        "empirical_ili_plus",
        ili_weighting,
        lab_scope,
        sep = "_"
      )
    ) %>%
    select(
      "draw", "state", "age_group", "season", "week", "burden",
      "posterior_ili_plus", "posterior_ili_prob",
      "sampled_nrevss_positive_prob", "model_source",
      "timing_shift_weeks"
    )

  transform_empirical_burden_draws(
    draws,
    timing_shift_sd = timing_shift_sd,
    curve_smoothing = if (is.null(smooth_spar)) "none" else "cubic_spline",
    spline_spar = if (is.null(smooth_spar)) 0.65 else smooth_spar
  )
}
