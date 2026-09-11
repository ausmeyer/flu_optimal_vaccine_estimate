library(dplyr)
library(tibble)

flu_ve_onset_to_enrollment_distribution <- function() {
  tribble(
    ~lower_day, ~upper_day, ~n_observed,
    0L, 2L, 8520L,
    3L, 4L, 10754L,
    5L, 7L, 7861L
  ) %>%
    mutate(
      source = "Balasubramani et al. 2020, doi:10.1111/irv.12741",
      source_population =
        "US Flu VE Network outpatients, 2011-2016; 45 of 27,180 missing interval",
      within_interval_assumption = "discrete uniform integer day"
    )
}

draw_outcome_date_offsets <- function(
    n_draws,
    seed,
    prior = c("none", "flu_ve_outpatient_onset_to_enrollment")) {
  prior <- match.arg(prior)
  if (prior == "none") {
    return(tibble(
      draw = seq_len(n_draws),
      outcome_date_delay_days = 0,
      outcome_date_offset_weeks = 0,
      outcome_date_source = "none"
    ))
  }

  distribution <- flu_ve_onset_to_enrollment_distribution()
  set.seed(seed)
  selected <- sample(
    seq_len(nrow(distribution)),
    size = n_draws,
    replace = TRUE,
    prob = distribution$n_observed
  )
  delay_days <- vapply(
    selected,
    function(i) sample(
      seq.int(distribution$lower_day[[i]], distribution$upper_day[[i]]),
      size = 1L
    ),
    integer(1)
  )

  tibble(
    draw = seq_len(n_draws),
    outcome_date_delay_days = delay_days,
    outcome_date_offset_weeks = -delay_days / 7,
    outcome_date_source =
      "Balasubramani_2020_onset_to_enrollment_discrete_within_category"
  )
}

apply_outcome_date_alignment <- function(
    burden_draws,
    prior = c("none", "flu_ve_outpatient_onset_to_enrollment"),
    seed = 20260531) {
  prior <- match.arg(prior)
  required <- c("draw", "state", "age_group", "season", "week", "burden")
  missing <- setdiff(required, names(burden_draws))
  if (length(missing) > 0L) {
    stop("burden_draws missing: ", paste(missing, collapse = ", "))
  }

  draw_ids <- sort(unique(burden_draws$draw))
  if (!identical(draw_ids, seq_len(length(draw_ids)))) {
    stop("Outcome-date alignment requires sequential draw identifiers.")
  }
  offsets <- draw_outcome_date_offsets(
    n_draws = length(draw_ids),
    seed = seed,
    prior = prior
  )

  burden_draws %>%
    left_join(offsets, by = "draw", relationship = "many-to-one") %>%
    group_by(.data$draw, .data$state, .data$age_group, .data$season) %>%
    group_modify(function(dat, key) {
      max_week <- if ("max_week" %in% names(dat)) {
        unique(dat$max_week)
      } else if (any(dat$week == 53L)) {
        53L
      } else {
        52L
      }
      if (length(max_week) != 1L) {
        stop("Each burden curve must have one maximum MMWR week.")
      }
      season_pos <- week_to_season_index(
        dat$week,
        start_week = 36L,
        end_week = 22L,
        max_week = max_week
      )
      ordered <- order(season_pos)
      offset <- dat$outcome_date_offset_weeks[[1]]
      dat$burden <- pmax(
        stats::approx(
          x = season_pos[ordered],
          y = dat$burden[ordered],
          xout = season_pos - offset,
          rule = 2,
          ties = mean
        )$y,
        0
      )
      dat
    }) %>%
    ungroup()
}
