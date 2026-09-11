# Helpers for ordering MMWR weeks within an influenza season.

season_week_order <- function(start_week = 23L, end_week = 22L, max_week = 52L) {
  c(seq.int(start_week, max_week), seq.int(1L, end_week))
}

week_to_season_index <- function(week, start_week = 23L, end_week = 22L, max_week = 52L) {
  week <- as.integer(week)
  max_week <- as.integer(max_week)
  if (length(max_week) == 1L) {
    max_week <- rep(max_week, length(week))
  }
  if (length(max_week) != length(week)) {
    stop("max_week must have length 1 or the same length as week.")
  }

  out <- ifelse(
    week >= start_week & week <= max_week,
    week - start_week + 1L,
    max_week - start_week + 1L + week
  )

  bad <- is.na(week) | week < 1L | week > max_week |
    (week < start_week & week > end_week)
  if (any(bad)) {
    stop("Week(s) outside modeled season order: ", paste(unique(week[bad]), collapse = ", "))
  }
  out
}

season_index_to_week <- function(index, start_week = 23L, end_week = 22L, max_week = 52L) {
  order <- season_week_order(start_week, end_week, max_week)
  if (any(index < 1L | index > length(order), na.rm = TRUE)) {
    stop("Season index outside 1:", length(order))
  }
  order[as.integer(index)]
}

weeks_since_full_response <- function(
    vaccination_week,
    burden_week,
    immune_lag_weeks = 2L,
    max_week = 52L) {
  week_to_season_index(burden_week, max_week = max_week) -
    week_to_season_index(vaccination_week, max_week = max_week) -
    immune_lag_weeks
}

is_in_burden_window <- function(week, start_week = 36L, end_week = 22L) {
  week >= start_week | week <= end_week
}

end_of_month_week <- function(month) {
  switch(
    as.character(month),
    "September" = 39L,
    "October" = 44L,
    "November" = 48L,
    stop("Supported months are September, October, and November.")
  )
}
