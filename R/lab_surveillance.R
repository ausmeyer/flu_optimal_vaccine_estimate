library(dplyr)
library(tidyr)
library(tibble)

source("R/calendar.R")

clean_column_names <- function(x) {
  out <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", out)
}

first_existing <- function(names, candidates) {
  hit <- candidates[candidates %in% names]
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  hit[[1]]
}

normalise_nrevss_clinical_labs <- function(x, region_label) {
  names(x) <- clean_column_names(names(x))

  total_col <- first_existing(names(x), c("total_specimens", "specimens", "total_tested", "no_of_specimens"))
  pct_col <- first_existing(names(x), c("percent_positive", "pct_positive", "percent_pos", "percent_positive_for_flu"))
  a_col <- first_existing(names(x), c("total_a", "positive_a", "a_total", "flu_a"))
  b_col <- first_existing(names(x), c("total_b", "positive_b", "b_total", "flu_b"))
  region_col <- first_existing(names(x), c("region", "statename", "state", "geography"))
  date_col <- first_existing(names(x), c("wk_date", "week_start", "week_start_date"))

  required <- c(total_col, pct_col, "year", "week")
  if (any(is.na(required))) {
    stop(
      "Could not identify required NREVSS clinical lab columns. Available columns: ",
      paste(names(x), collapse = ", ")
    )
  }

  positive_a <- rep(NA_real_, nrow(x))
  positive_b <- rep(NA_real_, nrow(x))
  positives <- rep(NA_real_, nrow(x))
  if (!is.na(a_col) && !is.na(b_col)) {
    positive_a <- suppressWarnings(as.numeric(x[[a_col]]))
    positive_b <- suppressWarnings(as.numeric(x[[b_col]]))
    positives <- positive_a + positive_b
  }
  combined_a_cols <- intersect(
    c("a_2009_h1n1", "a_h1", "a_h3", "a_subtyping_not_performed", "a_unable_to_subtype", "h3n2v", "a_h5"),
    names(x)
  )
  combined_b_cols <- intersect(c("b", "bvic", "byam"), names(x))
  combined_positive_cols <- c(combined_a_cols, combined_b_cols)
  if (length(combined_a_cols) > 0L) {
    positive_a <- dplyr::coalesce(
      positive_a,
      rowSums(
        as.data.frame(lapply(x[combined_a_cols], function(col) {
          suppressWarnings(as.numeric(col))
        })),
        na.rm = TRUE
      )
    )
  }
  if (length(combined_b_cols) > 0L) {
    positive_b <- dplyr::coalesce(
      positive_b,
      rowSums(
        as.data.frame(lapply(x[combined_b_cols], function(col) {
          suppressWarnings(as.numeric(col))
        })),
        na.rm = TRUE
      )
    )
  }
  if (length(combined_positive_cols) > 0L) {
    combined_positives <- rowSums(
      as.data.frame(lapply(x[combined_positive_cols], function(col) {
        suppressWarnings(as.numeric(col))
      })),
      na.rm = TRUE
    )
    positives <- dplyr::coalesce(positives, combined_positives)
  }
  fallback_positives <- suppressWarnings(as.numeric(x[[total_col]])) * suppressWarnings(as.numeric(x[[pct_col]])) / 100
  positives <- dplyr::coalesce(positives, fallback_positives)

  week_start <- if (!is.na(date_col)) {
    as.Date(x[[date_col]])
  } else {
    as.Date(NA)
  }

  tibble(
    source_region_type = region_label,
    source_region = if (!is.na(region_col)) as.character(x[[region_col]]) else region_label,
    mmwr_year = as.integer(x$year),
    mmwr_week = as.integer(x$week),
    week_start = week_start,
    positive_tests = pmax(positives, 0),
    positive_tests_a = pmax(positive_a, 0),
    positive_tests_b = pmax(positive_b, 0),
    total_tests = pmax(suppressWarnings(as.numeric(x[[total_col]])), 0),
    percent_positive = suppressWarnings(as.numeric(x[[pct_col]]))
  ) %>%
    filter(!is.na(.data$mmwr_year), !is.na(.data$mmwr_week), .data$total_tests > 0) %>%
    mutate(
      season_start_year = if_else(.data$mmwr_week >= 36L, .data$mmwr_year, .data$mmwr_year - 1L),
      season = paste0(.data$season_start_year, "/", substr(.data$season_start_year + 1L, 3, 4)),
      in_primary_burden_window = is_in_burden_window(.data$mmwr_week),
      state = case_when(
        .data$source_region_type == "national" ~ "National",
        .data$source_region == "District Of Columbia" ~ "District of Columbia",
        TRUE ~ .data$source_region
      ),
      ilitotal = .data$positive_tests,
      total_patients = .data$total_tests,
      unweighted_ili = .data$percent_positive,
      age_group = "all"
    )
}
