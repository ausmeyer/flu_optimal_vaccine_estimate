library(dplyr)
library(readr)
library(tidyr)

source("R/analysis_config.R")
source("R/fluview_download.R")

primary_season_start_years <- primary_season_start_years()
excluded_season_start_years <- excluded_season_start_years()
primary_jurisdictions <- c(state.name, "District of Columbia")

download_years <- fluview_download_seasons(primary_season_start_years)

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

redownload_ilinet <- tolower(Sys.getenv("REDOWNLOAD_ILINET", "false")) %in% c("1", "true", "yes")
raw_paths <- c(
  state = "data/raw/ilinet_state_raw.csv",
  hhs = "data/raw/ilinet_hhs_raw.csv",
  national = "data/raw/ilinet_national_raw.csv"
)

raw_files_cover_burden_window <- function(paths, years) {
  if (!all(file.exists(paths))) {
    return(FALSE)
  }
  coverage <- vapply(names(paths), function(region) {
    path <- paths[[region]]
    weeks <- suppressMessages(readr::read_csv(
      path,
      show_col_types = FALSE,
      col_select = c("year", "week")
    ))
    fluview_raw_covers_burden_window(
      weeks, years,
      first_available_date = if (region == "state") "2010-10-03" else NULL
    )
  }, logical(1))
  all(coverage)
}

if (!redownload_ilinet && raw_files_cover_burden_window(raw_paths, primary_season_start_years)) {
  state_raw <- read_csv(raw_paths[["state"]], show_col_types = FALSE)
  hhs_raw <- read_csv(raw_paths[["hhs"]], show_col_types = FALSE)
  national_raw <- read_csv(raw_paths[["national"]], show_col_types = FALSE)
} else {
  message(
    "Downloading ILINet seasons beginning in: ",
    paste(download_years, collapse = ", ")
  )
  state_raw <- fluview_download_ilinet(region = "state", years = download_years)
  hhs_raw <- fluview_download_ilinet(region = "hhs", years = download_years)
  national_raw <- fluview_download_ilinet(region = "national", years = download_years)

  write_csv(state_raw, raw_paths[["state"]])
  write_csv(hhs_raw, raw_paths[["hhs"]])
  write_csv(national_raw, raw_paths[["national"]])
}

sum_if_any_reported <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

add_season <- function(x) {
  x %>%
    mutate(
      season_start_year = if_else(.data$week >= 36L, .data$year, .data$year - 1L),
      season = paste0(.data$season_start_year, "/", substr(.data$season_start_year + 1L, 3, 4)),
      in_primary_burden_window = .data$week >= 36L | .data$week <= 22L,
      primary_analysis_season = .data$season_start_year %in% primary_season_start_years,
      excluded_covid_affected_season = .data$season_start_year %in% excluded_season_start_years
    )
}

state_aggregated <- state_raw %>%
  add_season() %>%
  filter(.data$in_primary_burden_window, .data$primary_analysis_season) %>%
  mutate(
    state = if_else(.data$region == "New York City", "New York", .data$region)
  ) %>%
  filter(.data$state %in% primary_jurisdictions) %>%
  group_by(
    .data$season,
    .data$season_start_year,
    .data$state,
    mmwr_year = .data$year,
    mmwr_week = .data$week,
    .data$week_start
  ) %>%
  summarise(
    ilitotal = sum_if_any_reported(.data$ilitotal),
    total_patients = sum_if_any_reported(.data$total_patients),
    num_providers = sum_if_any_reported(.data$num_of_providers),
    .groups = "drop"
  )

paired_missing <- is.na(state_aggregated$ilitotal) &
  is.na(state_aggregated$total_patients)
paired_zero <- !is.na(state_aggregated$ilitotal) &
  !is.na(state_aggregated$total_patients) &
  state_aggregated$ilitotal == 0 &
  state_aggregated$total_patients == 0
invalid_pairs <- xor(
  is.na(state_aggregated$ilitotal),
  is.na(state_aggregated$total_patients)
) |
  (!paired_missing & !paired_zero & (
    state_aggregated$ilitotal < 0 |
      state_aggregated$total_patients <= 0 |
      state_aggregated$ilitotal > state_aggregated$total_patients
  ))
if (any(invalid_pairs)) {
  stop(
    "State ILINet returned mismatched missing values or invalid outpatient counts.",
    call. = FALSE
  )
}

state_processed <- state_aggregated %>%
  mutate(
    missing_outpatient_counts = .env$paired_missing | .env$paired_zero,
    ilitotal = if_else(.data$missing_outpatient_counts, NA_real_, .data$ilitotal),
    total_patients = if_else(
      .data$missing_outpatient_counts,
      NA_real_,
      .data$total_patients
    ),
    num_providers = if_else(
      .data$missing_outpatient_counts,
      NA_real_,
      .data$num_providers
    )
  ) %>%
  transmute(
    .data$season,
    .data$season_start_year,
    .data$state,
    .data$mmwr_year,
    .data$mmwr_week,
    .data$week_start,
    ilitotal = .data$ilitotal,
    total_patients = .data$total_patients,
    unweighted_ili = 100 * .data$ilitotal / .data$total_patients,
    num_providers = .data$num_providers,
    age_group = "all"
  )

hhs_age_processed <- hhs_raw %>%
  add_season() %>%
  filter(.data$in_primary_burden_window, .data$primary_analysis_season) %>%
  select(
    season,
    season_start_year,
    hhs_region = .data$region,
    mmwr_year = .data$year,
    mmwr_week = .data$week,
    week_start = .data$week_start,
    total_patients = .data$total_patients,
    ilitotal = .data$ilitotal,
    age_0_4 = .data$age_0_4,
    age_5_24 = .data$age_5_24,
    age_25_49 = .data$age_25_49,
    age_50_64 = .data$age_50_64,
    age_65 = .data$age_65
  ) %>%
  pivot_longer(
    cols = starts_with("age_"),
    names_to = "age_group_raw",
    values_to = "ili_age_count"
  ) %>%
  mutate(
    age_group = recode(
      .data$age_group_raw,
      age_0_4 = "0-4",
      age_5_24 = "5-24",
      age_25_49 = "25-49",
      age_50_64 = "50-64",
      age_65 = "65+"
    )
  ) %>%
  select(-"age_group_raw")

national_processed <- national_raw %>%
  add_season() %>%
  filter(.data$in_primary_burden_window, .data$primary_analysis_season) %>%
  transmute(
    season,
    season_start_year,
    state = "US",
    mmwr_year = .data$year,
    mmwr_week = .data$week,
    week_start = .data$week_start,
    ilitotal = .data$ilitotal,
    total_patients = .data$total_patients,
    unweighted_ili = .data$unweighted_ili,
    weighted_ili = .data$weighted_ili,
    num_providers = .data$num_of_providers,
    age_group = "all"
  )

write_csv(state_processed, "data/processed/ilinet_state_all_age_primary_seasons.csv")
write_csv(hhs_age_processed, "data/processed/ilinet_hhs_age_primary_seasons.csv")
write_csv(national_processed, "data/processed/ilinet_national_all_age_primary_seasons.csv")

state_coverage <- state_processed %>%
  group_by(.data$season, .data$season_start_year) %>%
  summarise(
    n_states = n_distinct(.data$state),
    first_week_start = min(.data$week_start, na.rm = TRUE),
    last_week_start = max(.data$week_start, na.rm = TRUE),
    observed_weeks = n_distinct(paste(.data$mmwr_year, .data$mmwr_week)),
    complete_burden_window = any(.data$mmwr_week == 36L) &&
      any(.data$mmwr_week == 22L) &&
      all(diff(sort(unique(.data$week_start))) == 7),
    .groups = "drop"
  )

dir.create("outputs/diagnostics/data", recursive = TRUE, showWarnings = FALSE)
write_csv(state_coverage, "outputs/diagnostics/data/ilinet_state_coverage_by_season.csv")

print(state_coverage)
message("Wrote ILINet raw and processed files.")
